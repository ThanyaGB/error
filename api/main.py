"""FastAPI app. Serves the single-page frontend and the /api/* contract (127.0.0.1 only, no auth)."""
import hashlib
import logging
import shutil
import uuid
from contextlib import asynccontextmanager
from pathlib import Path

from fastapi import FastAPI, File, Form, Request, UploadFile
from fastapi.exceptions import RequestValidationError
from fastapi.responses import FileResponse, JSONResponse
from fastapi.staticfiles import StaticFiles
from starlette.exceptions import HTTPException as StarletteHTTPException

from api.jobs import manager
from config import settings
from db import mongo
from db.indexes import ensure_indexes
from pipeline import deeparg_calls, qc, runner
from pipeline.errors import (DB_UNAVAILABLE, HTTP_STATUS, INPUT_EMPTY, INPUT_NOT_FASTQ, MESSAGES, NOT_FOUND,
                             TOOL_FAILED, UPLOAD_TOO_LARGE, PipelineError, payload)

log = logging.getLogger("amr.api")


@asynccontextmanager
async def lifespan(app: FastAPI):
    if mongo.is_up():
        try:
            ensure_indexes()
            manager.recover()
        except Exception:
            log.warning("startup database setup skipped", exc_info=True)
    manager.start()
    yield


app = FastAPI(title="AMR Firewall", lifespan=lifespan, docs_url=None, redoc_url=None, openapi_url=None)


# ------------------------------------------------------------ error handling
def _err(code: str, message: str | None = None) -> JSONResponse:
    return JSONResponse(payload(code, message), status_code=HTTP_STATUS[code])


@app.exception_handler(PipelineError)
async def _pipeline_error(_: Request, exc: PipelineError):
    return _err(exc.code, exc.message)


@app.exception_handler(RequestValidationError)
async def _validation_error(_: Request, exc: RequestValidationError):
    return _err(INPUT_EMPTY, "Both read files (R1 and R2) are required.")


@app.exception_handler(StarletteHTTPException)
async def _http_error(_: Request, exc: StarletteHTTPException):
    if exc.status_code == 404:
        return _err(NOT_FOUND)
    return _err(TOOL_FAILED, "The request could not be processed.")


@app.exception_handler(Exception)
async def _unhandled(_: Request, exc: Exception):
    log.exception("unhandled error")
    return _err(TOOL_FAILED)


# -------------------------------------------------------------------- health
def _tools() -> dict:
    out = {}
    for name, path in (("fastp", settings.FASTP), ("spades", settings.SPADES), ("amrfinder", settings.AMRFINDER),
                       ("deeparg", settings.DEEPARG), ("prodigal", settings.PRODIGAL)):
        out[name] = "ok" if (shutil.which(path) or Path(path).exists()) else "missing"
    out["deeparg_model"] = "ok" if deeparg_calls.deeparg_model_snapshot() else "missing"
    out["trained_models"] = sum(1 for d in settings.DRUGS if (settings.ARTIFACT_DIR / f"{d}.json").exists())
    return out


@app.get("/api/health")
def health():
    return {"mongo": "ok" if mongo.is_up() else "down", "tools": _tools(),
            "supported_species": list(settings.SUPPORTED_SPECIES)}


# ------------------------------------------------------------------- uploads
def _save_upload(up: UploadFile, dest: Path) -> int:
    size = 0
    with open(dest, "wb") as out:
        while chunk := up.file.read(1 << 20):
            size += len(chunk)
            if size > settings.MAX_UPLOAD_BYTES:
                raise PipelineError(UPLOAD_TOO_LARGE, "upload exceeded limit")
            out.write(chunk)
    return size


def _suffix(filename: str) -> str:
    low = filename.lower()
    for s in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if low.endswith(s):
            return s
    return ""


@app.post("/api/isolates", status_code=202)
def create_isolate(r1: UploadFile = File(...), r2: UploadFile = File(...), species: str = Form("")):
    if not mongo.is_up():
        raise PipelineError(DB_UNAVAILABLE, "mongo ping failed")

    for up in (r1, r2):
        if not qc.has_fastq_extension(up.filename or ""):
            raise PipelineError(INPUT_NOT_FASTQ, "bad extension")

    isolate_id = uuid.uuid4().hex[:12]
    udir = runner.upload_dir(isolate_id)
    udir.mkdir(parents=True, exist_ok=True)
    try:
        p1 = udir / f"R1{_suffix(r1.filename)}"
        p2 = udir / f"R2{_suffix(r2.filename)}"
        _save_upload(r1, p1)
        _save_upload(r2, p2)
        qc.validate_pair(p1, p2)

        key = species.strip().lower().replace(" ", "_")
        warnings = []
        if key not in settings.SUPPORTED_SPECIES:
            warnings.append(
                f"Species '{species.strip() or 'not given'}' is not in the supported list "
                f"({', '.join(settings.SUPPORTED_SPECIES)}). Results were produced anyway but have not been "
                "validated for this species.")
        digest = hashlib.sha256((r1.filename + "|" + r2.filename).encode()).hexdigest()[:16]
        mongo.insert_isolate({"_id": isolate_id, "filename_hash": digest, "species": key, "status": "queued",
                              "warnings": warnings})
        mongo.init_runs(isolate_id, settings.STEPS)
    except BaseException:
        shutil.rmtree(udir, ignore_errors=True)
        raise
    manager.submit(isolate_id, p1, p2, key)
    return JSONResponse({"isolate_id": isolate_id}, status_code=202)


# ------------------------------------------------------------------- status
def _get_or_404(isolate_id: str) -> dict:
    iso = mongo.get_isolate(isolate_id)
    if iso is None:
        raise PipelineError(NOT_FOUND)
    return iso


@app.get("/api/isolates/{isolate_id}")
def get_isolate(isolate_id: str):
    iso = _get_or_404(isolate_id)
    runs = {r["step"]: r for r in mongo.get_runs(isolate_id)}
    steps = []
    for name in settings.STEPS:
        r = runs.get(name, {})
        steps.append({"name": name, "status": r.get("status", "pending"), "started_at": r.get("started_at"),
                      "ended_at": r.get("ended_at"), "error_code": r.get("error_code")})
    failed = iso.get("status") == "failed"
    return {"isolate_id": isolate_id, "status": iso.get("status"),
            "error_code": iso.get("error_code") if failed else None,
            "error_message": (iso.get("error_message") or MESSAGES.get(iso.get("error_code"))) if failed else None,
            "steps": steps}


@app.get("/api/isolates/{isolate_id}/results")
def get_results(isolate_id: str):
    iso = _get_or_404(isolate_id)
    if iso.get("status") != "done":
        raise PipelineError(NOT_FOUND, "results not ready", message="Results for this isolate are not available yet.")
    preds = {p["drug"]: p for p in mongo.get_predictions(isolate_id)}
    ordered = [preds[d] for d in settings.DRUGS if d in preds] + [p for d, p in preds.items() if d not in settings.DRUGS]
    return {"isolate_id": isolate_id, "species": iso.get("species", ""), "warnings": iso.get("warnings", []),
            "predictions": [{"drug": p["drug"], "drug_class": p.get("drug_class"), "call": p["call"],
                             "probability": p["probability"], "model_version": p["model_version"],
                             "evidence": p.get("evidence", [])} for p in ordered]}


@app.get("/api/benchmarks")
def get_benchmarks():
    rows = mongo.get_benchmarks()
    return {"benchmarks": [{"drug": r["drug"], "model": r["model"], "sensitivity": r.get("sensitivity"),
                            "specificity": r.get("specificity"), "n_test": r.get("n_test")} for r in rows]}


# ------------------------------------------------------------------- static
@app.get("/", include_in_schema=False)
def index():
    return FileResponse(settings.WEB_DIR / "index.html", media_type="text/html")


if (settings.WEB_DIR / "assets").is_dir():
    app.mount("/assets", StaticFiles(directory=settings.WEB_DIR / "assets"), name="assets")
