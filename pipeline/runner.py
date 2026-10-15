"""Runs the six steps in fixed order and records each one in `runs`.

qc -> assemble -> amr_calls -> deeparg_calls -> features -> predict

Any failure stops the run, marks that step and the isolate as failed, and leaves
later steps pending. Uploaded reads and scratch files are always deleted at the end;
only contigs.fasta (and the database records) are kept.
"""
import logging
import shutil
from pathlib import Path

from config import settings
from db import mongo
from pipeline import amr_calls, assemble, deeparg_calls, features, predict, qc
from pipeline.errors import MESSAGES, TOOL_FAILED, PipelineError
from pipeline.proc import ToolTimeout

log = logging.getLogger("amr.runner")


def upload_dir(isolate_id: str) -> Path:
    return settings.WORK_DIR / "uploads" / isolate_id


def contigs_path(isolate_id: str) -> Path:
    return settings.WORK_DIR / "contigs" / f"{isolate_id}.fasta"


def _step(isolate_id: str, name: str, fn):
    mongo.insert_run(isolate_id, name)
    try:
        tool_version, db_version = fn()
    except ToolTimeout:
        mongo.finish_run(isolate_id, name, "failed", error_code=TOOL_FAILED)
        raise PipelineError(TOOL_FAILED, f"{name} timed out") from None
    except PipelineError as exc:
        _try(mongo.finish_run, isolate_id, name, "failed", error_code=exc.code)
        raise
    except Exception:
        log.exception("unexpected error in step %s of %s", name, isolate_id)
        _try(mongo.finish_run, isolate_id, name, "failed", error_code=TOOL_FAILED)
        raise PipelineError(TOOL_FAILED, f"unexpected error in {name}") from None
    mongo.finish_run(isolate_id, name, "done", tool_version=tool_version, db_version=db_version)


def _try(fn, *args, **kwargs):
    try:
        fn(*args, **kwargs)
    except Exception:
        log.warning("secondary failure while recording an error", exc_info=True)


def run_isolate(isolate_id: str, r1: Path, r2: Path, species: str) -> None:
    work = settings.WORK_DIR / "runs" / isolate_id
    ctx: dict = {}
    warnings: list[str] = list((mongo.get_isolate(isolate_id) or {}).get("warnings", []))
    try:
        mongo.set_isolate_status(isolate_id, "running")
        work.mkdir(parents=True, exist_ok=True)

        def s_qc():
            ctx["qc"] = qc.run_qc(r1, r2, work / "qc")
            return ctx["qc"]["tool_version"], None

        def s_assemble():
            res = assemble.run_assembly(ctx["qc"]["r1"], ctx["qc"]["r2"], work / "assembly")
            dest = contigs_path(isolate_id)
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(res["contigs"], dest)
            ctx["contigs"] = dest
            mongo.insert_assembly(isolate_id, res["n_contigs"], res["total_length"], res["n50"], str(dest))
            return res["tool_version"], None

        def s_amr():
            res = amr_calls.run_amrfinder(ctx["contigs"], species, work / "amr")
            ctx["amr"] = res["calls"]
            mongo.insert_gene_calls(isolate_id, "amrfinder", res["calls"])
            return res["tool_version"], res["db_version"]

        def s_deeparg():
            res = deeparg_calls.run_deeparg(ctx["contigs"], work / "deeparg")
            ctx["deeparg"] = res["calls"]
            mongo.insert_gene_calls(isolate_id, "deeparg", res["calls"])
            return res["tool_version"], res["db_version"]

        def s_features():
            ctx["features"] = features.build_features(
                ctx["amr"], ctx["deeparg"], species if species in settings.SUPPORTED_SPECIES else None)
            mongo.set_isolate_status(isolate_id, "running", features=ctx["features"])
            return "features-v1", None

        def s_predict():
            preds = predict.predict_all(ctx["features"])
            for p in preds:
                mongo.upsert_prediction(isolate_id, p["drug"], p["call"], p["probability"],
                                        p["model_version"], p["evidence"], p["drug_class"])
            missing = [d for d in settings.DRUGS if d not in {p["drug"] for p in preds}]
            for d in missing:
                warnings.append(f"No trained model is installed for {d}, so it was not predicted.")
            return preds[0]["model_version"], None

        for name, fn in (("qc", s_qc), ("assemble", s_assemble), ("amr_calls", s_amr),
                         ("deeparg_calls", s_deeparg), ("features", s_features), ("predict", s_predict)):
            _step(isolate_id, name, fn)

        mongo.set_isolate_status(isolate_id, "done", warnings=warnings)
    except PipelineError as exc:
        log.error("run %s failed: %s (%s)", isolate_id, exc.code, exc.detail)
        _try(mongo.set_isolate_status, isolate_id, "failed", error_code=exc.code, error_message=exc.message)
    except Exception:
        log.exception("run %s crashed", isolate_id)
        _try(mongo.set_isolate_status, isolate_id, "failed", error_code=TOOL_FAILED,
             error_message=MESSAGES[TOOL_FAILED])
    finally:
        shutil.rmtree(upload_dir(isolate_id), ignore_errors=True)
        shutil.rmtree(work, ignore_errors=True)
