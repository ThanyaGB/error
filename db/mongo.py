"""All MongoDB access lives here. Nothing else in the project imports pymongo.

Any connectivity problem is turned into PipelineError(DB_UNAVAILABLE) so callers
never see driver exceptions (or connection strings) leak out.
"""
import functools
import threading
from datetime import datetime, timezone

from pymongo import MongoClient, ReturnDocument
from pymongo.errors import AutoReconnect, ConnectionFailure, ServerSelectionTimeoutError

from config import settings
from pipeline.errors import DB_UNAVAILABLE, PipelineError

_clients: dict[str, MongoClient] = {}
_lock = threading.Lock()
_CONN_ERRORS = (ServerSelectionTimeoutError, ConnectionFailure, AutoReconnect)


def now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _client() -> MongoClient:
    uri = settings.MONGO_URI
    with _lock:
        if uri not in _clients:
            _clients[uri] = MongoClient(
                uri,
                serverSelectionTimeoutMS=settings.MONGO_TIMEOUT_MS,
                connectTimeoutMS=settings.MONGO_TIMEOUT_MS,
                socketTimeoutMS=15000,
                tz_aware=True,
            )
        return _clients[uri]


def get_db():
    return _client()[settings.DB_NAME]


def guarded(fn):
    @functools.wraps(fn)
    def wrapper(*args, **kwargs):
        try:
            return fn(*args, **kwargs)
        except _CONN_ERRORS as exc:
            raise PipelineError(DB_UNAVAILABLE, type(exc).__name__) from None
    return wrapper


def is_up() -> bool:
    try:
        _client().admin.command("ping")
        return True
    except Exception:
        return False


def _clean(doc):
    if doc is not None:
        doc.pop("_id_internal", None)
    return doc


# ------------------------------------------------------------------ isolates
@guarded
def insert_isolate(doc: dict) -> str:
    doc = dict(doc)
    doc.setdefault("status", "queued")
    doc.setdefault("error_code", None)
    doc.setdefault("error_message", None)
    doc.setdefault("created_at", now())
    get_db().isolates.insert_one(doc)
    return doc["_id"]


@guarded
def upsert_isolate(doc: dict) -> None:
    doc = dict(doc)
    _id = doc.pop("_id")
    get_db().isolates.update_one({"_id": _id}, {"$set": doc}, upsert=True)


@guarded
def set_isolate_status(isolate_id: str, status: str, error_code=None, error_message=None, **extra) -> None:
    fields = {"status": status, "error_code": error_code, "error_message": error_message,
              "updated_at": now(), **extra}
    get_db().isolates.update_one({"_id": isolate_id}, {"$set": fields})


@guarded
def get_isolate(isolate_id: str):
    return get_db().isolates.find_one({"_id": isolate_id})


@guarded
def find_isolates_by_status(status: str) -> list[dict]:
    return list(get_db().isolates.find({"status": status}).sort("created_at", 1))


# ---------------------------------------------------------------------- runs
@guarded
def init_runs(isolate_id: str, steps: list[str]) -> None:
    db = get_db()
    db.runs.delete_many({"isolate_id": isolate_id})
    db.runs.insert_many([
        {"isolate_id": isolate_id, "step": s, "order": i, "status": "pending",
         "started_at": None, "ended_at": None, "tool_version": None,
         "db_version": None, "error_code": None}
        for i, s in enumerate(steps)
    ])


@guarded
def insert_run(isolate_id: str, step: str) -> None:
    """Mark a step as running (creates the run document if init_runs was skipped)."""
    get_db().runs.update_one(
        {"isolate_id": isolate_id, "step": step},
        {"$set": {"status": "running", "started_at": now(), "ended_at": None, "error_code": None}},
        upsert=True,
    )


@guarded
def finish_run(isolate_id: str, step: str, status: str, tool_version=None, db_version=None,
               error_code=None) -> None:
    get_db().runs.update_one(
        {"isolate_id": isolate_id, "step": step},
        {"$set": {"status": status, "ended_at": now(), "tool_version": tool_version,
                  "db_version": db_version, "error_code": error_code}},
        upsert=True,
    )


@guarded
def get_runs(isolate_id: str) -> list[dict]:
    return list(get_db().runs.find({"isolate_id": isolate_id}).sort("order", 1))


# ---------------------------------------------------------------- assemblies
@guarded
def insert_assembly(isolate_id: str, n_contigs: int, total_length: int, n50: int, contigs_path: str) -> None:
    get_db().assemblies.update_one(
        {"isolate_id": isolate_id},
        {"$set": {"isolate_id": isolate_id, "n_contigs": n_contigs, "total_length": total_length,
                  "n50": n50, "contigs_path": contigs_path, "created_at": now()}},
        upsert=True,
    )


@guarded
def get_assembly(isolate_id: str):
    return get_db().assemblies.find_one({"isolate_id": isolate_id})


# ---------------------------------------------------------------- gene calls
@guarded
def insert_gene_calls(isolate_id: str, source: str, calls: list[dict]) -> int:
    db = get_db()
    db.gene_calls.delete_many({"isolate_id": isolate_id, "source": source})
    if calls:
        db.gene_calls.insert_many([{**c, "isolate_id": isolate_id, "source": source} for c in calls])
    return len(calls)


@guarded
def get_gene_calls(isolate_id: str, source: str | None = None) -> list[dict]:
    q = {"isolate_id": isolate_id}
    if source:
        q["source"] = source
    return list(get_db().gene_calls.find(q, {"_id": 0}))


# --------------------------------------------------------------- predictions
@guarded
def upsert_prediction(isolate_id: str, drug: str, call: str, probability: float, model_version: str,
                      evidence: list[str], drug_class: str | None = None) -> None:
    get_db().predictions.update_one(
        {"isolate_id": isolate_id, "drug": drug},
        {"$set": {"isolate_id": isolate_id, "drug": drug, "drug_class": drug_class, "call": call,
                  "probability": float(probability), "model_version": model_version,
                  "evidence": list(evidence), "created_at": now()}},
        upsert=True,
    )


@guarded
def get_predictions(isolate_id: str) -> list[dict]:
    return list(get_db().predictions.find({"isolate_id": isolate_id}, {"_id": 0}))


# ---------------------------------------------------------------- benchmarks
@guarded
def insert_benchmark(drug: str, model: str, sensitivity, specificity, n_test: int, **extra) -> None:
    get_db().benchmarks.update_one(
        {"drug": drug, "model": model},
        {"$set": {"drug": drug, "model": model, "sensitivity": sensitivity,
                  "specificity": specificity, "n_test": int(n_test), "created_at": now(), **extra}},
        upsert=True,
    )


@guarded
def get_benchmarks() -> list[dict]:
    return list(get_db().benchmarks.find({}, {"_id": 0}).sort([("drug", 1), ("model", 1)]))


# -------------------------------------------------------------------- models
@guarded
def upsert_model(drug: str, version: str, features: list[str], trained_on_n: int, artifact_path: str,
                 **extra) -> None:
    get_db().models.update_one(
        {"drug": drug, "version": version},
        {"$set": {"drug": drug, "version": version, "features": features, "trained_on_n": trained_on_n,
                  "artifact_path": artifact_path, "created_at": now(), **extra}},
        upsert=True,
    )


@guarded
def get_latest_model(drug: str):
    return get_db().models.find_one({"drug": drug}, {"_id": 0}, sort=[("created_at", -1)])


# --------------------------------------------------------- training isolates
@guarded
def upsert_training_isolate(accession: str, species: str, phenotypes: dict, split: str, **extra) -> None:
    get_db().training_isolates.update_one(
        {"accession": accession},
        {"$set": {"accession": accession, "species": species, "phenotypes": phenotypes,
                  "split": split, **extra}},
        upsert=True,
    )


@guarded
def count_training_isolates(split: str | None = None) -> int:
    return get_db().training_isolates.count_documents({"split": split} if split else {})
