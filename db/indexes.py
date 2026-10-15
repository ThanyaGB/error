"""Index definitions. ensure_indexes() is idempotent and called from seed_mongo.py and API startup."""
from pymongo import ASCENDING

from db.mongo import get_db, guarded


@guarded
def ensure_indexes() -> None:
    db = get_db()
    db.isolates.create_index([("status", ASCENDING)])
    db.runs.create_index([("isolate_id", ASCENDING)])
    db.assemblies.create_index([("isolate_id", ASCENDING)], unique=True)
    db.gene_calls.create_index([("isolate_id", ASCENDING)])
    db.gene_calls.create_index([("isolate_id", ASCENDING), ("gene", ASCENDING)])
    db.predictions.create_index([("isolate_id", ASCENDING), ("drug", ASCENDING)], unique=True)
    db.training_isolates.create_index([("accession", ASCENDING)], unique=True)
    db.training_isolates.create_index([("split", ASCENDING)])
    db.benchmarks.create_index([("drug", ASCENDING)])
    db.models.create_index([("drug", ASCENDING)])
