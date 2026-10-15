"""Idempotent seeding (upserts only; safe to run repeatedly).

  1. creates indexes
  2. loads data/splits.csv into training_isolates (if the file exists)
  3. loads data/fallback_example.json (an export of a REAL pipeline run, made by
     scripts/fallback_run.py) as the isolate "fallback-example"

Nothing is invented: if a source file is missing, that part is skipped and reported.
"""
import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from config import settings  # noqa: E402
from db import mongo  # noqa: E402
from db.indexes import ensure_indexes  # noqa: E402

FALLBACK_JSON = settings.DATA_DIR / "fallback_example.json"
FALLBACK_ID = "fallback-example"


def seed_splits() -> int:
    if not settings.SPLITS_CSV.exists():
        print(f"skip: {settings.SPLITS_CSV} not found (no training data yet)")
        return 0
    n = 0
    with open(settings.SPLITS_CSV, newline="") as fh:
        for row in csv.DictReader(fh):
            pheno = {d: row[d] for d in settings.DRUGS if row.get(d) in ("R", "S")}
            mongo.upsert_training_isolate(row["accession"], row["species"], pheno, row["split"])
            n += 1
    print(f"training_isolates upserted: {n}")
    return n


def seed_fallback() -> bool:
    if not FALLBACK_JSON.exists():
        print(f"skip: {FALLBACK_JSON} not found -- run scripts/fallback_run.py after a real pipeline run")
        return False
    doc = json.loads(FALLBACK_JSON.read_text())
    mongo.upsert_isolate({"_id": FALLBACK_ID, "species": doc["species"], "status": "done",
                          "error_code": None, "error_message": None, "created_at": doc.get("created_at"),
                          "warnings": doc["warnings"]})
    mongo.get_db().runs.delete_many({"isolate_id": FALLBACK_ID})
    mongo.get_db().runs.insert_many([{**r, "isolate_id": FALLBACK_ID} for r in doc["runs"]])
    for p in doc["predictions"]:
        mongo.upsert_prediction(FALLBACK_ID, p["drug"], p["call"], p["probability"], p["model_version"],
                                p["evidence"], p.get("drug_class"))
    print(f"{FALLBACK_ID}: {len(doc['predictions'])} predictions, {len(doc['runs'])} runs")
    return True


def main() -> int:
    if not mongo.is_up():
        print("MongoDB is not reachable; nothing seeded.")
        return 1
    ensure_indexes()
    print("indexes ensured")
    seed_splits()
    seed_fallback()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
