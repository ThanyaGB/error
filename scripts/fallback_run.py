"""Export a finished (done) real pipeline run as data/fallback_example.json for scripts/seed_mongo.py.

usage: python scripts/fallback_run.py <isolate_id>

It refuses to export anything that is not a completed run with predictions, and it tags the
warnings so the UI says the example is precomputed.
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from config import settings  # noqa: E402
from db import mongo  # noqa: E402

NOTE = "This is a precomputed example from a real earlier run, not a live analysis."


def main(isolate_id: str) -> int:
    iso = mongo.get_isolate(isolate_id)
    preds = mongo.get_predictions(isolate_id)
    runs = [{k: v for k, v in r.items() if k not in ("_id", "isolate_id")} for r in mongo.get_runs(isolate_id)]
    if not iso or iso.get("status") != "done" or not preds or [r["status"] for r in runs] != ["done"] * 6:
        print("refusing: isolate is not a completed run with 6 done steps and predictions")
        return 1
    out = {"species": iso["species"], "created_at": iso.get("created_at"),
           "warnings": [NOTE] + [w for w in iso.get("warnings", []) if w != NOTE],
           "runs": runs, "predictions": preds, "source_isolate_id": isolate_id}
    dest = settings.DATA_DIR / "fallback_example.json"
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps(out, indent=2, default=str))
    print(f"wrote {dest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1]) if len(sys.argv) > 1 else 2)
