"""Each endpoint returns exactly the shape in the contract (uses a fixture under 'fallback-example' in the test DB)."""
import re

from config import settings
from db import mongo

STEPS = ["qc", "assemble", "amr_calls", "deeparg_calls", "features", "predict"]


def test_root_serves_frontend(client):
    r = client.get("/")
    assert r.status_code == 200 and "text/html" in r.headers["content-type"]
    assert "api/isolates" in r.text


def test_assets_served(client):
    r = client.get("/assets/blood-slide-prep.webp")
    assert r.status_code == 200 and r.headers["content-type"] == "image/webp"


def test_health_shape(client):
    r = client.get("/api/health")
    assert r.status_code == 200
    body = r.json()
    assert set(body) == {"mongo", "tools", "supported_species"}
    assert body["mongo"] == "ok"
    assert isinstance(body["tools"], dict)
    assert body["supported_species"] == ["escherichia_coli", "klebsiella_pneumoniae"]
    assert all(re.fullmatch(r"[a-z]+(_[a-z]+)*", s) for s in body["supported_species"])


def test_isolate_status_shape(client, example_fixture):
    r = client.get("/api/isolates/fallback-example")
    assert r.status_code == 200
    b = r.json()
    assert set(b) == {"isolate_id", "status", "error_code", "error_message", "steps"}
    assert b["isolate_id"] == "fallback-example" and b["status"] == "done"
    assert b["error_code"] is None and b["error_message"] is None
    assert [s["name"] for s in b["steps"]] == STEPS
    for s in b["steps"]:
        assert set(s) == {"name", "status", "started_at", "ended_at", "error_code"}
        assert s["status"] in {"pending", "running", "done", "failed"}


def test_results_shape(client, example_fixture):
    r = client.get("/api/isolates/fallback-example/results")
    assert r.status_code == 200
    b = r.json()
    assert set(b) == {"isolate_id", "species", "warnings", "predictions"}
    assert b["species"] == "escherichia_coli" and isinstance(b["warnings"], list)
    assert [p["drug"] for p in b["predictions"]] == list(settings.DRUGS)         # 6 drugs
    for p in b["predictions"]:
        assert set(p) == {"drug", "drug_class", "call", "probability", "model_version", "evidence"}
        assert p["call"] in {"R", "S", "uncertain"}
        assert isinstance(p["probability"], float) and 0.0 <= p["probability"] <= 1.0
        assert isinstance(p["evidence"], list) and all(isinstance(g, str) for g in p["evidence"])
        expect = "R" if p["probability"] >= 0.7 else "S" if p["probability"] <= 0.3 else "uncertain"
        assert p["call"] == expect


def test_benchmarks_shape_and_empty(client):
    mongo.get_db().benchmarks.delete_many({})
    assert client.get("/api/benchmarks").json() == {"benchmarks": []}
    mongo.insert_benchmark("ciprofloxacin", "custom", 0.9, 0.8, 40, threshold=0.5)
    mongo.insert_benchmark("ciprofloxacin", "amrfinder", 0.85, 0.9, 40)
    rows = client.get("/api/benchmarks").json()["benchmarks"]
    assert len(rows) == 2
    for row in rows:
        assert set(row) == {"drug", "model", "sensitivity", "specificity", "n_test"}
        assert row["model"] in {"custom", "amrfinder"} and isinstance(row["n_test"], int)
    mongo.get_db().benchmarks.delete_many({})


def test_failed_isolate_carries_code_and_message(client):
    mongo.upsert_isolate({"_id": "t-failed", "species": "escherichia_coli", "status": "failed",
                          "error_code": "ASSEMBLY_TIMEOUT", "error_message": None})
    b = client.get("/api/isolates/t-failed").json()
    assert b["status"] == "failed" and b["error_code"] == "ASSEMBLY_TIMEOUT" and b["error_message"]
    assert [s["status"] for s in b["steps"]] == ["pending"] * 6


def test_upload_accepted_returns_202_and_id(client, tmp, monkeypatch):
    from api import main
    from tests.helpers import write_fastq
    monkeypatch.setattr(main.manager, "submit", lambda *a, **k: None)      # don't start the pipeline here
    a, b = write_fastq(tmp / "a.fastq.gz"), write_fastq(tmp / "b.fastq.gz")
    r = client.post("/api/isolates", files={"r1": ("a.fastq.gz", a.read_bytes()), "r2": ("b.fastq.gz", b.read_bytes())},
                    data={"species": "escherichia_coli"})
    assert r.status_code == 202
    iid = r.json()["isolate_id"]
    assert set(r.json()) == {"isolate_id"}
    st = client.get(f"/api/isolates/{iid}").json()
    assert st["status"] == "queued" and [s["status"] for s in st["steps"]] == ["pending"] * 6


def test_unsupported_species_runs_with_warning(client, tmp, monkeypatch):
    from api import main
    from tests.helpers import write_fastq
    monkeypatch.setattr(main.manager, "submit", lambda *a, **k: None)
    a, b = write_fastq(tmp / "a.fastq"), write_fastq(tmp / "b.fastq")
    r = client.post("/api/isolates", files={"r1": ("a.fastq", a.read_bytes()), "r2": ("b.fastq", b.read_bytes())},
                    data={"species": "salmonella_enterica"})
    assert r.status_code == 202
    iso = mongo.get_isolate(r.json()["isolate_id"])
    assert iso["warnings"] and "not in the supported list" in iso["warnings"][0]
