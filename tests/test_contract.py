"""Every error code is produced by at least one real path, and every message is plain language."""
import gzip
import re

import pytest

from config import settings
from db import mongo
from pipeline import runner
from pipeline.errors import (ALL_CODES, ASSEMBLY_TIMEOUT, DB_UNAVAILABLE, INPUT_CORRUPT, INPUT_EMPTY,
                             INPUT_NOT_FASTQ, INPUT_TOO_FEW_READS, MESSAGES, MODEL_MISSING, NOT_FOUND,
                             TOOL_FAILED, UPLOAD_TOO_LARGE)
from tests.helpers import fake_tool, write_fastq

seen: dict[str, dict] = {}


def post(client, name1, data1, name2=None, data2=None, species="escherichia_coli"):
    r = client.post("/api/isolates", files={"r1": (name1, data1), "r2": (name2 or name1, data2 if data2 is not None else data1)},
                    data={"species": species})
    return r


def record(r, code):
    body = r.json()
    assert body["error_code"] == code, body
    assert set(body) == {"error_code", "error_message"}
    seen[code] = body
    return r


def test_input_not_fastq(client, tmp):
    record(post(client, "a.fastq", b">seq\nACGT\n"), INPUT_NOT_FASTQ)
    r = record(post(client, "reads.txt", write_fastq(tmp / "ok.fastq").read_bytes()), INPUT_NOT_FASTQ)
    assert r.status_code == 400


def test_input_empty(client):
    assert record(post(client, "a.fastq", b""), INPUT_EMPTY).status_code == 400


def test_input_corrupt(client, tmp):
    data = write_fastq(tmp / "g.fastq.gz").read_bytes()
    assert record(post(client, "a.fastq.gz", data[: len(data) // 2]), INPUT_CORRUPT).status_code == 400


def test_input_too_few_reads(client, tmp):
    assert record(post(client, "a.fastq", write_fastq(tmp / "s.fastq", n_reads=500).read_bytes()),
                  INPUT_TOO_FEW_READS).status_code == 400


def test_upload_too_large(client, tmp, monkeypatch):
    monkeypatch.setattr(settings, "MAX_UPLOAD_BYTES", 1000)
    assert record(post(client, "a.fastq", write_fastq(tmp / "b.fastq").read_bytes()), UPLOAD_TOO_LARGE).status_code == 413


def test_not_found(client):
    assert record(client.get("/api/isolates/nope"), NOT_FOUND).status_code == 404
    assert record(client.get("/api/isolates/nope/results"), NOT_FOUND).status_code == 404
    assert client.get("/api/does-not-exist").json()["error_code"] == NOT_FOUND


def test_db_unavailable_and_health_down(client, tmp, monkeypatch):
    monkeypatch.setattr(settings, "MONGO_URI", "mongodb://127.0.0.1:1")
    monkeypatch.setattr(settings, "MONGO_TIMEOUT_MS", 300)
    h = client.get("/api/health")
    assert h.status_code == 200 and h.json()["mongo"] == "down"
    assert record(post(client, "a.fastq", write_fastq(tmp / "b.fastq").read_bytes()), DB_UNAVAILABLE).status_code == 503
    assert client.get("/api/isolates/fallback-example").json()["error_code"] == DB_UNAVAILABLE   # no crash


def _queue_run(tmp, species="escherichia_coli"):
    iid = "t-" + species[:3] + str(abs(hash(str(tmp))) % 10**6)
    up = runner.upload_dir(iid)
    up.mkdir(parents=True, exist_ok=True)
    r1, r2 = write_fastq(up / "R1.fastq"), write_fastq(up / "R2.fastq")
    mongo.insert_isolate({"_id": iid, "species": species, "warnings": []})
    mongo.init_runs(iid, settings.STEPS)
    return iid, r1, r2, up


def test_tool_failed_marks_step_and_isolate_failed(tmp, monkeypatch):
    monkeypatch.setattr(settings, "FASTP", str(fake_tool(tmp / "fastp", "echo /secret/path/trace >&2; exit 1")))
    iid, r1, r2, up = _queue_run(tmp)
    runner.run_isolate(iid, r1, r2, "escherichia_coli")
    iso = mongo.get_isolate(iid)
    assert iso["status"] == "failed" and iso["error_code"] == TOOL_FAILED
    assert "/secret" not in iso["error_message"]
    runs = {r["step"]: r["status"] for r in mongo.get_runs(iid)}
    assert runs["qc"] == "failed" and all(runs[s] == "pending" for s in settings.STEPS[1:])
    assert not up.exists()                                      # uploaded reads removed
    seen[TOOL_FAILED] = {"error_message": iso["error_message"]}


def test_assembly_timeout_through_runner(tmp, monkeypatch):
    # a fastp stand-in that succeeds, then a SPAdes stand-in that hangs
    def fastp_body():
        return ('while [ $# -gt 0 ]; do case "$1" in -i) i=$2;; -I) I=$2;; -o) o=$2;; -O) O=$2;; -j) j=$2;; esac; shift; done\n'
                'cp "$i" "$o"; cp "$I" "$O"; mkdir -p "$(dirname "$j")"\n'
                'echo \'{"summary":{"before_filtering":{"total_reads":24000},"after_filtering":{"total_reads":24000}}}\' > "$j"')
    monkeypatch.setattr(settings, "FASTP", str(fake_tool(tmp / "fastp", fastp_body())))
    monkeypatch.setattr(settings, "SPADES", str(fake_tool(tmp / "spades.py", "sleep 30")))
    monkeypatch.setattr(settings, "ASSEMBLY_TIMEOUT_S", 1)
    iid, r1, r2, up = _queue_run(tmp, "klebsiella_pneumoniae")
    runner.run_isolate(iid, r1, r2, "klebsiella_pneumoniae")
    iso = mongo.get_isolate(iid)
    assert iso["status"] == "failed" and iso["error_code"] == ASSEMBLY_TIMEOUT
    runs = {r["step"]: r["status"] for r in mongo.get_runs(iid)}
    assert runs["qc"] == "done" and runs["assemble"] == "failed" and runs["amr_calls"] == "pending"
    seen[ASSEMBLY_TIMEOUT] = {"error_message": iso["error_message"]}


def test_model_missing_through_predict(tmp, monkeypatch):
    from pipeline import predict
    monkeypatch.setattr(settings, "ARTIFACT_DIR", tmp / "none")
    with pytest.raises(Exception) as ei:
        predict.predict_all({"values": {}, "sources": {}})
    assert ei.value.code == MODEL_MISSING
    seen[MODEL_MISSING] = {"error_message": ei.value.message}


def test_every_code_is_covered_and_plain_language():
    assert set(ALL_CODES) == {"INPUT_NOT_FASTQ", "INPUT_EMPTY", "INPUT_CORRUPT", "INPUT_TOO_FEW_READS",
                              "UPLOAD_TOO_LARGE", "ASSEMBLY_TIMEOUT", "TOOL_FAILED", "MODEL_MISSING",
                              "DB_UNAVAILABLE", "NOT_FOUND"}
    # all codes in the contract were produced by a test above (run this file as a whole)
    assert set(seen) == set(ALL_CODES), f"not exercised: {set(ALL_CODES) - set(seen)}"
    for code, msg in list(MESSAGES.items()) + [(c, b["error_message"]) for c, b in seen.items()]:
        assert msg and msg[0].isupper() and msg.rstrip().endswith("."), (code, msg)
        assert not re.search(r"Traceback|Exception|\.py\b|File \"|/home/|/tmp/|/mnt/|[A-Z]:\\\\|\bline \d+", msg), (code, msg)
