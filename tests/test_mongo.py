from db import mongo


def test_insert_isolate_and_read_back():
    mongo.insert_isolate({"_id": "t-roundtrip", "species": "escherichia_coli", "filename_hash": "abc"})
    got = mongo.get_isolate("t-roundtrip")
    assert got["species"] == "escherichia_coli"
    assert got["status"] == "queued"
    mongo.set_isolate_status("t-roundtrip", "failed", error_code="TOOL_FAILED", error_message="x")
    assert mongo.get_isolate("t-roundtrip")["status"] == "failed"


def test_prediction_unique_per_isolate_and_drug():
    mongo.upsert_prediction("t-pred", "meropenem", "S", 0.1, "v", [], "carbapenem")
    mongo.upsert_prediction("t-pred", "meropenem", "R", 0.9, "v", ["blaKPC-2"], "carbapenem")
    preds = mongo.get_predictions("t-pred")
    assert len(preds) == 1 and preds[0]["call"] == "R"


def test_benchmark_upsert_is_idempotent():
    mongo.insert_benchmark("tetracycline", "custom", 0.9, 0.8, 10)
    mongo.insert_benchmark("tetracycline", "custom", 0.91, 0.8, 10)
    rows = [b for b in mongo.get_benchmarks() if b["drug"] == "tetracycline"]
    assert len(rows) == 1 and rows[0]["sensitivity"] == 0.91
