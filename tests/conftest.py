import os
import tempfile
from pathlib import Path

# Must be set before `config.settings` is imported anywhere: tests use their own database and work dir.
_TMP = Path(tempfile.mkdtemp(prefix="amr_tests_"))
os.environ["AMR_DB_NAME"] = "amr_firewall_test"
os.environ["AMR_WORK_DIR"] = str(_TMP / "work")
os.environ["AMR_ARTIFACT_DIR"] = str(_TMP / "artifacts")

import pytest  # noqa: E402

from config import settings  # noqa: E402
from db import mongo  # noqa: E402
from db.indexes import ensure_indexes  # noqa: E402


@pytest.fixture(scope="session", autouse=True)
def test_db():
    if not mongo.is_up():
        pytest.exit("MongoDB must be running on 127.0.0.1:27017 for the test-suite", returncode=2)
    assert settings.DB_NAME == "amr_firewall_test"
    mongo._client().drop_database(settings.DB_NAME)
    ensure_indexes()
    yield
    mongo._client().drop_database(settings.DB_NAME)


@pytest.fixture()
def tmp(tmp_path):
    return tmp_path


@pytest.fixture()
def client():
    from fastapi.testclient import TestClient
    from api.main import app
    with TestClient(app, raise_server_exceptions=False) as c:
        yield c


@pytest.fixture()
def example_fixture():
    """A TEST FIXTURE under the id 'fallback-example' in the *test* database only.

    It exists so API shapes can be tested. It is not real data and is never copied to amr_firewall.
    """
    iid = "fallback-example"
    mongo.upsert_isolate({"_id": iid, "species": "escherichia_coli", "status": "done", "error_code": None,
                          "error_message": None, "warnings": ["test fixture"]})
    mongo.init_runs(iid, settings.STEPS)
    for s in settings.STEPS:
        mongo.insert_run(iid, s)
        mongo.finish_run(iid, s, "done", tool_version="t", db_version="d")
    for i, (drug, cls) in enumerate(settings.DRUGS.items()):
        p = [0.95, 0.5, 0.1, 0.8, 0.2, 0.65][i]
        mongo.upsert_prediction(iid, drug, "R" if p >= 0.7 else "S" if p <= 0.3 else "uncertain", p, "test-v0",
                                ["geneA"] if p > 0.3 else [], cls)
    return iid
