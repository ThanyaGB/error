"""End-to-end check against a real server. Prints PASS / FAIL for every step; exit code 1 if any FAIL.

usage:  python scripts/e2e_check.py [--url http://127.0.0.1:8000]   (default: starts its own server on :8765)

Step 3 uses a read pair from testFiles/ (*R1*/*R2* FASTQ) if that folder exists, otherwise it generates a
small synthetic pair. Synthetic random reads cannot be assembled, so there the *expected* outcome is a clean
terminal state (done, or failed with a documented error code and a plain message), never a hang or a crash.
"""
import argparse
import gzip
import json
import os
import random
import subprocess
import sys
import tempfile
import time
import urllib.error
import urllib.request
import uuid
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
CODES = {"INPUT_NOT_FASTQ", "INPUT_EMPTY", "INPUT_CORRUPT", "INPUT_TOO_FEW_READS", "UPLOAD_TOO_LARGE",
         "ASSEMBLY_TIMEOUT", "TOOL_FAILED", "MODEL_MISSING", "DB_UNAVAILABLE", "NOT_FOUND"}
DRUGS = ["ciprofloxacin", "ceftazidime", "meropenem", "gentamicin", "tetracycline", "trimethoprim-sulfamethoxazole"]
failures = 0


def report(name, ok, detail=""):
    global failures
    failures += 0 if ok else 1
    print(f"{'PASS' if ok else 'FAIL'}  {name}" + (f"  -- {detail}" if detail else ""), flush=True)


def http(method, url, body=None, headers=None):
    req = urllib.request.Request(url, data=body, method=method, headers=headers or {})
    try:
        with urllib.request.urlopen(req, timeout=300) as r:
            return r.status, r.read()
    except urllib.error.HTTPError as e:
        return e.code, e.read()


def getj(url):
    status, raw = http("GET", url)
    return status, json.loads(raw)


def multipart(files: dict, fields: dict):
    boundary = uuid.uuid4().hex
    parts = []
    for k, v in fields.items():
        parts.append(f'--{boundary}\r\nContent-Disposition: form-data; name="{k}"\r\n\r\n{v}\r\n'.encode())
    for k, (fname, data) in files.items():
        parts.append(f'--{boundary}\r\nContent-Disposition: form-data; name="{k}"; filename="{fname}"\r\n'
                     f'Content-Type: application/octet-stream\r\n\r\n'.encode() + data + b"\r\n")
    parts.append(f"--{boundary}--\r\n".encode())
    return b"".join(parts), {"Content-Type": f"multipart/form-data; boundary={boundary}"}


def upload(base, r1, r2, species="escherichia_coli"):
    body, hdr = multipart({"r1": r1, "r2": r2}, {"species": species})
    status, raw = http("POST", base + "/api/isolates", body, hdr)
    return status, json.loads(raw)


def synthetic_pair(n=12_000, length=100):
    rng = random.Random(7)
    def one(tag):
        lines = []
        for i in range(n):
            lines += [f"@s{i}/{tag}", "".join(rng.choice("ACGT") for _ in range(length)), "+", "I" * length]
        return gzip.compress(("\n".join(lines) + "\n").encode())
    return ("R1.fastq.gz", one(1)), ("R2.fastq.gz", one(2))


def find_test_pair():
    d = ROOT / "testFiles"
    if not d.is_dir():
        return None
    r1 = sorted(p for p in d.iterdir() if "R1" in p.name.upper() or "_1." in p.name)
    r2 = sorted(p for p in d.iterdir() if "R2" in p.name.upper() or "_2." in p.name)
    return ((r1[0].name, r1[0].read_bytes()), (r2[0].name, r2[0].read_bytes())) if r1 and r2 else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--url")
    ap.add_argument("--timeout", type=int, default=2700, help="seconds to wait for the pipeline")
    args = ap.parse_args()

    proc = None
    base = args.url
    if not base:
        port = 8765
        base = f"http://127.0.0.1:{port}"
        proc = subprocess.Popen([sys.executable, "-m", "uvicorn", "api.main:app", "--host", "127.0.0.1",
                                 "--port", str(port), "--log-level", "warning"], cwd=ROOT)
    try:
        # 1. server up
        up = False
        for _ in range(60):
            try:
                getj(base + "/api/health")
                up = True
                break
            except Exception:
                time.sleep(0.5)
        report("1 server starts and answers /api/health", up)
        if not up:
            return
        status, health = getj(base + "/api/health")
        report("2 health shape + mongo ok + species list",
               status == 200 and health.get("mongo") == "ok" and health.get("supported_species")
               == ["escherichia_coli", "klebsiella_pneumoniae"] and isinstance(health.get("tools"), dict),
               json.dumps({k: health[k] for k in ("mongo", "supported_species")}))
        status, page = http("GET", base + "/")
        report("2b GET / serves the page", status == 200 and b"api/isolates" in page)

        # 3. upload a small pair and follow it to a terminal state
        pair = find_test_pair()
        src = "testFiles/" if pair else "synthetic reads (no testFiles/ folder)"
        r1, r2 = pair or synthetic_pair()
        status, body = upload(base, r1, r2)
        if status != 202 or "isolate_id" not in body:
            report(f"3 upload {src}", False, f"HTTP {status} {body}")
        else:
            iid, t0, last = body["isolate_id"], time.time(), None
            while time.time() - t0 < args.timeout:
                _, st = getj(f"{base}/api/isolates/{iid}")
                last = st
                if st["status"] in ("done", "failed"):
                    break
                time.sleep(2)
            steps = " ".join(f"{s['name']}={s['status']}" for s in last["steps"])
            if last["status"] == "done":
                _, res = getj(f"{base}/api/isolates/{iid}/results")
                report(f"3 upload {src} -> done", len(res["predictions"]) > 0, f"{len(res['predictions'])} predictions; {steps}")
            elif last["status"] == "failed":
                ok = last["error_code"] in CODES and bool(last["error_message"]) and "/" not in last["error_message"]
                report(f"3 upload {src} -> clean failure (not 'done')", ok,
                       f"{last['error_code']}: {last['error_message']} | {steps}")
            else:
                report(f"3 upload {src} reaches a terminal state", False, f"still {last['status']} | {steps}")

        # 4. example
        status, res = getj(base + "/api/isolates/fallback-example/results")
        ok = status == 200 and [p["drug"] for p in res.get("predictions", [])] == DRUGS
        report("4 fallback-example results (6 drugs)", ok,
               "" if ok else f"HTTP {status} {res} -- run a real pipeline, scripts/fallback_run.py, scripts/seed_mongo.py")

        # 5. non-FASTQ
        status, body = upload(base, ("notes.fastq", b">this is a FASTA file\nACGT\n"), ("notes2.fastq", b">x\nACGT\n"))
        report("5 non-FASTQ upload -> INPUT_NOT_FASTQ", status == 400 and body.get("error_code") == "INPUT_NOT_FASTQ"
               and bool(body.get("error_message")), f"HTTP {status} {body}")

        # 6. benchmarks
        status, b = getj(base + "/api/benchmarks")
        rows = b.get("benchmarks")
        ok = status == 200 and isinstance(rows, list) and all(
            set(r) == {"drug", "model", "sensitivity", "specificity", "n_test"} for r in rows)
        report("6 /api/benchmarks shape", ok, f"{len(rows or [])} rows" + (" (empty: no benchmark run yet)" if not rows else ""))
    finally:
        if proc:
            proc.terminate()
            proc.wait(timeout=15)
    print(f"\n{'ALL STEPS PASSED' if not failures else str(failures) + ' STEP(S) FAILED'}")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
