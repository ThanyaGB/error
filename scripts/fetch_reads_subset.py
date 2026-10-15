"""Fetch the first N read pairs of a public ENA run as two local `.fastq.gz` files.

Used to produce the downsampled isolate for the live demo without downloading a whole run (a whole
ENA/SRA FASTQ file is hundreds of megabytes). Public data, used as-is:
https://www.ebi.ac.uk/ena/portal/api  (see DECLARATIONS.md)

usage:
  python scripts/fetch_reads_subset.py --run DRR002018 --pairs 500000 --out DIR
or, with the FASTQ URLs already known (skips the metadata call):
  python scripts/fetch_reads_subset.py --fastq1 https://.../X_1.fastq.gz \
                                       --fastq2 https://.../X_2.fastq.gz --pairs 500000 --out DIR

How the subset is taken: the run's FASTQ.gz files support HTTP byte ranges, so the script streams a
prefix of each compressed file, decompresses it on the fly and stops once --pairs complete records
have been written. That is the first N read pairs of the run, not a random sample; the reads
themselves are untouched.
"""
import argparse
import gzip
import json
import sys
import time
import urllib.request
import zlib
from pathlib import Path

FILEREPORT = ("https://www.ebi.ac.uk/ena/portal/api/filereport"
              "?accession={run}&result=read_run"
              "&fields=run_accession,fastq_ftp,fastq_bytes,read_count&format=json")
HEADERS = {"User-Agent": "amr-firewall/1.0 (research build)"}


def _get(url: str, *, timeout: int = 120, attempts: int = 3, **headers):
    """GET with retries: the ENA portal API answers 500 / times out intermittently."""
    last = None
    for attempt in range(1, attempts + 1):
        req = urllib.request.Request(url, headers={**HEADERS, **headers})
        try:
            return urllib.request.urlopen(req, timeout=timeout)
        except Exception as exc:                     # noqa: BLE001 - retried below
            last = exc
            if attempt < attempts:
                time.sleep(2 * attempt)
    raise SystemExit(f"GET {url} failed after {attempts} attempts: {last}")


def filereport(run: str) -> dict:
    # No "Accept: application/json": the ENA portal API answers that with HTTP 500.
    with _get(FILEREPORT.format(run=run)) as resp:
        rows = json.loads(resp.read())
    if not rows:
        raise SystemExit(f"{run}: no ENA record (is the accession a run accession?)")
    return rows[0]


RANGE_CHUNK = 4 << 20          # bytes per ranged request; a fresh connection keeps the rate up


def stream_pairs(url: str, out_gz: Path, pairs: int) -> int:
    """Write the first `pairs` complete FASTQ records of a (possibly truncated gzip) stream.

    The file is read as a series of byte ranges, one small request each: a single long-lived
    connection to the archive is throttled to a crawl, while fresh ranged requests run at full
    speed. All the ranges form one continuous byte stream for the decompressor.
    """
    out_gz.parent.mkdir(parents=True, exist_ok=True)
    kept = 0
    buf = b""                      # decompressed bytes not yet written out
    dec = zlib.decompressobj(31)
    offset = 0
    with gzip.open(out_gz, "wb") as out:
        while kept < pairs:
            end = offset + RANGE_CHUNK - 1
            with _get(url, timeout=180, **{"Range": f"bytes={offset}-{end}"}) as resp:
                raw = resp.read()
            if not raw:
                break
            offset += len(raw)
            # Concatenate, then split on newlines with a limit: a line split across two network
            # chunks must stay one line, so the pieces are never split separately.
            buf += dec.decompress(raw)
            while kept < pairs:
                parts = buf.split(b"\n", 4)
                if len(parts) < 5:
                    break
                record, buf = parts[:4], parts[4]
                if not record[0].startswith(b"@"):
                    raise SystemExit("unexpected FASTQ structure while subsetting")
                out.write(b"\n".join(record) + b"\n")
                kept += 1
            if len(raw) < RANGE_CHUNK:
                break                              # end of file
            print(f"    ... {offset / 1e6:.0f} MB read, {kept} reads", file=sys.stderr, flush=True)
    return kept


def first_read_length(out_gz: Path) -> int:
    with gzip.open(out_gz, "rt") as fh:
        fh.readline()
        return len(fh.readline().rstrip("\n"))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", help="ENA run accession, e.g. DRR002018")
    ap.add_argument("--fastq1", help="FASTQ.gz URL for R1 (skips the metadata lookup)")
    ap.add_argument("--fastq2", help="FASTQ.gz URL for R2 (skips the metadata lookup)")
    ap.add_argument("--pairs", type=int, default=500_000)
    ap.add_argument("--out", default=".")
    args = ap.parse_args()

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.fastq1 and args.fastq2:
        urls, stem = [args.fastq1, args.fastq2], args.run or Path(args.fastq1).name.split("_")[0]
    else:
        if not args.run:
            raise SystemExit("give --run, or both --fastq1 and --fastq2")
        rec = filereport(args.run)
        urls = [u for u in (rec.get("fastq_ftp") or "").split(";") if u]
        if len(urls) < 2:
            raise SystemExit(f"{args.run}: expected a paired run, got {len(urls)} FASTQ file(s)")
        urls = ["https://" + u for u in urls[:2]]
        stem = args.run
        print(f"{args.run}: {rec.get('read_count')} reads, {rec.get('fastq_bytes')} bytes total",
              file=sys.stderr)

    for side, url in zip(("R1", "R2"), urls):
        dest = out_dir / f"{stem}_{side}.fastq.gz"
        print(f"  {side}: {url}", file=sys.stderr)
        kept = stream_pairs(url, dest, args.pairs)
        print(f"  {side}: wrote {dest} ({dest.stat().st_size / 1e6:.1f} MB, {kept} reads)", file=sys.stderr)

    length = first_read_length(out_dir / f"{stem}_R1.fastq.gz")
    print(f"done: {out_dir}  (first read length {length} bases)")
    if not 50 <= length <= 300:
        print(f"warning: read length {length} is outside the pipeline's 50-300 base window", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
