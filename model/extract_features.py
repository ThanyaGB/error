"""Run AMRFinderPlus + DeepARG on the genome FASTA of every isolate listed in data/splits.csv.

Expects  <AMR_TRAIN_CACHE>/<accession>/genome.fna  and writes  calls.json  next to it
(re-runnable: isolates with an existing calls.json are skipped).

usage: python -m model.extract_features [--no-deeparg]
"""
import argparse
import json
import sys
from pathlib import Path

import pandas as pd

from config import settings
from model.train import read_splits
from pipeline import amr_calls, deeparg_calls


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-deeparg", action="store_true")
    args = ap.parse_args()
    df = read_splits()
    done = failed = 0
    for row in df.itertuples():
        d = settings.TRAIN_CACHE_DIR / row.accession
        fasta, out = d / "genome.fna", d / "calls.json"
        if out.exists():
            continue
        if not fasta.exists():
            print(f"missing genome for {row.accession}", file=sys.stderr)
            failed += 1
            continue
        try:
            amr = amr_calls.run_amrfinder(fasta, row.species, d / "amr")["calls"]
            dp = [] if args.no_deeparg else deeparg_calls.run_deeparg(fasta, d / "deeparg")["calls"]
        except Exception as exc:                       # keep going; report at the end
            print(f"{row.accession}: {type(exc).__name__} {getattr(exc, 'code', '')}", file=sys.stderr)
            failed += 1
            continue
        out.write_text(json.dumps({"species": row.species, "amr": amr, "deeparg": dp,
                                   "deeparg_used": not args.no_deeparg}))
        done += 1
    print(f"extracted {done}, failed/missing {failed}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
