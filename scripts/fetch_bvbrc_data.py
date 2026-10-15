"""Build data/splits.csv and the training cache from public BV-BRC data.

Everything here is public data used as-is (declared in DECLARATIONS.md):

  * phenotypes  -- BV-BRC `genome_amr` API, evidence == "Laboratory Method" only,
                   resistant_phenotype in {Resistant, Susceptible}
  * sequences   -- BV-BRC `genome_sequence` API (no FTP: ftp.bvbrc.org is unreachable here)

usage:
  python scripts/fetch_bvbrc_data.py [--out DIR] [--splits FILE] [--seed 42]
                                     [--min-drugs 4] [--max-per-species 40] [--test-frac 0.3]

Writes:
  <splits>                      columns: accession,species,split,<one column per drug with R/S>
  <out>/<accession>/genome.fna  one FASTA per isolate (headers kept from BV-BRC)

`accession` is the BV-BRC genome id (e.g. 562.7707). Nothing is invented: a drug column is empty
when BV-BRC has no laboratory measurement for that isolate, and such cells are ignored downstream.
"""
import argparse
import csv
import json
import random
import sys
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from config import settings  # noqa: E402

API = "https://www.bv-brc.org/api"
# species key used by our pipeline -> BV-BRC taxon id
TAXONS = {"escherichia_coli": 562, "klebsiella_pneumoniae": 573}
# our drug name -> the antibiotic name BV-BRC writes in `genome_amr.antibiotic`
BVBRC_DRUG_NAMES = {
    "ciprofloxacin": "ciprofloxacin",
    "ceftazidime": "ceftazidime",
    "meropenem": "meropenem",
    "gentamicin": "gentamicin",
    "tetracycline": "tetracycline",
    "trimethoprim-sulfamethoxazole": "trimethoprim/sulfamethoxazole",
}
# The BV-BRC data API has no offset/start parameter, so one large request per taxon is used.
PAGE = 50_000


def api(endpoint: str, params: list[tuple[str, str]]) -> list[dict]:
    url = f"{API}/{endpoint}/?" + "&".join(urllib.parse.quote(p, safe="(),") for p in params)
    # BV-BRC rejects the default Python-urllib agent with 403, so send a real one.
    req = urllib.request.Request(url, headers={"Accept": "application/json",
                                               "User-Agent": "amr-firewall/1.0 (research build)"})
    with urllib.request.urlopen(req, timeout=180) as resp:
        return json.loads(resp.read())


def fetch_phenotypes(taxon_id: int) -> dict[str, dict[str, str]]:
    """genome_id -> {drug: 'R'|'S'} using laboratory evidence only."""
    by_antibiotic = defaultdict(set)      # bvbrc antibiotic name -> (genome id, R|S)
    rows = api("genome_amr", [
        f"eq(taxon_id,{taxon_id})", "eq(evidence,Laboratory Method)",
        "select(genome_id,antibiotic,resistant_phenotype,evidence)",
        f"limit({PAGE})",
    ])
    for r in rows:
        if r.get("resistant_phenotype") in ("Resistant", "Susceptible"):
            by_antibiotic[r["antibiotic"]].add((r["genome_id"], r["resistant_phenotype"][0]))
    print(f"  genome_amr taxon {taxon_id}: {len(rows)} laboratory rows", file=sys.stderr)

    per_genome: dict[str, dict[str, str]] = defaultdict(dict)
    for drug, bvbrc_name in BVBRC_DRUG_NAMES.items():
        for genome_id, call in by_antibiotic.get(bvbrc_name, ()):
            per_genome[genome_id][drug] = call
    return dict(per_genome)


def fetch_genome_fasta(genome_id: str) -> str:
    """Whole genome as FASTA text (all contigs, in BV-BRC order)."""
    params = [f"eq(genome_id,{genome_id})", "select(sequence_id,sequence,length)", "limit(20000)"]
    rows = api("genome_sequence", params)
    out = []
    for r in rows:
        seq = r.get("sequence") or ""
        if not seq:
            continue
        out.append(f">{r.get('sequence_id') or genome_id}\n")
        for i in range(0, len(seq), 70):
            out.append(seq[i:i + 70] + "\n")
    return "".join(out)


def choose(candidates: dict[str, dict[str, str]], min_drugs: int, cap: int, exclude: set[str]):
    picked = [(g, p) for g, p in candidates.items() if len(p) >= min_drugs and g not in exclude]
    # deterministic: most-labelled first, then by id
    picked.sort(key=lambda x: (-len(x[1]), x[0]))
    return picked[:cap]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(settings.TRAIN_CACHE_DIR))
    ap.add_argument("--splits", default=str(settings.SPLITS_CSV))
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--min-drugs", type=int, default=4)
    ap.add_argument("--max-per-species", type=int, default=40)
    ap.add_argument("--test-frac", type=float, default=0.3)
    args = ap.parse_args()

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    splits_path = Path(args.splits)
    splits_path.parent.mkdir(parents=True, exist_ok=True)

    rows, used = [], set()
    for species, taxon in TAXONS.items():
        print(f"{species} (taxon {taxon}): fetching laboratory phenotypes", file=sys.stderr)
        pheno = fetch_phenotypes(taxon)
        picked = choose(pheno, args.min_drugs, args.max_per_species, used)
        print(f"  {len(pheno)} genomes with any laboratory call, {len(picked)} kept", file=sys.stderr)
        rng = random.Random(args.seed)
        rng.shuffle(picked)
        n_test = int(round(len(picked) * args.test_frac))
        for i, (genome_id, calls) in enumerate(picked):
            used.add(genome_id)
            rows.append({"accession": genome_id, "species": species,
                         "split": "test" if i < n_test else "train", **calls})

    if not rows:
        print("no isolates selected; nothing written", file=sys.stderr)
        return 1

    fields = ["accession", "species", "split", *settings.DRUGS]
    with open(splits_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        for r in rows:
            writer.writerow({f: r.get(f, "") for f in fields})
    print(f"wrote {splits_path} ({len(rows)} isolates)")

    missing = 0
    for i, r in enumerate(rows, 1):
        dest = out_dir / r["accession"] / "genome.fna"
        if dest.exists() and dest.stat().st_size > 0:
            continue
        try:
            fasta = fetch_genome_fasta(r["accession"])
        except Exception as exc:
            print(f"{r['accession']}: download failed ({type(exc).__name__})", file=sys.stderr)
            missing += 1
            continue
        if not fasta:
            print(f"{r['accession']}: empty sequence", file=sys.stderr)
            missing += 1
            continue
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_text(fasta)
        print(f"[{i}/{len(rows)}] {r['accession']} {len(fasta) // 1000} kb")
    print(f"done: {len(rows) - missing} genomes cached, {missing} failed")

    labelled = defaultdict(lambda: defaultdict(int))
    for r in rows:
        for drug in settings.DRUGS:
            if r.get(drug) in ("R", "S"):
                labelled[r["split"]][f"{drug}:{r[drug]}"] += 1
    for split in ("train", "test"):
        print(f"{split}: " + ", ".join(f"{k}={v}" for k, v in sorted(labelled[split].items())))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
