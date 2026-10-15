# AMR Firewall — backend

Local, offline pipeline that takes one bacterial isolate's paired Illumina reads and returns
per-drug resistance calls with the gene evidence behind them. FastAPI + MongoDB + one XGBoost
classifier per drug, with SPAdes, AMRFinderPlus and DeepARG called as external command-line tools.

> **Scoping statement (team-owned, not written by the build agent).** Two lines, locked once
> submitted. Placeholder until the team pastes the final wording here:
> `SCOPE_LINE_1` · `SCOPE_LINE_2`

---

## 1. Start it

Everything runs **inside WSL (Ubuntu)**. One command starts MongoDB if needed, seeds the database,
and serves the app:

```bash
cd ~/proj                 # /mnt/c/Users/Thanya/Desktop/error-chaos
bash scripts/serve.sh
```

Then open **http://localhost:8000/** in a Windows browser (WSL port 8000 is forwarded to Windows).
The server binds `127.0.0.1:8000` only and has no authentication — local demo use only.

Other entry points:

```bash
python -m pytest -q                 # tests (uses the separate DB amr_firewall_test)
python scripts/e2e_check.py --url http://127.0.0.1:8000   # end-to-end check against a running server
python scripts/seed_mongo.py        # re-seed (idempotent upserts)
```

Environment: micromamba env `amr` (`export PATH=$HOME/micromamba/envs/amr/bin:$PATH`),
Python 3.11.17. `requirements.txt` pins every Python package; the external programs are not pip
packages (see `DECLARATIONS.md`).

---

## 2. What the pipeline does

`POST /api/isolates` with R1, R2 and a species key. The run goes through six steps, in this order,
one assembly at a time:

| # | Step | Tool | What is stored |
|---|---|---|---|
| 1 | `qc` | fastp 1.4.0 | trimmed pair, read counts before/after |
| 2 | `assemble` | SPAdes 4.3.0 (`--isolate -t 4 -m 6`) | `assemblies` (n_contigs, total_length, N50), contigs kept at `~/amr_work/contigs/<id>.fasta` |
| 3 | `amr_calls` | AMRFinderPlus 4.2.7 (DB 2026-08-07.1) | `gene_calls` with source `amrfinder` |
| 4 | `deeparg_calls` | prodigal → DeepARG 2.0.0 (`--model LS`) | `gene_calls` with source `deeparg` |
| 5 | `features` | our code (`pipeline/features.py`) | feature vector on the isolate document |
| 6 | `predict` | our XGBoost models, one per drug | `predictions` (call, probability, evidence genes) |

Uploaded reads and all scratch files are deleted when the run ends; only the contigs and the
database records survive.

### DeepARG input decision (recorded, as required)

DeepARG's `predict` does **not** take assembled contigs: it expects reads (model `SS`) or
gene-length sequences (model `LS`). This build therefore runs

```
contigs.fasta --prodigal--> orfs.fna (nucleotide ORFs) --deeparg predict --model LS --type nucl--> *.mapping.ARG
```

and reads only the high-confidence `*.mapping.ARG` output (probability ≥ 0.8). DeepARG runs on CPU
(`--hf-device cpu`) from the local Hugging Face cache; our child-process environment sets
`HF_HUB_OFFLINE=1` and `TRANSFORMERS_OFFLINE=1`, so the demo needs no network.

---

## 3. API contract

| Method | Path | Returns |
|---|---|---|
| GET | `/` | `web/index.html` (assets under `/assets/`) |
| GET | `/api/health` | `{"mongo": "ok"\|"down", "tools": {...}, "supported_species": [...]}` |
| POST | `/api/isolates` | `202 {"isolate_id": "..."}`, or `400 {"error_code", "error_message"}` |
| GET | `/api/isolates/{id}` | `{"isolate_id", "status", "error_code", "error_message", "steps": [...]}`; status is `queued\|running\|done\|failed`; steps in the fixed order above |
| GET | `/api/isolates/{id}/results` | `{"isolate_id", "species", "warnings", "predictions": [{"drug", "drug_class", "call", "probability", "model_version", "evidence"}]}` |
| GET | `/api/benchmarks` | `{"benchmarks": [{"drug", "model", "sensitivity", "specificity", "n_test"}]}` |

Supported species: `escherichia_coli` (AMRFinderPlus `--organism Escherichia`),
`klebsiella_pneumoniae` (`Klebsiella_pneumoniae`). Any other species is accepted: the pipeline runs
and a warning is added to `warnings` instead of refusing the upload.

Error codes are the fixed contract — no stack trace and no file path is ever returned:

`INPUT_NOT_FASTQ`, `INPUT_EMPTY`, `INPUT_CORRUPT`, `INPUT_TOO_FEW_READS`, `UPLOAD_TOO_LARGE`,
`ASSEMBLY_TIMEOUT`, `TOOL_FAILED`, `MODEL_MISSING`, `DB_UNAVAILABLE`, `NOT_FOUND`.

Every code has a plain-language message (`pipeline/errors.py`). A failure stops the run, marks that
step failed, leaves later steps `pending`, and sets the isolate to `failed`.

---

## 4. MongoDB

**This is a user-space MongoDB in WSL, not a Windows service** (deviation from the build plan, which
assumed the Windows MSI):

```bash
~/mongodb/bin/mongod --dbpath ~/mongodb_data --bind_ip 127.0.0.1 --port 27017 \
                     --logpath ~/mongodb_logs/mongod.log --fork
```

`scripts/serve.sh` starts it automatically when port 27017 is closed. Version 8.0.4, SSPL v1.

Database `amr_firewall` (tests use `amr_firewall_test`), collections:
`isolates`, `runs`, `assemblies`, `gene_calls`, `predictions`, `benchmarks`, `models`,
`training_isolates`. Indexes are created by `db/indexes.py` at startup and by the seed script.

If MongoDB is down, `/api/health` reports `{"mongo": "down"}` and uploads return `DB_UNAVAILABLE` —
the server itself keeps running.

---

## 5. Tests and the end-to-end check (measured)

```
$ python -m pytest -q
......................................                                   [100%]
38 passed in 32.72s
```

Coverage: Mongo round-trip, input validation, **every error code** (including `ASSEMBLY_TIMEOUT`,
`TOOL_FAILED`, `MODEL_MISSING`, `DB_UNAVAILABLE`, `UPLOAD_TOO_LARGE`), the API response shapes, the
runner failure path, and upload cleanup.

```
$ python scripts/e2e_check.py --url http://127.0.0.1:8000
PASS  1 server starts and answers /api/health
PASS  2 health shape + mongo ok + species list
PASS  2b GET / serves the page
PASS  3 upload synthetic reads -> clean failure (not 'done')  -- TOOL_FAILED: An analysis tool
      stopped with an error, so the run was halted. ... | qc=done assemble=failed amr_calls=pending
      deeparg_calls=pending features=pending predict=pending
FAIL  4 fallback-example results (6 drugs)  -- HTTP 404 NOT_FOUND
PASS  5 non-FASTQ upload -> INPUT_NOT_FASTQ
PASS  6 /api/benchmarks shape -- 0 rows
1 STEP(S) FAILED
```

Step 4 is the one honest failure: `fallback-example` does not exist yet because it must come from a
**real** completed run (see §7). Nothing is seeded to make this step pass.

`scripts/e2e_check.py` uses `testFiles/*R1*`/`*R2*` if that folder exists; it does not exist here, so
step 3 falls back to 12,000 synthetic read pairs. Random reads cannot be assembled, so the expected
outcome there is a clean terminal state — which is what happens (`TOOL_FAILED` after a successful
`qc`), not a hang or a crash.

### Browser checks (measured, against `http://localhost:8000/`)

| Check | Result |
|---|---|
| Species dropdown | filled from `/api/health`: *Escherichia coli*, *Klebsiella pneumoniae* |
| Empty benchmark state | renders “No benchmarks yet. Run the benchmark script to fill this table.” |
| **Use the example** with no example seeded | failure panel `THE ANALYSIS COULD NOT FINISH`, message “No isolate with that id was found.”, code `NOT_FOUND` |
| **Try again** | clears the panel, both file inputs and the species select (verified: 0 files, select back to “Select species”) |
| Valid-looking FASTQ pair | POST → `202`, then `GET /api/isolates/<id>` polled **every 2 s**; all six steps render in order (`qc → assemble → amr_calls → deeparg_calls → features → predict`), ending in a clean `failed` state with `INPUT_TOO_FEW_READS` and its plain message; polling stops at the terminal state |
| 375 px width | no page-level horizontal scroll (`scrollWidth` 359 ≤ 375); the form stacks cleanly and the benchmark table scrolls sideways inside its own box |
| Offline | 13 resources loaded, **all same-origin** (`http://localhost:8000`), zero external requests or tags |

---

## 6. Model and benchmarks

One XGBoost binary classifier per drug (`max_depth=4`, 200 trees, defaults otherwise, no
hyper-parameter search). Features per isolate (`pipeline/features.py`): gene presence from the union
of AMRFinderPlus + DeepARG (`g:`), gene family (`f:`), AMRFinderPlus drug-class/subclass counts
(`c:`, `s:`), best DeepARG probability per class (`d:`), species one-hot (`sp:`).

Rules enforced in code: training uses **only** rows with `split == "train"` in `data/splits.csv`
(`model/train.py`), benchmarking uses **only** `split == "test"` (`model/benchmark.py`). The
AMRFinderPlus baseline in the benchmark is subclass-aware (`pipeline/rules.py`) so it is not an
artificially weak comparison.

### Benchmark results (measured)

Produced from the `benchmarks` collection: 28 public BV-BRC isolates, 20 used for training and the
8 held-out ones for this table, laboratory phenotypes only. `n_test` is the number of test isolates
that had a laboratory call for that drug.

| Drug | Model | Sensitivity | Specificity | n_test |
|---|---|---|---|---|
| ciprofloxacin | amrfinder | 1.000 | 1.000 | 7 |
| ciprofloxacin | custom | 1.000 | 1.000 | 7 |
| ceftazidime | amrfinder | 1.000 | 1.000 | 8 |
| ceftazidime | custom | 1.000 | 1.000 | 8 |
| meropenem | amrfinder | 1.000 | 0.800 | 6 |
| meropenem | custom | **0.000** | 1.000 | 6 |
| gentamicin | amrfinder | 1.000 | 1.000 | 7 |
| gentamicin | custom | 1.000 | 1.000 | 7 |
| tetracycline | amrfinder | 1.000 | **0.333** | 7 |
| tetracycline | custom | 1.000 | 1.000 | 7 |
| trimethoprim-sulfamethoxazole | amrfinder | 1.000 | 1.000 | 5 |
| trimethoprim-sulfamethoxazole | custom | 1.000 | 1.000 | 5 |

**How to read this honestly.** The custom model does **not** beat the AMRFinderPlus rule overall:
the two are identical on four drugs, the model is clearly better on tetracycline specificity
(1.000 vs 0.333), and it is clearly worse on meropenem (sensitivity 0.000 vs 1.000) — meropenem
has only **3** resistant training isolates, so that model is not usable and the rule should be
preferred there. Every test set is small (5–8 isolates), so the intervals are wide and these
numbers are indicative, not validated. The only defensible claim today is the table above.

To reproduce them:

```bash
python scripts/fetch_bvbrc_data.py                 # public BV-BRC phenotypes + genome FASTA
python -m model.extract_features                   # AMRFinderPlus + DeepARG on every listed genome
python -m model.train                              # train on split == "train" only
python -m model.benchmark                          # score on split == "test" only
python scripts/seed_mongo.py
```

Both data sources are public and declared in `DECLARATIONS.md`. Only BV-BRC records with
`evidence == "Laboratory Method"` are used; the computational predictions in the same table are
discarded. The split is fixed by seed 42 in `fetch_bvbrc_data.py`, and the accession column is kept
as text (`562.7707` would otherwise be read as a float and could not be used as a path or key).

---

## 7. The example shown on the page

`fallback-example` is a precomputed copy of an earlier **pipeline run**, not a mock. The path is:

```bash
bash scripts/build_example.sh        # extract features -> train -> benchmark -> seed -> demo run -> export -> seed
```

or by hand:

```bash
# 1. upload a read pair through the page (or POST /api/isolates) and let it finish
# 2. export that finished run:
python scripts/fallback_run.py <isolate_id>        # -> data/fallback_example.json
python scripts/seed_mongo.py                        # -> isolate "fallback-example"
```

`scripts/fallback_run.py` refuses to export anything that is not a completed run with six `done`
steps and predictions, and the seeded copy is tagged with the warning that it is precomputed. Until
a run is exported, the page's example flow shows `NOT_FOUND`.

### The demo read pair is simulated, and this is deliberate

A real archive run could not be used here. Downloading a public FASTQ from ENA works for a couple of
megabytes and is then throttled to a few kilobytes per second on this network: a measured full-file
transfer of `DRR002018` (~121 MB per mate) settled at **~2.6–21 KB/s**, which is hours for one
isolate, and `ftp.bvbrc.org` was unreachable. `scripts/fetch_reads_subset.py` is therefore kept for
the case where the network behaves (it takes the first N read pairs of an ENA run with HTTP byte
ranges), and `scripts/build_example.sh` falls back to **ART-simulated reads from a real BV-BRC
genome** so that the demo path is exercised end to end.

What that means for honesty: the example's reads are not a clinical sample — they are simulated
from a real genome whose laboratory phenotypes are in `data/splits.csv`. The gene calls, the features
and the model prediction itself are all real; only the sequencing reads are synthetic, and the
example carries a warning saying so. Nothing about the **benchmark** depends on this: benchmark rows
come from real BV-BRC genomes and their laboratory phenotypes on the held-out split.

---

## 8. Known limits (all true today)

- **Models and benchmarks exist now** (six XGBoost artefacts in `model/artifacts/`, 12 rows in
  `benchmarks`); see §6 for the measured numbers and their caveats. What is still missing is
  **`fallback-example`**, which needs a completed demo run — until `scripts/build_example.sh`
  finishes it (or it is seeded from any finished run), `/api/isolates/fallback-example/results`
  answers `NOT_FOUND`.
- **Read length outside 50–300 bases is reported as `INPUT_CORRUPT`** with a specific message about
  the observed lengths; the agreed code list has no dedicated code for it.
- **AMRFinderPlus TSV parsing** (`pipeline/amr_calls.py`) follows the AMRFinderPlus 4.x column names
  and accepts the older spellings as well (`Element symbol`/`Gene symbol`, …). Verified on real
  assemblies: over the first 14 BV-BRC genomes in `data/splits.csv` it produced 124 gene calls with
  correct fields (e.g. `gyrA_S83L`, type `AMR`, class `QUINOLONE`, coverage 99.66, identity 98.86) and
  the DeepARG mapping parser produced 975 calls (gene, class, probability). Both are the same code
  paths the live pipeline uses.
- **The test suite is timing-sensitive under heavy CPU load.** With four `tblastn` processes and
  DeepARG saturating the machine, one `pytest -q` run reported 8 failures in the timeout-based
  contract tests; an idle-machine run is 38 passed. Run `pytest -q` when the pipeline is not busy.
- **RAM is the binding constraint.** The host has 7.8 GB and WSL sees 3.8 GB, so SPAdes `-m 6` is a
  cap, not a guarantee. `.wslconfig` was not changed.
- **`testFiles/` does not exist**, so the e2e script uses synthetic reads and cannot exercise
  `assemble` onward without a real pair.
- The frontend is served unmodified from the team's own repository; the API was written to the fixed
  contract rather than reshaping the page to the API.

### Deviations from the brief

| Item | Brief | What was done |
|---|---|---|
| MongoDB | Windows service via MSI | user-space `~/mongodb` inside WSL, started by `scripts/serve.sh` |
| RAM | 16 GB, WSL capped at 10 GB | host has 7.8 GB, WSL sees 3.8 GB; no `.wslconfig` change |
| Dataset | training data assumed present | nothing was on disk; `scripts/fetch_bvbrc_data.py` was written to build it from public BV-BRC data |
| ftp.bvbrc.org | genomes from the FTP site | unreachable from this network; genomes are fetched through the BV-BRC `genome_sequence` API instead |

---

## 9. Repository layout

```
config/settings.py        all paths, limits, timeouts, drugs, species, tool locations
db/mongo.py               every pymongo call; connectivity errors -> DB_UNAVAILABLE
db/indexes.py             idempotent index creation
pipeline/errors.py        the error-code contract and its plain-language messages
pipeline/proc.py          run_tool(): subprocess, timeout kills the process group, TOOL_FAILED
pipeline/qc.py            FASTQ/gzip validation + fastp
pipeline/assemble.py      SPAdes + contig stats
pipeline/amr_calls.py     AMRFinderPlus (+ TSV parser)
pipeline/deeparg_calls.py prodigal ORFs + DeepARG (+ mapping parser)
pipeline/features.py      sparse feature builder
pipeline/rules.py         AMRFinderPlus baseline rule for the benchmark
pipeline/predict.py       model loading, probability -> R/S/uncertain, evidence genes
pipeline/runner.py        the six steps, status bookkeeping, cleanup
api/main.py               FastAPI app, endpoints, static site, error handlers
api/jobs.py               single worker queue + startup recovery
model/extract_features.py AMRFinderPlus + DeepARG over data/splits.csv genomes
model/train.py            per-drug XGBoost, train split only
model/benchmark.py        test split only, custom vs AMRFinderPlus
scripts/serve.sh          one command: mongod + seed + uvicorn
scripts/seed_mongo.py     idempotent seeding
scripts/fallback_run.py   export a real finished run as the page example
scripts/fetch_bvbrc_data.py  public BV-BRC data -> data/splits.csv + training cache
scripts/fetch_reads_subset.py public ENA run -> first N read pairs (demo reads; see §7)
scripts/build_example.sh  unattended: features -> train -> benchmark -> seed -> demo run -> example
scripts/e2e_check.py      end-to-end PASS/FAIL check
web/index.html            frontend (copied, unmodified)
DECLARATIONS.md           every external component: version, licence, exact command
```
