# DECLARATIONS

Every component below is **external work used as-is**. None of it is our code. Each tool is called
as a command-line program (`pipeline/proc.py: run_tool`); no tool is linked or imported as a library.

Verified on this machine on 2026-10-15 (WSL2 Ubuntu, micromamba env `amr`, Python 3.11.17).

---

## 1. External programs and data

### fastp — read quality control
| | |
|---|---|
| Version | 1.4.0 |
| URL | https://github.com/OpenGene/fastp |
| Licence | MIT |
| Purpose | Adapter detection, quality/length trimming of the read pair (step `qc`) |
| Called from | `pipeline/qc.py` → `run_qc()` |

```bash
fastp -i <R1> -I <R2> -o <out>/clean_R1.fastq -O <out>/clean_R2.fastq \
      -j <out>/fastp.json -h <out>/fastp.html --thread 4 --detect_adapter_for_pe
```

### SPAdes — genome assembly
| | |
|---|---|
| Version | 4.3.0 (`spades.py --version`) |
| URL | https://github.com/ablab/spades |
| Licence | GPL-2.0 |
| Purpose | De novo assembly of the trimmed pair, isolate mode (step `assemble`) |
| Called from | `pipeline/assemble.py` → `run_assembly()` |

```bash
spades.py --isolate -1 <clean_R1> -2 <clean_R2> -t 4 -m 6 -o <out>/spades
```
`-t 4` and `-m 6` are the build-plan caps (4 threads, 6 GB). The 40-minute timeout is enforced by
our own process handler, not by SPAdes.

### AMRFinderPlus (+ its database) — resistance gene detection
| | |
|---|---|
| Version | software 4.2.7 · database `2026-08-07.1` |
| URL | https://github.com/ncbi/amr |
| Licence | Public domain (United States Government work, NCBI) |
| Database licence | AMRFinderPlus database is distributed by NCBI under the same public-domain terms; the `--database` directory shipped by the bioconda package is used unmodified |
| Purpose | Curated gene/point-mutation calls on the contigs (step `amr_calls`) |
| Called from | `pipeline/amr_calls.py` → `run_amrfinder()` |

```bash
amrfinder -n <contigs.fasta> -o <out>/amrfinder.tsv --threads 4 \
          [--organism Escherichia | --organism Klebsiella_pneumoniae] \
          [--database <AMRFINDER_DB>]
```
`--organism` is only added when the uploaded species is in our supported list. Version/DB strings are
read with `amrfinder --version` and `amrfinder --database_version`.

### prodigal — gene prediction (feeds DeepARG)
| | |
|---|---|
| Version | 2.6.3 |
| URL | https://github.com/hyattpd/Prodigal |
| Licence | GPL-3.0 |
| Purpose | Turn assembled contigs into gene-length nucleotide ORFs, which is the input DeepARG actually accepts (step `deeparg_calls`) |
| Called from | `pipeline/deeparg_calls.py` → `predict_orfs()` |

```bash
prodigal -i <contigs.fasta> -d <out>/orfs.fna -a <out>/orfs.faa \
         -o <out>/prodigal.gff -f gff -q
```

### DeepARG — deep-learning resistance gene calls
| | |
|---|---|
| Version | `deeparg-modern` 2.0.0 (command `deeparg`) |
| URL | https://github.com/gaarangoa/deeparg (modernised fork: `gaarangoa/deeparg-modern`) |
| Licence | MIT |
| Model bundle | Hugging Face `gaarangoa/deeparg`, snapshot `07d2eaa78502f41725c4222de1288b8b2f0d6259`, model `LS` (long sequences) |
| Purpose | Second, independent resistance-gene source, run on the prodigal ORFs (step `deeparg_calls`) |
| Called from | `pipeline/deeparg_calls.py` → `run_deeparg()` |

```bash
deeparg predict --model LS --type nucl -i <out>/orfs.fna -o <out>/deeparg --hf-device cpu
```
Runs on CPU, offline: our child-process environment sets `HF_HUB_OFFLINE=1` and `TRANSFORMERS_OFFLINE=1`,
so the cached bundle is used and no network access is attempted. Only the high-confidence
`*.mapping.ARG` output is parsed.

### DIAMOND — sequence aligner used by DeepARG
| | |
|---|---|
| Version | 2.2.8 (`diamond --version`) |
| URL | https://github.com/bbuchfink/diamond |
| Licence | GPL-3.0 |
| Purpose | Alignment engine invoked **inside** `deeparg predict`; a copy is bundled with the DeepARG install |
| Called from | not called by our code directly — declared because the DeepARG step depends on it |

### MongoDB Community Server — storage
| | |
|---|---|
| Version | 8.0.4 (`mongod --version`) |
| URL | https://www.mongodb.com/try/download/community |
| Licence | Server Side Public License (SSPL) v1 |
| Purpose | Stores isolates, runs, assemblies, gene_calls, predictions, benchmarks, models, training_isolates |
| Called from | `db/mongo.py` (PyMongo driver); started by `scripts/serve.sh` |

```bash
# user-space install in WSL, NOT a Windows service (see README "Deviations")
~/mongodb/bin/mongod --dbpath ~/mongodb_data --bind_ip 127.0.0.1 --port 27017 \
                     --logpath ~/mongodb_logs/mongod.log --fork
```

### Frontend page (copied, not written in this build)
| | |
|---|---|
| Files | `web/index.html`, `web/assets/blood-slide-prep.webp` |
| Source | the team's own repository, https://github.com/ThanyaGB/error |
| Purpose | Static single page served at `/`; talks only to `/api/*` on the same origin |
| Note | Copied in unmodified; the API was built to the fixed contract, the page was not reshaped to the API |

---

## 1b. Public data sources (used as-is, never re-labelled as our own)

### BV-BRC — laboratory resistance phenotypes and genome sequences
| | |
|---|---|
| URL | https://www.bv-brc.org (data APIs `genome_amr`, `genome_sequence`, `genome`) |
| Version | live service, queried 2026-10-15 |
| Terms of use | BV-BRC is an NIAID-funded integrated resource, publicly accessible; per its Database Commons entry it is *free for academic use only, with commercial use requiring a licence* (https://ngdc.cncb.ac.cn/databasecommons/database/id/8261). The AMR phenotypes it republishes come from NCBI BioSample / antibiogram submissions and the genomes from the INSDC archives. |
| Purpose | `data/splits.csv` (laboratory R/S calls for the six drugs, evidence `Laboratory Method` only) and the genome FASTA used to train and test the models |
| Called from | `scripts/fetch_bvbrc_data.py` (`api()`), no CLI — HTTP GET only |

```
https://www.bv-brc.org/api/genome_amr/?eq(taxon_id,562)&eq(evidence,Laboratory Method)
                                  &select(genome_id,antibiotic,resistant_phenotype,evidence)&limit(50000)
https://www.bv-brc.org/api/genome_sequence/?eq(genome_id,<id>)&select(sequence_id,sequence,length)&limit(20000)
```

`ftp.bvbrc.org` was **not** used: it is unreachable from this network, so sequences come from the
`genome_sequence` API instead. Only laboratory measurements are taken; computational predictions in
the same table are discarded.

### ENA (European Nucleotide Archive / INSDC) — the demo read pair
| | |
|---|---|
| URL | https://www.ebi.ac.uk/ena |
| Version | live service, queried 2026-10-15 |
| Terms of use | INSDC "uniform policy of free and unrestricted access to all of the data records"; no use restrictions or licensing requirements are placed on sequence records (https://www.ebi.ac.uk/ena/browser/about/policies) |
| Purpose | the downsampled paired reads used for the one real end-to-end run (first 500,000 read pairs of run `DRR002018`) |
| Called from | `scripts/fetch_reads_subset.py` (`filereport()` then HTTP byte-range GETs) |

The subset is taken as **the first N read pairs** of the run, by reading a prefix of the official
FASTQ.gz files with HTTP range requests; it is not a random sample and the reads are not modified.

---

## 2. Python packages (exact installed versions)

Installed in the micromamba env `amr`; frozen in `requirements.txt`.

| Package | Version | Licence |
|---|---|---|
| fastapi | 0.143.1 | MIT |
| starlette | 1.8.0 | BSD-3-Clause |
| uvicorn | 0.54.0 | BSD-3-Clause |
| pydantic | 2.14.1 | MIT |
| python-multipart | 0.0.32 | Apache-2.0 |
| pymongo | 4.19.0 | Apache-2.0 |
| pandas | 3.0.6 | BSD-3-Clause |
| numpy | 2.4.6 | BSD-3-Clause |
| xgboost | 3.2.0 | Apache-2.0 |
| scikit-learn | 1.9.1 | BSD-3-Clause |
| scipy | 1.17.1 | BSD-3-Clause |
| pytest | 9.1.1 | MIT |
| httpx | 0.28.1 | BSD-3-Clause |
| torch | 2.14.1+cpu | BSD-3-Clause |
| transformers | 5.19.0 | Apache-2.0 |
| biopython | 1.88 | Biopython Licence (BSD-like) |

`python-multipart` is the multipart parser behind `POST /api/isolates`. `torch`/`transformers`/`biopython`
are pulled in by DeepARG.

---

## 3. Cited as background, **not used**

| Project | Why it is not used |
|---|---|
| ARIBA | GPL-3.0 and unmaintained; excluded by the build plan |
| PIMA (its pipeline) | Development stage, no licence; excluded by the build plan |
| wgs2amr | Background reference only |
| BALROG | Background reference only |
| Reads2Resistome | Background reference only |
| ARGNet | Background reference only |

No code from any of the projects above (or from any other project) was copied into this codebase.

---

## 4. What is our own work

Written during this build, with no code copied from elsewhere: `config/`, `db/`, `pipeline/` (process
handling, input validation, QC/assembly/gene-call glue, feature building, rules, prediction, runner),
`api/`, `model/` (feature extraction, training, benchmarking), `scripts/`, `tests/`, and the decision to
feed DeepARG prodigal ORFs rather than raw contigs.
