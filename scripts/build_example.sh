#!/usr/bin/env bash
# Unattended build of the data, the models and the demo example. Run inside WSL:
#
#   cd ~/proj && bash scripts/build_example.sh
#
# Steps (each is safe to re-run; extraction skips genomes that already have calls.json):
#   1. model.extract_features   AMRFinderPlus + DeepARG on every genome in data/splits.csv
#   2. model.train              one XGBoost model per drug, split == "train" only
#   3. model.benchmark          custom vs AMRFinderPlus rule, split == "test" only
#   4. scripts/seed_mongo.py    load splits + models into MongoDB
#   5. demo reads               ART-simulated pair from a real BV-BRC genome (see README)
#   6. one real pipeline run    POST the pair to the running API and follow it to the end
#   7. scripts/fallback_run.py  export that finished run
#   8. scripts/seed_mongo.py    seed it as "fallback-example"
#
# Log: ~/amr_work/build_example.log
set -u
cd "$(dirname "$0")/.."
export PATH="$HOME/micromamba/envs/amr/bin:$PATH"
export HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1
WORK="$HOME/amr_work"
READS="$WORK/demo_reads"
API="http://127.0.0.1:8000"
mkdir -p "$WORK" "$READS"

echo "=== 1. extract features $(date -Is) ==="
python -m model.extract_features || echo "extract_features reported problems (see above)"

echo "=== 2. train $(date -Is) ==="
python -m model.train || echo "train reported problems (see above)"

echo "=== 3. benchmark $(date -Is) ==="
python -m model.benchmark || echo "benchmark reported problems (see above)"

echo "=== 4. seed $(date -Is) ==="
python scripts/seed_mongo.py || echo "seed reported problems (see above)"

echo "=== 5. demo reads $(date -Is) ==="
GENOME=$(python - <<'PY'
import csv, sys
from config import settings
rows = list(csv.DictReader(open(settings.SPLITS_CSV)))
drugs = list(settings.DRUGS)
def labelled(r): return sum(1 for d in drugs if r.get(d) in ("R", "S"))
# prefer the held-out split, and the isolate with the most laboratory calls
rows.sort(key=lambda r: (r["split"] != "test", -labelled(r), r["accession"]))
print(rows[0]["accession"])
PY
)
echo "simulating reads from genome $GENOME"
FASTA="$WORK/training/$GENOME/genome.fna"
if [ ! -s "$READS/demo_R1.fastq.gz" ]; then
  art_illumina -ss HS25 -i "$FASTA" -p -l 150 -f 29 -m 400 -s 10 -o "$READS/demo_" -na > "$READS/art.log" 2>&1
  mv -f "$READS/demo_1.fq" "$READS/demo_R1.fastq"
  mv -f "$READS/demo_2.fq" "$READS/demo_R2.fastq"
  gzip -f "$READS/demo_R1.fastq" "$READS/demo_R2.fastq"
  echo "generated:"; ls -la "$READS"/demo_R*.fastq.gz
  echo "read pairs: $(( $(zcat "$READS/demo_R1.fastq.gz" | wc -l) / 4 ))"
fi

echo "=== 6. one pipeline run $(date -Is) ==="
if ! curl -fsS "$API/api/health" > /dev/null; then
  echo "API is not answering on $API; start it with: bash scripts/serve.sh" >&2
  exit 1
fi
RESP=$(curl -sS -X POST "$API/api/isolates" \
  -F "r1=@$READS/demo_R1.fastq.gz;filename=demo_R1.fastq.gz" \
  -F "r2=@$READS/demo_R2.fastq.gz;filename=demo_R2.fastq.gz" \
  -F "species=escherichia_coli")
echo "upload: $RESP"
ID=$(printf '%s' "$RESP" | python -c "import json,sys; print(json.load(sys.stdin).get('isolate_id',''))")
[ -n "$ID" ] || { echo "upload rejected; stopping" >&2; exit 1; }

for _ in $(seq 1 900); do
  STATE=$(curl -sS "$API/api/isolates/$ID")
  STATUS=$(printf '%s' "$STATE" | python -c "import json,sys; print(json.load(sys.stdin)['status'])")
  [ "$STATUS" = "done" ] || [ "$STATUS" = "failed" ] && break
  sleep 4
done
echo "isolate $ID finished as $STATUS"
printf '%s\n' "$STATE" | python -m json.tool

echo "=== 7. export the example $(date -Is) ==="
python scripts/fallback_run.py "$ID" || echo "export refused (the run must be done with 6 done steps)"

echo "=== 8. seed the example $(date -Is) ==="
python scripts/seed_mongo.py

echo "=== BUILD DONE $(date -Is) ==="
