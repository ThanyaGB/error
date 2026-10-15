#!/usr/bin/env bash
# One command inside WSL:  bash scripts/serve.sh
# Starts the (user-space) MongoDB if it is not running, then the API on 127.0.0.1:8000.
set -e
cd "$(dirname "$0")/.."
export PATH="$HOME/micromamba/envs/amr/bin:$PATH"
if ! (echo > /dev/tcp/127.0.0.1/27017) 2>/dev/null; then
  if [ -x "$HOME/mongodb/bin/mongod" ]; then
    mkdir -p "$HOME/mongodb_data" "$HOME/mongodb_logs"
    "$HOME/mongodb/bin/mongod" --dbpath "$HOME/mongodb_data" --bind_ip 127.0.0.1 --port 27017 \
        --logpath "$HOME/mongodb_logs/mongod.log" --fork >/dev/null
    sleep 2
  else
    echo "MongoDB is not running on 127.0.0.1:27017 (the page will show a database banner)."
  fi
fi
python scripts/seed_mongo.py || true
exec python -m uvicorn api.main:app --host 127.0.0.1 --port 8000
