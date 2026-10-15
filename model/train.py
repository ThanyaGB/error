"""Train one XGBoost classifier per drug using ONLY rows with split == "train" in data/splits.csv.

usage: python -m model.train
Needs model.extract_features to have produced <cache>/<accession>/calls.json for the training isolates.
"""
import json
import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd
import xgboost as xgb

from config import settings
from db import mongo
from pipeline.features import build_features, vectorize

PARAMS = {"max_depth": 4, "n_estimators": 200}      # as per the build plan; no hyper-parameter search
# The contract asks for one model per drug, and the public laboratory-phenotype data available for
# some drugs is thin (a few resistant isolates). Train whenever there are at least this many isolates
# in each class, and record how many were used (n_resistant / n_susceptible) so the weak ones are
# visible in the benchmark rather than silently missing.
MIN_PER_CLASS = 3


def load_features(accession: str):
    f = settings.TRAIN_CACHE_DIR / accession / "calls.json"
    if not f.exists():
        return None
    c = json.loads(f.read_text())
    sp = c["species"] if c["species"] in settings.SUPPORTED_SPECIES else None
    return build_features(c["amr"], c["deeparg"], sp)


def read_splits() -> pd.DataFrame:
    """Load data/splits.csv with the identifiers kept as text.

    Accessions are text by contract, but a value like "562.7707" would otherwise be inferred as a
    float and could no longer be used as a path or a database key.
    """
    return pd.read_csv(settings.SPLITS_CSV, dtype={"accession": str, "species": str, "split": str})


def labelled_rows(split: str) -> pd.DataFrame:
    df = read_splits()
    return df[df["split"] == split]


def main() -> int:
    train = labelled_rows("train")                 # the ONLY rows this script may use
    feats = {r.accession: load_features(r.accession) for r in train.itertuples()}
    feats = {k: v for k, v in feats.items() if v is not None}
    if not feats:
        print("no extracted training features found; run model.extract_features first", file=sys.stderr)
        return 1
    version = "xgb-" + datetime.now(timezone.utc).strftime("%Y%m%d-%H%M")
    settings.ARTIFACT_DIR.mkdir(parents=True, exist_ok=True)
    trained = 0
    for drug in settings.DRUGS:
        if drug not in train.columns:
            print(f"{drug}: no phenotype column, skipped")
            continue
        rows = train[train[drug].isin(["R", "S"]) & train["accession"].isin(feats)]
        n_r, n_s = int((rows[drug] == "R").sum()), int((rows[drug] == "S").sum())
        if min(n_r, n_s) < MIN_PER_CLASS:
            print(f"{drug}: only R={n_r} S={n_s} training labels, skipped")
            continue
        names = sorted({k for acc in rows["accession"] for k in feats[acc]["values"]})
        X = np.array([vectorize(feats[a]["values"], names) for a in rows["accession"]])
        y = (rows[drug] == "R").astype(int).to_numpy()
        clf = xgb.XGBClassifier(**PARAMS, eval_metric="logloss", n_jobs=settings.THREADS, random_state=42)
        clf.fit(X, y)
        path = settings.ARTIFACT_DIR / f"{drug}.json"
        clf.get_booster().save_model(str(path))
        (settings.ARTIFACT_DIR / f"{drug}.meta.json").write_text(json.dumps(
            {"drug": drug, "version": version, "features": names, "trained_on_n": len(rows),
             "n_resistant": n_r, "n_susceptible": n_s, "params": PARAMS}))
        if mongo.is_up():
            mongo.upsert_model(drug, version, names, len(rows), str(path), params=PARAMS,
                               n_resistant=n_r, n_susceptible=n_s)
        print(f"{drug}: trained on {len(rows)} isolates (R={n_r}, S={n_s}), {len(names)} features")
        trained += 1
    return 0 if trained else 1


if __name__ == "__main__":
    raise SystemExit(main())
