"""Benchmark on rows with split == "test" ONLY: custom XGBoost model vs the AMRFinderPlus rule.

usage: python -m model.benchmark
Sensitivity = TP/(TP+FN) for resistant isolates; specificity = TN/(TN+FP) for susceptible ones.
The custom model is scored at probability >= 0.5 (R) so the comparison is binary; the share of test
isolates that the product would have labelled "uncertain" (0.3 < p < 0.7) is stored as `fraction_uncertain`.
"""
import json
import sys

import numpy as np
import xgboost as xgb

from config import settings
from db import mongo
from model.train import labelled_rows, load_features
from pipeline import predict
from pipeline.errors import PipelineError
from pipeline.features import vectorize
from pipeline.rules import amrfinder_rule


def rates(y_true, y_pred):
    y_true, y_pred = np.array(y_true), np.array(y_pred)
    tp = int(((y_true == 1) & (y_pred == 1)).sum())
    fn = int(((y_true == 1) & (y_pred == 0)).sum())
    tn = int(((y_true == 0) & (y_pred == 0)).sum())
    fp = int(((y_true == 0) & (y_pred == 1)).sum())
    sens = tp / (tp + fn) if tp + fn else None
    spec = tn / (tn + fp) if tn + fp else None
    return sens, spec


def main() -> int:
    test = labelled_rows("test")                  # the ONLY rows this script may use
    calls = {}
    for r in test.itertuples():
        f = settings.TRAIN_CACHE_DIR / r.accession / "calls.json"
        if f.exists():
            calls[r.accession] = json.loads(f.read_text())
    if not calls:
        print("no extracted test features found; run model.extract_features first", file=sys.stderr)
        return 1
    written = 0
    for drug in settings.DRUGS:
        if drug not in test.columns:
            continue
        rows = test[test[drug].isin(["R", "S"]) & test["accession"].isin(calls)]
        if rows.empty:
            continue
        y = [1 if v == "R" else 0 for v in rows[drug]]
        rule_pred = [int(amrfinder_rule(drug, calls[a]["amr"])) for a in rows["accession"]]
        s, p = rates(y, rule_pred)
        mongo.insert_benchmark(drug, "amrfinder", s, p, len(rows), n_resistant=sum(y), n_susceptible=len(y) - sum(y))
        written += 1
        try:
            booster, meta = predict.load_model(drug)
        except PipelineError:
            print(f"{drug}: no trained model; only the AMRFinderPlus rule was benchmarked")
            continue
        feats = [load_features(a) for a in rows["accession"]]
        X = np.array([vectorize(f["values"], meta["features"]) for f in feats])
        prob = booster.predict(xgb.DMatrix(X, feature_names=meta["features"]))
        s2, p2 = rates(y, (prob >= 0.5).astype(int))
        unc = float(((prob > settings.S_THRESHOLD) & (prob < settings.R_THRESHOLD)).mean())
        mongo.insert_benchmark(drug, "custom", s2, p2, len(rows), n_resistant=sum(y),
                               n_susceptible=len(y) - sum(y), threshold=0.5, fraction_uncertain=unc,
                               model_version=meta["version"])
        print(f"{drug}: n_test={len(rows)}  rule sens={s} spec={p} | custom sens={s2} spec={p2}")
    return 0 if written else 1


if __name__ == "__main__":
    raise SystemExit(main())
