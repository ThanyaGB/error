"""Step 6: per-drug XGBoost prediction with an R / S / uncertain call."""
import json
from pathlib import Path

from config import settings
from pipeline.errors import MODEL_MISSING, PipelineError
from pipeline.features import vectorize

_cache: dict[str, tuple] = {}


def call_from_probability(p: float) -> str:
    if p >= settings.R_THRESHOLD:
        return "R"
    if p <= settings.S_THRESHOLD:
        return "S"
    return "uncertain"


def artifact_paths(drug: str) -> tuple[Path, Path]:
    return settings.ARTIFACT_DIR / f"{drug}.json", settings.ARTIFACT_DIR / f"{drug}.meta.json"


def load_model(drug: str):
    """Return (booster, meta). Missing or unreadable artifact -> MODEL_MISSING."""
    model_path, meta_path = artifact_paths(drug)
    key = str(model_path)
    if not model_path.exists() or not meta_path.exists():
        raise PipelineError(MODEL_MISSING, f"no artifact for {drug}")
    cached = _cache.get(key)
    if cached and cached[2] == model_path.stat().st_mtime:
        return cached[0], cached[1]
    try:
        import xgboost as xgb
        booster = xgb.Booster()
        booster.load_model(str(model_path))
        meta = json.loads(meta_path.read_text())
    except Exception:
        raise PipelineError(MODEL_MISSING, f"artifact for {drug} unreadable") from None
    _cache[key] = (booster, meta, model_path.stat().st_mtime)
    return booster, meta


def predict_drug(drug: str, feats: dict) -> dict:
    import numpy as np
    import xgboost as xgb

    booster, meta = load_model(drug)
    names = meta["features"]
    x = np.array([vectorize(feats["values"], names)], dtype=float)
    dm = xgb.DMatrix(x, feature_names=names)
    prob = float(booster.predict(dm)[0])
    prob = min(1.0, max(0.0, prob))
    call = call_from_probability(prob)

    evidence: list[str] = []
    if call != "S":
        contribs = booster.predict(dm, pred_contribs=True)[0][:-1]
        order = sorted(range(len(names)), key=lambda i: contribs[i], reverse=True)
        for i in order:
            if contribs[i] <= 0 or x[0][i] == 0:
                continue
            for gene in feats["sources"].get(names[i], []):
                if gene not in evidence:
                    evidence.append(gene)
            if len(evidence) >= 6:
                break
    return {"drug": drug, "drug_class": settings.DRUGS.get(drug), "call": call, "probability": prob,
            "model_version": meta["version"], "evidence": evidence[:6]}


def predict_all(feats: dict) -> list[dict]:
    """Predict every configured drug that has a trained model. Raises MODEL_MISSING if none do."""
    out = []
    for drug in settings.DRUGS:
        model_path, _ = artifact_paths(drug)
        if model_path.exists():
            out.append(predict_drug(drug, feats))
    if not out:
        raise PipelineError(MODEL_MISSING, "no trained models in artifact dir")
    return out
