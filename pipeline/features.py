"""Step 5: gene presence/absence feature table from AMRFinderPlus + DeepARG calls.

Feature names (all numeric, sparse -- absent means 0):
  g:<key>   1 if a gene is present (union of AMRFinderPlus and DeepARG, names canonicalised)
  f:<fam>   1 if any gene of that family is present (CTX-M-14 and CTX-M-15 both -> f:ctxm)
  c:<class> number of AMRFinderPlus AMR calls in a drug class
  s:<sub>   number of AMRFinderPlus AMR calls in a drug subclass
  d:<class> highest DeepARG probability for a drug class
  sp:<species> one-hot species
"""
import re


def _key(gene: str) -> str:
    k = re.sub(r"[^a-z0-9]", "", gene.lower())
    return re.sub(r"^bla(?=[a-z])", "", k)


def _family(key: str) -> str:
    return re.sub(r"\d+$", "", key)


def _norm_class(c: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", (c or "").lower()).strip("_")


def build_features(amr_calls: list[dict], deeparg_calls: list[dict], species: str | None = None) -> dict:
    """Return {'values': {feature: float}, 'sources': {feature: [display gene, ...]}}."""
    values: dict[str, float] = {}
    sources: dict[str, list[str]] = {}

    def add_source(feat, gene):
        lst = sources.setdefault(feat, [])
        if gene not in lst:
            lst.append(gene)

    for c in amr_calls:
        if (c.get("element_type") or "AMR").upper() != "AMR":
            continue
        gene = c["gene"]
        key = _key(gene)
        values[f"g:{key}"] = 1.0
        add_source(f"g:{key}", gene)
        fam = _family(key)
        if fam and fam != key:
            values[f"f:{fam}"] = 1.0
            add_source(f"f:{fam}", gene)
        for prefix, field in (("c", "drug_class"), ("s", "subclass")):
            for part in re.split(r"[/;,]", c.get(field) or ""):
                nm = _norm_class(part)
                if nm:
                    values[f"{prefix}:{nm}"] = values.get(f"{prefix}:{nm}", 0.0) + 1.0
                    add_source(f"{prefix}:{nm}", gene)

    for c in deeparg_calls:
        gene = c["gene"]
        key = _key(gene)
        values[f"g:{key}"] = 1.0
        add_source(f"g:{key}", gene)
        fam = _family(key)
        if fam and fam != key:
            values[f"f:{fam}"] = 1.0
            add_source(f"f:{fam}", gene)
        cls = _norm_class(c.get("drug_class"))
        if cls:
            p = c.get("probability") or 0.0
            values[f"d:{cls}"] = max(values.get(f"d:{cls}", 0.0), float(p))
            add_source(f"d:{cls}", gene)

    if species:
        values[f"sp:{species}"] = 1.0
    return {"values": values, "sources": sources}


def vectorize(values: dict[str, float], feature_list: list[str]) -> list[float]:
    """Align an isolate's sparse features to a model's fixed feature order (unseen features dropped)."""
    return [float(values.get(f, 0.0)) for f in feature_list]
