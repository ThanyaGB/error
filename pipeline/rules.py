"""AMRFinderPlus baseline rule used by the benchmark.

A drug is called resistant if AMRFinderPlus reports any AMR element whose class
(and, where the drug is specific, subclass) matches. The rule is subclass-aware
where AMRFinderPlus provides the information, so the baseline is not artificially weak.
"""

# drug -> list of (class substrings, subclass substrings or None)
RULES = {
    "ciprofloxacin": [({"QUINOLONE"}, None)],
    "ceftazidime": [({"BETA-LACTAM"}, {"CEPHALOSPORIN", "CARBAPENEM"})],
    "meropenem": [({"BETA-LACTAM"}, {"CARBAPENEM"})],
    "gentamicin": [({"AMINOGLYCOSIDE"}, {"GENTAMICIN"})],
    "tetracycline": [({"TETRACYCLINE"}, None)],
    "trimethoprim-sulfamethoxazole": [({"TRIMETHOPRIM", "SULFONAMIDE"}, None)],
}


def rule_genes(drug: str, amr_calls: list[dict]) -> list[str]:
    hits = []
    for c in amr_calls:
        if (c.get("element_type") or "AMR").upper() != "AMR":
            continue
        cls, sub = (c.get("drug_class") or ""), (c.get("subclass") or "")
        for classes, subclasses in RULES.get(drug, []):
            if any(k in cls for k in classes) and (subclasses is None or any(k in sub for k in subclasses)):
                if c["gene"] not in hits:
                    hits.append(c["gene"])
    return hits


def amrfinder_rule(drug: str, amr_calls: list[dict]) -> bool:
    return bool(rule_genes(drug, amr_calls))
