"""Step 4: DeepARG on gene-length ORFs.

Input-format decision (from DeepARG's README): `deeparg predict` takes FASTA of
reads (model SS) or of gene-length sequences (model LS) -- it does NOT take
assembled contigs. We therefore:

    contigs.fasta --prodigal--> ORFs (nucleotide FASTA) --deeparg predict --model LS--> <prefix>.mapping.ARG

Only the high-confidence `*.mapping.ARG` file is used (probability >= 0.8).
"""
import csv
import os
from importlib import metadata
from pathlib import Path

from config import settings
from pipeline.errors import MODEL_MISSING, TOOL_FAILED, PipelineError
from pipeline.proc import run_tool


def hf_cache_dir() -> Path:
    hub = os.environ.get("HF_HUB_CACHE") or os.path.join(
        os.environ.get("HF_HOME", str(Path.home() / ".cache" / "huggingface")), "hub")
    return Path(hub) / "models--gaarangoa--deeparg"


def deeparg_model_snapshot() -> Path | None:
    """Newest cached snapshot that contains the configured model folder (LS/SS), else None."""
    snaps = hf_cache_dir() / "snapshots"
    if not snaps.is_dir():
        return None
    for snap in sorted(snaps.iterdir(), key=lambda p: p.stat().st_mtime, reverse=True):
        if (snap / settings.DEEPARG_MODEL).is_dir():
            return snap
    return None


def deeparg_versions() -> tuple[str, str]:
    try:
        sw = metadata.version("deeparg-modern")
    except metadata.PackageNotFoundError:
        sw = "unknown"
    snap = deeparg_model_snapshot()
    return sw, (f"gaarangoa/deeparg@{snap.name[:8]} (v2 {settings.DEEPARG_MODEL})" if snap else "missing")


def predict_orfs(contigs: Path, outdir: Path) -> Path:
    orfs = outdir / "orfs.fna"
    run_tool([settings.PRODIGAL, "-i", contigs, "-d", orfs, "-a", outdir / "orfs.faa",
              "-o", outdir / "prodigal.gff", "-f", "gff", "-q"],
             timeout=settings.TOOL_TIMEOUT_S, log_path=outdir / "prodigal.log", name="prodigal")
    if not orfs.exists() or orfs.stat().st_size == 0:
        raise PipelineError(TOOL_FAILED, "prodigal found no ORFs")
    return orfs


def parse_mapping_arg(path: Path) -> list[dict]:
    """Parse DeepARG `*.mapping.ARG` into gene_calls rows.

    Columns: #ARG, query-start, query-end, read_id, predicted_ARG-class, best-hit,
    probability, identity, alignment-length, ... ; best-hit looks like
    `ACCESSION|FEATURES|gene|class|gene`.
    """
    calls = []
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            best = row.get("best-hit") or ""
            parts = best.split("|")
            gene = parts[2] if len(parts) >= 3 and parts[2] else (row.get("#ARG") or "")
            if not gene:
                continue
            query = row.get("read_id") or ""

            def num(k):
                try:
                    return float(row[k])
                except (KeyError, TypeError, ValueError):
                    return None

            calls.append({
                "gene": gene,
                "drug_class": (row.get("predicted_ARG-class") or "").lower() or None,
                "identity": num("identity"),
                "probability": num("probability"),
                "contig": query.rsplit("_", 1)[0] if "_" in query else query,
                "start": num("query-start"),
                "stop": num("query-end"),
                "best_hit": best,
            })
    return calls


def run_deeparg(contigs: Path, outdir: Path) -> dict:
    outdir.mkdir(parents=True, exist_ok=True)
    if deeparg_model_snapshot() is None:
        raise PipelineError(MODEL_MISSING, "DeepARG model not in local Hugging Face cache")
    orfs = predict_orfs(contigs, outdir)
    prefix = outdir / "deeparg"
    run_tool([settings.DEEPARG, "predict", "--model", settings.DEEPARG_MODEL, "--type", "nucl",
              "-i", orfs, "-o", prefix, "--hf-device", "cpu"],
             timeout=settings.DEEPARG_TIMEOUT_S, log_path=outdir / "deeparg.log", name="DeepARG")
    mapping = Path(str(prefix) + ".mapping.ARG")
    if not mapping.exists():
        raise PipelineError(TOOL_FAILED, "DeepARG wrote no mapping file")
    sw, db = deeparg_versions()
    return {"calls": parse_mapping_arg(mapping), "tool_version": sw, "db_version": db}
