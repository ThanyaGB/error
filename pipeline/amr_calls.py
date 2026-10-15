"""Step 3: resistance gene calls with AMRFinderPlus (command line, TSV output)."""
import csv
import re
from pathlib import Path

from config import settings
from pipeline.errors import TOOL_FAILED, PipelineError
from pipeline.proc import run_tool, tool_version

# AMRFinderPlus renamed some columns between releases; accept either spelling.
_COLS = {
    "gene": ("Element symbol", "Gene symbol"),
    "name": ("Element name", "Sequence name"),
    "element_type": ("Type", "Element type"),
    "subtype": ("Subtype", "Element subtype"),
    "drug_class": ("Class",),
    "subclass": ("Subclass",),
    "method": ("Method",),
    "contig": ("Contig id",),
    "start": ("Start",),
    "stop": ("Stop",),
    "coverage": ("% Coverage of reference", "% Coverage of reference sequence"),
    "identity": ("% Identity to reference", "% Identity to reference sequence"),
}


def _pick(row: dict, keys: tuple):
    for k in keys:
        if k in row and row[k] not in (None, "", "NA"):
            return row[k]
    return None


def _num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def parse_amrfinder_tsv(path: Path) -> list[dict]:
    """Parse an AMRFinderPlus report into gene_calls rows (no isolate_id/source yet)."""
    calls = []
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            gene = _pick(row, _COLS["gene"])
            if not gene:
                continue
            calls.append({
                "gene": gene,
                "name": _pick(row, _COLS["name"]),
                "element_type": _pick(row, _COLS["element_type"]),
                "subtype": _pick(row, _COLS["subtype"]),
                "drug_class": (_pick(row, _COLS["drug_class"]) or "").upper() or None,
                "subclass": (_pick(row, _COLS["subclass"]) or "").upper() or None,
                "method": _pick(row, _COLS["method"]),
                "contig": _pick(row, _COLS["contig"]),
                "start": _num(_pick(row, _COLS["start"])),
                "stop": _num(_pick(row, _COLS["stop"])),
                "identity": _num(_pick(row, _COLS["identity"])),
                "coverage": _num(_pick(row, _COLS["coverage"])),
            })
    return calls


def amrfinder_versions() -> tuple[str, str]:
    sw = tool_version([settings.AMRFINDER, "--version"])
    db = tool_version([settings.AMRFINDER, "--database_version"], r"Database version:\s*(\S+)")
    return sw, db


def run_amrfinder(contigs: Path, species: str, outdir: Path) -> dict:
    """Run AMRFinderPlus on contigs. `species` is the lower_snake key; unsupported -> no --organism."""
    outdir.mkdir(parents=True, exist_ok=True)
    tsv = outdir / "amrfinder.tsv"
    cmd = [settings.AMRFINDER, "-n", contigs, "-o", tsv, "--threads", str(settings.THREADS)]
    organism = settings.SUPPORTED_SPECIES.get(species)
    if organism:
        cmd += ["--organism", organism]
    if settings.AMRFINDER_DB:
        cmd += ["--database", settings.AMRFINDER_DB]
    run_tool(cmd, timeout=settings.TOOL_TIMEOUT_S, log_path=outdir / "amrfinder.log", name="AMRFinderPlus")
    if not tsv.exists():
        raise PipelineError(TOOL_FAILED, "AMRFinderPlus wrote no report")
    sw, db = amrfinder_versions()
    return {"calls": parse_amrfinder_tsv(tsv), "tool_version": sw, "db_version": db,
            "organism": organism}
