"""Step 2: genome assembly with SPAdes (--isolate), 40 minute timeout."""
from pathlib import Path

from config import settings
from pipeline.errors import ASSEMBLY_TIMEOUT, TOOL_FAILED, PipelineError
from pipeline.proc import ToolTimeout, run_tool, tool_version


def contig_stats(fasta: Path) -> dict:
    """n_contigs, total_length, N50 from a FASTA file."""
    lengths, cur = [], 0
    with open(fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur:
                    lengths.append(cur)
                cur = 0
            else:
                cur += len(line.strip())
    if cur:
        lengths.append(cur)
    if not lengths:
        return {"n_contigs": 0, "total_length": 0, "n50": 0}
    total = sum(lengths)
    acc = 0
    n50 = 0
    for ln in sorted(lengths, reverse=True):
        acc += ln
        if acc * 2 >= total:
            n50 = ln
            break
    return {"n_contigs": len(lengths), "total_length": total, "n50": n50}


def run_assembly(r1: Path, r2: Path, outdir: Path) -> dict:
    """Assemble paired reads. Returns {'contigs': Path, stats..., 'tool_version'}."""
    spades_dir = outdir / "spades"
    cmd = [settings.SPADES, "--isolate", "-1", r1, "-2", r2,
           "-t", str(settings.THREADS), "-m", str(settings.SPADES_MEM_GB), "-o", spades_dir]
    try:
        run_tool(cmd, timeout=settings.ASSEMBLY_TIMEOUT_S, log_path=outdir / "spades.log", name="SPAdes")
    except ToolTimeout:
        raise PipelineError(ASSEMBLY_TIMEOUT, f"SPAdes exceeded {settings.ASSEMBLY_TIMEOUT_S}s") from None
    contigs = spades_dir / "contigs.fasta"
    if not contigs.exists() or contigs.stat().st_size == 0:
        raise PipelineError(TOOL_FAILED, "SPAdes produced no contigs")
    stats = contig_stats(contigs)
    if stats["n_contigs"] == 0:
        raise PipelineError(TOOL_FAILED, "SPAdes contigs empty")
    return {"contigs": contigs, **stats,
            "tool_version": tool_version([settings.SPADES, "--version"], r"v(\d+(?:\.\d+)+)")}
