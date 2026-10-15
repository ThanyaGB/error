"""Step 1: input validation (server-side, at upload) and read QC with fastp."""
import json
import re
import statistics
import zlib
from pathlib import Path

from config import settings
from pipeline.errors import (INPUT_CORRUPT, INPUT_EMPTY, INPUT_NOT_FASTQ, INPUT_TOO_FEW_READS,
                             PipelineError)
from pipeline.proc import run_tool, tool_version

ALLOWED_SUFFIXES = (".fastq", ".fq", ".fastq.gz", ".fq.gz")
HEAD_BYTES = 8 * 1024 * 1024        # structure + read-length sample
SAMPLE_RECORDS = 20_000
_SEQ_OK = re.compile(rb"^[ACGTNacgtnRYKMSWBDHVrykmswbdhv.\-]+$")


def has_fastq_extension(filename: str) -> bool:
    return (filename or "").lower().endswith(ALLOWED_SUFFIXES)


def _chunks(path: Path):
    """Yield decompressed bytes. Detects gzip by magic number, not by extension.

    A truncated or damaged gzip stream raises INPUT_CORRUPT.
    """
    with open(path, "rb") as fh:
        gz = fh.read(2) == b"\x1f\x8b"
        fh.seek(0)
        if not gz:
            while chunk := fh.read(1 << 20):
                yield chunk
            return
        dec = zlib.decompressobj(31)
        try:
            while raw := fh.read(1 << 20):
                buf = raw
                while buf:
                    if dec.eof:                      # next member of a multi-member gzip
                        dec = zlib.decompressobj(31)
                    out = dec.decompress(buf)
                    if out:
                        yield out
                    buf = dec.unused_data if dec.eof else b""
        except zlib.error:
            raise PipelineError(INPUT_CORRUPT, "gzip stream error") from None
        if not dec.eof:
            raise PipelineError(INPUT_CORRUPT, "gzip stream truncated")


def validate_fastq(path) -> dict:
    """Validate one FASTQ(.gz). Returns {'n_reads', 'median_len', 'max_len'} or raises PipelineError.

    * whole file is read once at C speed (gzip integrity + line count = 4 * reads)
    * the first 8 MB are parsed record by record (header, '+', seq/qual length, bases)
    * read length is judged on a sample of the first 20,000 reads
    """
    path = Path(path)
    if path.stat().st_size == 0:
        raise PipelineError(INPUT_EMPTY, "zero-byte file")

    head = bytearray()
    newlines = 0
    total = 0
    last = b""
    for chunk in _chunks(path):
        if len(head) < HEAD_BYTES:
            head += chunk[: HEAD_BYTES - len(head)]
        newlines += chunk.count(b"\n")
        total += len(chunk)
        last = chunk[-1:]
    if total == 0:
        raise PipelineError(INPUT_EMPTY, "no content after decompression")

    if head[:1] != b"@":
        raise PipelineError(INPUT_NOT_FASTQ, "first byte is not '@'")

    lines = bytes(head).split(b"\n")
    complete_file = total <= HEAD_BYTES
    if not complete_file:
        lines = lines[:-1]                            # last line may be cut by the sample boundary
    lengths = []
    for i in range(0, len(lines) - len(lines) % 4, 4):
        l1, l2, l3, l4 = (x.rstrip(b"\r") for x in lines[i:i + 4])
        if not l1.startswith(b"@") or not l3.startswith(b"+"):
            raise PipelineError(INPUT_NOT_FASTQ if i == 0 else INPUT_CORRUPT, "bad record markers")
        if len(l2) != len(l4) or not l2 or not _SEQ_OK.match(l2):
            raise PipelineError(INPUT_NOT_FASTQ if i == 0 else INPUT_CORRUPT, "bad sequence/quality")
        lengths.append(len(l2))
        if len(lengths) >= SAMPLE_RECORDS:
            break

    n_lines = newlines + (0 if last == b"\n" else 1)
    if n_lines % 4 != 0:
        raise PipelineError(INPUT_CORRUPT, "line count is not a multiple of 4")
    n_reads = n_lines // 4
    if n_reads < settings.MIN_READS:
        raise PipelineError(INPUT_TOO_FEW_READS, f"{n_reads} reads")

    med, mx = statistics.median(lengths), max(lengths)
    if mx > settings.MAX_READ_LEN or med < settings.MIN_READ_LEN:
        raise PipelineError(
            INPUT_CORRUPT,
            f"read length median={med} max={mx}",
            message=f"Read lengths must be between {settings.MIN_READ_LEN} and {settings.MAX_READ_LEN} bases. "
                    f"These reads look like {int(med)}-{mx} bases, which is not supported.",
        )
    return {"n_reads": n_reads, "median_len": int(med), "max_len": mx}


def validate_pair(r1, r2) -> dict:
    a, b = validate_fastq(r1), validate_fastq(r2)
    if a["n_reads"] != b["n_reads"]:
        raise PipelineError(
            INPUT_CORRUPT, "R1/R2 read counts differ",
            message="R1 and R2 do not contain the same number of reads, so they are not a matching pair.",
        )
    return {"r1": a, "r2": b}


# ------------------------------------------------------------------ fastp step
def run_qc(r1: Path, r2: Path, outdir: Path) -> dict:
    """Trim/filter with fastp. Returns the trimmed file paths and read counts."""
    outdir.mkdir(parents=True, exist_ok=True)
    o1, o2 = outdir / "clean_R1.fastq", outdir / "clean_R2.fastq"
    report = outdir / "fastp.json"
    run_tool(
        [settings.FASTP, "-i", r1, "-I", r2, "-o", o1, "-O", o2, "-j", report,
         "-h", outdir / "fastp.html", "--thread", str(settings.THREADS), "--detect_adapter_for_pe"],
        timeout=settings.TOOL_TIMEOUT_S, log_path=outdir / "fastp.log", name="fastp",
    )
    try:
        summary = json.loads(report.read_text())["summary"]
        before = summary["before_filtering"]["total_reads"]
        after = summary["after_filtering"]["total_reads"]
    except (OSError, KeyError, ValueError):
        from pipeline.errors import TOOL_FAILED
        raise PipelineError(TOOL_FAILED, "fastp report unreadable") from None
    if after < settings.MIN_READS:
        raise PipelineError(INPUT_TOO_FEW_READS, f"{after} reads left after QC")
    return {"r1": o1, "r2": o2, "reads_before": before, "reads_after": after,
            "tool_version": tool_version([settings.FASTP, "--version"])}
