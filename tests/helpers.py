import gzip
import random
from pathlib import Path


def write_fastq(path: Path, n_reads: int = 12_000, length: int = 100, seed: int = 1, gz: bool | None = None) -> Path:
    rng = random.Random(seed)
    lines = []
    for i in range(n_reads):
        seq = "".join(rng.choice("ACGT") for _ in range(length))
        lines += [f"@r{i}/1", seq, "+", "I" * length]
    data = ("\n".join(lines) + "\n").encode()
    if gz is None:
        gz = str(path).endswith(".gz")
    path.write_bytes(gzip.compress(data) if gz else data)
    return path


def fake_tool(path: Path, body: str) -> Path:
    path.write_text("#!/usr/bin/env bash\n" + body + "\n")
    path.chmod(0o755)
    return path
