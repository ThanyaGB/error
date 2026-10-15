"""Input-validation failures and tool failures, each mapped to the right error code."""
import gzip

import pytest

from config import settings
from pipeline import assemble, qc
from pipeline.errors import (ASSEMBLY_TIMEOUT, INPUT_CORRUPT, INPUT_EMPTY, INPUT_NOT_FASTQ,
                             INPUT_TOO_FEW_READS, MODEL_MISSING, TOOL_FAILED, PipelineError)
from pipeline.proc import ToolTimeout
from tests.helpers import fake_tool, write_fastq


def code_of(fn, *a):
    with pytest.raises(PipelineError) as ei:
        fn(*a)
    return ei.value.code


def test_valid_fastq_and_gz(tmp):
    assert qc.validate_fastq(write_fastq(tmp / "a.fastq"))["n_reads"] == 12_000
    assert qc.validate_fastq(write_fastq(tmp / "a.fastq.gz"))["n_reads"] == 12_000


def test_not_fastq(tmp):
    f = tmp / "x.fastq"
    f.write_text(">seq1\nACGT\n")
    assert code_of(qc.validate_fastq, f) == INPUT_NOT_FASTQ
    f.write_bytes(b"\x89PNG\r\n\x1a\n" + b"\x00" * 50)
    assert code_of(qc.validate_fastq, f) == INPUT_NOT_FASTQ


def test_empty(tmp):
    f = tmp / "e.fastq"
    f.write_bytes(b"")
    assert code_of(qc.validate_fastq, f) == INPUT_EMPTY
    g = tmp / "e.fastq.gz"
    g.write_bytes(gzip.compress(b""))
    assert code_of(qc.validate_fastq, g) == INPUT_EMPTY


def test_truncated_gzip(tmp):
    good = write_fastq(tmp / "t.fastq.gz")
    data = good.read_bytes()
    good.write_bytes(data[: len(data) // 2])
    assert code_of(qc.validate_fastq, good) == INPUT_CORRUPT


def test_corrupt_record_and_cut_file(tmp):
    f = write_fastq(tmp / "c.fastq")
    f.write_bytes(f.read_bytes()[:-30])                   # last record cut off
    assert code_of(qc.validate_fastq, f) == INPUT_CORRUPT


def test_too_few_reads(tmp):
    assert code_of(qc.validate_fastq, write_fastq(tmp / "s.fastq", n_reads=500)) == INPUT_TOO_FEW_READS


def test_read_length_out_of_range(tmp):
    assert code_of(qc.validate_fastq, write_fastq(tmp / "l.fastq", length=20)) == INPUT_CORRUPT
    assert code_of(qc.validate_fastq, write_fastq(tmp / "m.fastq", length=400)) == INPUT_CORRUPT


def test_pair_mismatch(tmp):
    a = write_fastq(tmp / "a.fastq", n_reads=12_000)
    b = write_fastq(tmp / "b.fastq", n_reads=13_000)
    assert code_of(qc.validate_pair, a, b) == INPUT_CORRUPT


def test_extension_check():
    assert qc.has_fastq_extension("R1.fq.gz") and qc.has_fastq_extension("a.FASTQ")
    assert not qc.has_fastq_extension("reads.txt") and not qc.has_fastq_extension("reads.fasta")


def test_assembly_timeout(tmp, monkeypatch):
    monkeypatch.setattr(settings, "SPADES", str(fake_tool(tmp / "slow.sh", "sleep 30")))
    monkeypatch.setattr(settings, "ASSEMBLY_TIMEOUT_S", 1)
    assert code_of(assemble.run_assembly, tmp / "r1", tmp / "r2", tmp / "out") == ASSEMBLY_TIMEOUT


def test_tool_failed_on_nonzero_exit(tmp, monkeypatch):
    monkeypatch.setattr(settings, "SPADES", str(fake_tool(tmp / "bad.sh", "echo boom >&2; exit 3")))
    assert code_of(assemble.run_assembly, tmp / "r1", tmp / "r2", tmp / "out") == TOOL_FAILED


def test_tool_failed_when_tool_cannot_start(tmp, monkeypatch):
    monkeypatch.setattr(settings, "SPADES", str(tmp / "does-not-exist"))
    assert code_of(assemble.run_assembly, tmp / "r1", tmp / "r2", tmp / "out") == TOOL_FAILED


def test_model_missing(monkeypatch, tmp):
    from pipeline import predict
    monkeypatch.setattr(settings, "ARTIFACT_DIR", tmp)
    assert code_of(predict.predict_all, {"values": {}, "sources": {}}) == MODEL_MISSING


def test_assembly_stats(tmp):
    f = tmp / "c.fa"
    f.write_text(">a\n" + "A" * 100 + "\n>b\n" + "A" * 300 + "\n>c\n" + "A" * 600 + "\n")
    assert assemble.contig_stats(f) == {"n_contigs": 3, "total_length": 1000, "n50": 600}


def test_tool_timeout_exception_type(tmp):
    from pipeline.proc import run_tool
    with pytest.raises(ToolTimeout):
        run_tool([str(fake_tool(tmp / "s.sh", "sleep 30"))], timeout=1)
