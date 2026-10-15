"""The only error codes the system may emit, with plain-language messages.

Messages never contain stack traces or file paths.
"""

INPUT_NOT_FASTQ = "INPUT_NOT_FASTQ"
INPUT_EMPTY = "INPUT_EMPTY"
INPUT_CORRUPT = "INPUT_CORRUPT"
INPUT_TOO_FEW_READS = "INPUT_TOO_FEW_READS"
UPLOAD_TOO_LARGE = "UPLOAD_TOO_LARGE"
ASSEMBLY_TIMEOUT = "ASSEMBLY_TIMEOUT"
TOOL_FAILED = "TOOL_FAILED"
MODEL_MISSING = "MODEL_MISSING"
DB_UNAVAILABLE = "DB_UNAVAILABLE"
NOT_FOUND = "NOT_FOUND"

MESSAGES = {
    INPUT_NOT_FASTQ: "One of the files is not a FASTQ file. Upload paired .fastq, .fq, .fastq.gz or .fq.gz files.",
    INPUT_EMPTY: "One of the files is empty. Check the file and upload it again.",
    INPUT_CORRUPT: "One of the files is damaged or cut off and could not be read. Re-export or re-download it and try again.",
    INPUT_TOO_FEW_READS: "There are too few reads to assemble a genome. At least 10,000 reads per file are needed.",
    UPLOAD_TOO_LARGE: "One of the files is larger than 1 GB. Downsample the reads or upload a smaller file.",
    ASSEMBLY_TIMEOUT: "Genome assembly took longer than 40 minutes and was stopped. Try a smaller (downsampled) read set.",
    TOOL_FAILED: "An analysis tool stopped with an error, so the run was halted. Check that the reads are from a bacterial isolate and try again.",
    MODEL_MISSING: "A required model is not installed on this server, so the run was halted.",
    DB_UNAVAILABLE: "The database is not reachable right now. Start MongoDB and try again.",
    NOT_FOUND: "No isolate with that id was found.",
}

HTTP_STATUS = {
    INPUT_NOT_FASTQ: 400,
    INPUT_EMPTY: 400,
    INPUT_CORRUPT: 400,
    INPUT_TOO_FEW_READS: 400,
    UPLOAD_TOO_LARGE: 413,
    ASSEMBLY_TIMEOUT: 500,
    TOOL_FAILED: 500,
    MODEL_MISSING: 500,
    DB_UNAVAILABLE: 503,
    NOT_FOUND: 404,
}

ALL_CODES = tuple(MESSAGES)


class PipelineError(Exception):
    """Raised anywhere in the system. Carries a code from ALL_CODES only."""

    def __init__(self, code: str, detail: str = "", message: str | None = None):
        if code not in MESSAGES:
            raise ValueError(f"unknown error code: {code}")
        super().__init__(code)
        self.code = code
        self.message = message or MESSAGES[code]
        self.detail = detail  # for server logs only, never sent to clients


def payload(code: str, message: str | None = None) -> dict:
    """The JSON body used by every error response."""
    return {"error_code": code, "error_message": message or MESSAGES[code]}
