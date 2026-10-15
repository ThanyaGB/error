"""Central configuration. Every value can be overridden by an environment variable."""
import os
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

# ---------------------------------------------------------------- MongoDB
MONGO_URI = os.environ.get("AMR_MONGO_URI", "mongodb://127.0.0.1:27017")
DB_NAME = os.environ.get("AMR_DB_NAME", "amr_firewall")
MONGO_TIMEOUT_MS = int(os.environ.get("AMR_MONGO_TIMEOUT_MS", "1500"))

# ---------------------------------------------------------------- Paths
# Heavy files (uploads, SPAdes scratch, contigs) live outside the repo, on the
# WSL ext4 disk, because /mnt/c is slow for the assembler.
WORK_DIR = Path(os.environ.get("AMR_WORK_DIR", str(Path.home() / "amr_work")))
DATA_DIR = Path(os.environ.get("AMR_DATA_DIR", str(ROOT / "data")))
SPLITS_CSV = DATA_DIR / "splits.csv"
TRAIN_CACHE_DIR = Path(os.environ.get("AMR_TRAIN_CACHE", str(WORK_DIR / "training")))
ARTIFACT_DIR = Path(os.environ.get("AMR_ARTIFACT_DIR", str(ROOT / "model" / "artifacts")))
WEB_DIR = ROOT / "web"

# ---------------------------------------------------------------- Limits
MAX_UPLOAD_BYTES = 1024 ** 3          # 1 GB per file
MIN_READS = 10_000
MIN_READ_LEN = 50
MAX_READ_LEN = 300
ASSEMBLY_TIMEOUT_S = int(os.environ.get("AMR_ASSEMBLY_TIMEOUT_S", str(40 * 60)))
DEEPARG_TIMEOUT_S = int(os.environ.get("AMR_DEEPARG_TIMEOUT_S", str(20 * 60)))
TOOL_TIMEOUT_S = int(os.environ.get("AMR_TOOL_TIMEOUT_S", str(30 * 60)))
THREADS = int(os.environ.get("AMR_THREADS", "4"))
SPADES_MEM_GB = int(os.environ.get("AMR_SPADES_MEM_GB", "6"))

# ---------------------------------------------------------------- Calls
R_THRESHOLD = 0.7
S_THRESHOLD = 0.3

# ---------------------------------------------------------------- Domain
STEPS = ["qc", "assemble", "amr_calls", "deeparg_calls", "features", "predict"]

# drug -> drug class shown in the UI
DRUGS = {
    "ciprofloxacin": "fluoroquinolone",
    "ceftazidime": "cephalosporin",
    "meropenem": "carbapenem",
    "gentamicin": "aminoglycoside",
    "tetracycline": "tetracycline",
    "trimethoprim-sulfamethoxazole": "folate pathway inhibitor",
}

# species key (lowercase, underscore) -> AMRFinderPlus --organism value
SUPPORTED_SPECIES = {
    "escherichia_coli": "Escherichia",
    "klebsiella_pneumoniae": "Klebsiella_pneumoniae",
}

# ---------------------------------------------------------------- External tools
def _tool(env_name: str, binary: str) -> str:
    """Resolve a tool: explicit env var, else next to this Python, else PATH."""
    if os.environ.get(env_name):
        return os.environ[env_name]
    beside = Path(sys.executable).parent / binary
    if beside.exists():
        return str(beside)
    return shutil.which(binary) or binary


FASTP = _tool("AMR_FASTP", "fastp")
SPADES = _tool("AMR_SPADES", "spades.py")
AMRFINDER = _tool("AMR_AMRFINDER", "amrfinder")
DEEPARG = _tool("AMR_DEEPARG", "deeparg")
PRODIGAL = _tool("AMR_PRODIGAL", "prodigal")
DIAMOND = _tool("AMR_DIAMOND", "diamond")
AMRFINDER_DB = os.environ.get("AMR_AMRFINDER_DB")  # None -> tool default
DEEPARG_MODEL = os.environ.get("AMR_DEEPARG_MODEL", "LS")


def tool_env() -> dict:
    """Environment for child processes: the tool dir first on PATH, HF offline."""
    env = dict(os.environ)
    env["PATH"] = str(Path(sys.executable).parent) + os.pathsep + env.get("PATH", "")
    env["HF_HUB_OFFLINE"] = "1"           # the demo must run without internet
    env["TRANSFORMERS_OFFLINE"] = "1"
    env.setdefault("CONDA_PREFIX", str(Path(sys.executable).parent.parent))
    return env
