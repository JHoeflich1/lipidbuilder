from pathlib import Path

PACKAGE_DIR = Path(__file__).resolve().parent
PROJECT_DIR = PACKAGE_DIR.parent

DATA_DIR = PACKAGE_DIR / "data"
RESULTS_DIR = PACKAGE_DIR / "results"
DEFAULT_RESULTS_DIR = RESULTS_DIR / "lipid_system"

AVAILABLE_LIPIDS_DIR = DATA_DIR / "available-lipids"
LIPID_LIBRARY_PATH = AVAILABLE_LIPIDS_DIR / "PulledLipid.csv"
LIPID_PDB_DIR = AVAILABLE_LIPIDS_DIR / "lipid_pdbs"
SOLVENT_PDB_DIR = AVAILABLE_LIPIDS_DIR / "solvent_pdbs"
FORCEFIELD_DIR = DATA_DIR / "forcefields"

# Backward-compatible aliases for older scripts.
ROOT_DIR = PACKAGE_DIR
LIP_DIR = AVAILABLE_LIPIDS_DIR

SCRATCH_DIR = AVAILABLE_LIPIDS_DIR / "scripts" / "scratch"


def resolve_force_field_path(force_field_file: str | Path) -> Path:
    """Resolve a force-field file from an absolute path or packaged data."""
    path = Path(force_field_file).expanduser()
    if path.is_absolute():
        return path
    return FORCEFIELD_DIR / path
