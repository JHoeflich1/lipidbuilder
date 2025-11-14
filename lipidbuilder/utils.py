from pathlib import Path

# Root directory of your lipidbuilder project (folder containing utils.py)
ROOT_DIR = Path(__file__).resolve().parent

DATA_DIR = ROOT_DIR / "data"

LIP_DIR = ROOT_DIR / "data" / "available-lipids"

# Path to the scratch directory
SCRATCH_DIR = LIP_DIR / "scripts" / "scratch"

# Path to the lipid_pdbs directory
LIPID_PDB_DIR = LIP_DIR / "lipid_pdbs"

