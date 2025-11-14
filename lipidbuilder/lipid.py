from openff.toolkit import Molecule
import pandas as pd
from pathlib import Path
import os

#from .utils import DATA_DIR
from lipidbuilder.utils import DATA_DIR

class Lipid:
    def __init__(self, name: str):
        """
        Represents a lipid and loads its properties from a lipid library CSV.

        Parameters
        ----------
        name : str
            The name of the lipid (must match the "Name" column in the CSV).
        """
        self.name = name

        # Ensure path is a string or Path
        library_path = DATA_DIR / "available-lipids" / "PulledLipid.csv"

        # Load info from lipid library CSV
        try:
            lipid_library = pd.read_csv(library_path)
        except FileNotFoundError:
            raise FileNotFoundError(f"Lipid library file not found: {library_path}")

        info = lipid_library[lipid_library["Name"] == name]
        if info.empty:
            raise ValueError(f"Lipid '{name}' not found in {library_path}")

        # Required
        self.smiles = info["Smiles String"].values[0]
        self.head_to_tail_distance = float(info["HG/TG distance"].values[0])
        self.headgroup_atom_index = int(info["Headgroup Atom Index"].values[0])
        self.tailgroup_atom_index = int(info["Tailgroup Atom Index"].values[0])
        
        # # Optional fields
        # self.volume = float(info["Volume"].values[0]) if "Volume" in info.columns else None
        # self.charge = float(info["Net Charge"].values[0]) if "Net Charge" in info.columns else None


        # Create OpenFF molecule (optional, but nice to cache) 
        # self.molecule = Molecule.from_smiles(self.smiles)
        # self.molecule.name = self.name

    def to_dict(self):
        return {
            "name": self.name,
            "smiles": self.smiles,
            "head_to_tail_distance": self.head_to_tail_distance,
            "volume": self.volume,
            "charge": self.charge,
        }

    def __repr__(self):
        return f"<Lipid {self.name}: Δz={self.head_to_tail_distance} Å, charge={self.charge}>"


def list_available_lipids():
    """
    Lists all lipid names available in the lipid library CSV.

    Returns
    -------
    list of str
        Names of all available lipids.
    """
    library_path = DATA_DIR / "available-lipids" / "PulledLipid.csv"

    if not library_path.exists():
        raise FileNotFoundError(f"Lipid library file not found: {library_path}")

    lipid_library = pd.read_csv(library_path)

    if "Name" not in lipid_library.columns:
        raise ValueError(f"'Name' column not found in {library_path}")

    lipid_names = lipid_library["Name"].dropna().unique().tolist()

    print(f"\n {len(lipid_names)} lipids available:\n" + "\n".join(f" - {n}" for n in lipid_names))

    return lipid_names


if __name__ == "__main__":
    import sys

    # If you ran: python lipid.py list_available_lipids
    if len(sys.argv) > 1 and sys.argv[1] == "list_available_lipids":
        list_available_lipids()
    else:
        print("Usage: python lipid.py list_available_lipids")
