from openff.toolkit import Molecule
import pandas as pd
import os

class Lipid:
    def __init__(self, name: str, library_path: str = "data/available-lipids/PulledLipid.csv"):
        self.name = name

        # Load info from lipid library CSV
        lipid_library = pd.read_csv(library_path)
        info = lipid_library[lipid_library["Name"] == name]
        if info.empty:
            raise ValueError(f"Lipid '{name}' not found in {library_path}")

        self.smiles = info["Lipid Smiles String"].values[0]
        self.head_to_tail_distance = float(info["HG/TG distance"].values[0])
        self.headgroup_atom_index = int(info["Headgroup Atom Index"].values[0])
        self.tailgroup_atom_index = int(info["Tailgroup Atom Index"].values[0])
        self.volume = float(info["Approximate Volume"].values[0]) if "Approximate Volume" in info else None
        self.charge = float(info["Net Charge"].values[0]) if "Net Charge" in info else None

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
