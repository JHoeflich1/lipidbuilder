import pandas as pd

from .utils import LIPID_LIBRARY_PATH

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

        # Load info from lipid library CSV
        try:
            lipid_library = pd.read_csv(LIPID_LIBRARY_PATH)
        except FileNotFoundError:
            raise FileNotFoundError(f"Lipid library file not found: {LIPID_LIBRARY_PATH}")

        info = lipid_library[lipid_library["Name"] == name]
        if info.empty:
            raise ValueError(f"Lipid '{name}' not found in {LIPID_LIBRARY_PATH}")

        # Required
        self.smiles = info["Smiles String"].values[0]
        self.head_to_tail_distance = float(info["HG/TG distance"].values[0])
        self.headgroup_atom_index = int(info["Headgroup Atom Index"].values[0])
        self.tailgroup_atom_index = int(info["Tailgroup Atom Index"].values[0])
        self.experimental_density = (
            float(info["Experimental Density"].values[0])
            if "Experimental Density" in info.columns and pd.notna(info["Experimental Density"].values[0])
            else None
        )

    def to_dict(self):
        return {
            "name": self.name,
            "smiles": self.smiles,
            "head_to_tail_distance": self.head_to_tail_distance,
            "headgroup_atom_index": self.headgroup_atom_index,
            "tailgroup_atom_index": self.tailgroup_atom_index,
            "experimental_density": self.experimental_density,
        }

    def __repr__(self):
        return f"<Lipid {self.name}: head_to_tail_distance={self.head_to_tail_distance} A>"


def list_available_lipids():
    """
    Lists all lipid names available in the lipid library CSV.

    Returns
    -------
    list of str
        Names of all available lipids.
    """
    if not LIPID_LIBRARY_PATH.exists():
        raise FileNotFoundError(f"Lipid library file not found: {LIPID_LIBRARY_PATH}")

    lipid_library = pd.read_csv(LIPID_LIBRARY_PATH)

    if "Name" not in lipid_library.columns:
        raise ValueError(f"'Name' column not found in {LIPID_LIBRARY_PATH}")

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
