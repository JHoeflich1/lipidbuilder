from pathlib import Path
from typing import List, Tuple

import numpy as np
from lipid import Lipid
import os
import shutil
import json
from datetime import datetime


from openff.toolkit import ForceField, Molecule, Topology
from openff.units import unit


# functions for seting up packmol input file, running and minimizing 

def copy_lipid_files(lipids: List[Lipid], source_dir: Path, dest_dir: Path):
    """Copy each lipid's .pdb and .top file to the working directory."""
    for lipid in lipids:
        for ext in [".pdb", ".top"]:
            src = source_dir / lipid.name / f"{lipid.name}{ext}"
            dst = dest_dir / f"{lipid.name}{ext}"
            if not src.exists():
                raise FileNotFoundError(f"{src} does not exist.")
            shutil.copy(src, dst)


def init_config_file() -> dict:
    """Create and return a JSON metadata dictionary for the build."""
    current_date = datetime.now().isoformat()
    return {
        "experiment": "Bilayer Build",
        "date": current_date,
        "parameters": {}
    }


def save_config(file_path: Path, data: dict):
    with open(file_path, "w") as f:
        json.dump(data, f, indent=4)


def build_packmol_input(
    lipid_names: List[str],
    lipid_counts: List[int],
    solvent_name: str,
    solvent_count: int,
    tolerance: float = 2.0,
    output_name: str = "packmol_input.inp",
    lipid_library_csv: str = "data/available-lipids/PulledLipid.csv",
) -> Tuple[str, str, List[float]]:

    if len(lipid_names) != len(lipid_counts):
        raise ValueError("lipid_names and lipid_counts must be the same length.")

    # Create Lipid objects
    lipids = [Lipid(name, lipid_library_csv) for name in lipid_names]
    solvent = Lipid(solvent_name, lipid_library_csv) if solvent_name else None

    # Copy files to working directory
    cwd = Path.cwd()
    copy_lipid_files(lipids, cwd / "data/available-lipids/lipids_parameterized", cwd)
    if solvent:
        copy_lipid_files([solvent], cwd / "data/available-lipids/lipids_parameterized", cwd)

    # Determine box size
    n = 7  # spacing factor
    total_count = sum(lipid_counts)
    xy = np.sqrt(total_count) * n
    max_z = max([lipid.head_to_tail_distance for lipid in lipids]) * 1.2
    z = 2 * max_z

    if solvent:
        density_water = 1e-21 / 18.01  # mol/nm^3
        water_per_layer = solvent_count / 2
        volume_per_layer = water_per_layer / (density_water * 6.02e23) * 3.5 * 1000
        z += volume_per_layer / (xy ** 2)

    box_dims = [xy, xy, z]

    # Build input
    lines = [
        f"tolerance {tolerance}",
        "filetype pdb",
        f"output packmol_output.pdb",
        ""
    ]

    for lipid, count in zip(lipids, lipid_counts):
        for leaflet, z_sign in [("top", 1), ("bottom", -1)]:
            z_low = 0 if leaflet == "top" else -max_z
            z_high = max_z if leaflet == "top" else 0
            z_center = (z_high + z_low) / 2
            lines.extend([
                f"structure {lipid.name}.pdb",
                f"  number {int(count / 2)}",
                f"  inside box 0. 0. {z_low:.2f} {xy:.2f} {xy:.2f} {z_high:.2f}",
                f"  atoms {lipid.tailgroup_atom_index}",
                f"    below plane 0. 0. 1. {z_center:.2f}",
                f"  end atoms",
                f"  atoms {lipid.headgroup_atom_index}",
                f"    over plane 0. 0. 1. {z_center:.2f}",
                f"  end atoms",
                "end structure",
                ""
            ])

    if solvent:
        for z_start in [max_z, -max_z - 20]:
            lines.extend([
                f"structure {solvent.name}.pdb",
                f"  number {int(solvent_count / 2)}",
                f"  inside box 0. 0. {z_start:.2f} {xy:.2f} {xy:.2f} {z_start + 20:.2f}",
                "end structure",
                ""
            ])

    # Write Packmol input
    packmol_file = cwd / output_name
    with open(packmol_file, "w") as f:
        f.write("\n".join(lines))

    # Save config
    config = init_config_file()
    config["parameters"] = {
        "lipid_names": lipid_names,
        "lipid_counts": lipid_counts,
        "solvent_name": solvent_name,
        "solvent_count": solvent_count,
        "box_dimensions": box_dims,
        "output_name": output_name
    }
    save_config(cwd / "config.json", config)

    return str(packmol_file), "packmol_output.pdb", box_dims
