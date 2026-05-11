import json
import time
import shutil
import subprocess
from datetime import datetime
from pathlib import Path
from typing import List, Optional, Tuple

import mdtraj
import numpy as np
import pandas as pd

from openff.toolkit import Molecule, Topology, ForceField
from openff.interchange import Interchange
from openff.toolkit.utils.nagl_wrapper import NAGLToolkitWrapper
from openff.units import unit
from .lipid import Lipid
from .utils import (
    FORCEFIELD_DIR,
    ION_PDB_DIR,
    LIPID_LIBRARY_PATH,
    LIPID_PDB_DIR,
    SOLVENT_PDB_DIR,
    resolve_force_field_path,
)

# functions for seting up packmol input file, running and minimizing 

def copy_lipid_files(lipids: List[Lipid], source_dir: Path = None, dest_dir: Path = None):
    """Copy each lipid's .pdb and .top file to the working directory."""
    if source_dir is None:
        source_dir = LIPID_PDB_DIR
    if dest_dir is None:
        dest_dir = Path.cwd()

    for lipid in lipids:
        for ext in [".pdb"]:
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
    force_field_file: str,
    charge_model: str,
    hmr: bool,
    ion_type: Optional[str] = None,
    ion_count: int = 0,
    tolerance: float = 2.0,
    input_name: str = "packmol_input.inp",
    ) -> Tuple[str, List[float]]:

    lipids = [Lipid(name) for name in lipid_names]

    cwd = Path.cwd()

    # Copy lipid files from absolute repo data dir into current working directory
    copy_lipid_files(lipids, source_dir=LIPID_PDB_DIR, dest_dir=cwd)

    # Copy solvent files from absolute repo data dir into cwd
    if solvent_name:
        for ext in [".pdb"]:
            src = SOLVENT_PDB_DIR / f"{solvent_name}{ext}"
            dst = cwd / f"{solvent_name}{ext}"
            if not src.exists():
                raise FileNotFoundError(f"Solvent file {src} does not exist.")
            shutil.copy(src, dst)

    if ion_type is not None:
        ion_type = ion_type.strip().capitalize()
        if ion_type not in {"Na", "Cl"}:
            raise ValueError("ion_type must be either 'Na' or 'Cl'.")
    if ion_count < 0:
        raise ValueError("ion_count must be non-negative.")
    ion_count = int(ion_count)
    ion_file_stem = ion_type.lower() if ion_type else None
    if ion_type and ion_count > 0:
        src = ION_PDB_DIR / f"{ion_file_stem}.pdb"
        dst = cwd / f"{ion_file_stem}.pdb"
        if not src.exists():
            raise FileNotFoundError(f"Ion file {src} does not exist.")
        shutil.copy(src, dst)



    area_per_lipid = 10  # not sure how I came on this value, but it works
    N_top = sum(lipid_counts[::2])
    N_bottom = sum(lipid_counts[1::2])

    # Take the max for xy (largest leaflet sets the box)
    N_leaflet_max = max(N_top, N_bottom)
    xy = np.sqrt(N_leaflet_max) * area_per_lipid # ideally this is sqrt(N * area) but cant figure out the right area per lipid to use, so just fudging with a constant for now.


    output_name = "packmol_output.pdb"
    # Build input
    lines = [
        f"tolerance {tolerance}",
        "filetype pdb",
        f"output {output_name}",
        ""
    ]

    #Could first determine the height of both leaflets so that would be my maz z for each leaflet
    max_z_top = 0.0
    max_z_bottom = 0.0

    for i, lipid in enumerate(lipids):
        top_count = lipid_counts[2*i]      # top leaflet: the user-provided number of lipids for the top leaflet, where i is the index of the lipid. Top are even 
        bottom_count = lipid_counts[2*i+1] # bottom leaflet: user-provided number of lipids for the bottom leaflet. bottom are odd indicies
        # Loop over top and bottom leaflets for this lipid
        for leaflet, count in [("top", top_count), ("bottom", bottom_count)]:
                if count <= 0:
                    continue  # skip this block entirely if no lipids, packmol does not accept 0 lipids

                # Begin Packmol structure block for this lipid/leaflet
                if leaflet == "top":
                    # top leaflet: tails below plane, heads above plane
                    lines.extend([
                        f"structure {lipid.name}.pdb",
                        f"  number {count}",
                        f"  inside box 0. 0. 0. {xy:.2f} {xy:.2f} {lipid.head_to_tail_distance * 1.1:.2f}",
                        f"  atoms {lipid.tailgroup_atom_index}",
                        f"    below plane 0. 0. 1. {lipid.head_to_tail_distance * 0.1:.2f}",
                        f"  end atoms",
                        f"  atoms {lipid.headgroup_atom_index}",
                        f"    over plane 0. 0. 1. {lipid.head_to_tail_distance * 1.1 - lipid.head_to_tail_distance * 0.1:.2f}",  
                        f"  end atoms",
                        "end structure",
                        ""
                    ])
                else:
                    # bottom leaflet
                    lines.extend([
                        f"structure {lipid.name}.pdb",
                        f"  number {count}",
                        f"  inside box 0. 0. -{lipid.head_to_tail_distance * 1.1:.2f} {xy:.2f} {xy:.2f} 0.",
                        f"  atoms {lipid.tailgroup_atom_index}",
                        f"    over plane 0. 0. 1. -{lipid.head_to_tail_distance * 0.1:.2f}",
                        f"  end atoms",
                        f"  atoms {lipid.headgroup_atom_index}",
                        f"    below plane 0. 0. 1. -{lipid.head_to_tail_distance * 1.1 - lipid.head_to_tail_distance * 0.1:.2f}",
                        f"  end atoms",
                        "end structure",
                        ""
                    ])

    max_z_top = max(
        (lipid.head_to_tail_distance * 1.1 for i, lipid in enumerate(lipids) if lipid_counts[2*i] > 0),
        default=0.0,
    )
    max_z_bottom = max(
        (lipid.head_to_tail_distance * 1.1 for i, lipid in enumerate(lipids) if lipid_counts[2*i+1] > 0),
        default=0.0,
    )
    solvent_layer_thickness = 0.0
    
    if solvent_name and solvent_count > 0:
        density_water = 1e-21 / 18.01  # mol/nm^3  experimental density of water
        top_solvent_count = solvent_count // 2 + solvent_count % 2
        bottom_solvent_count = solvent_count // 2
        water_per_layer = solvent_count / 2
        fudge_factor_packing = 2
        #estimate volume needed to contain a number of water molecuels per layer
        volume_per_layer = water_per_layer / (density_water * 6.02e23) * fudge_factor_packing * 1000 # 3.5 is a fudge factor to give packmol room for packing
        solvent_layer_thickness = volume_per_layer / (xy ** 2)  # adaptive thickness based on number of waters

        # Place water above top leaflet
        z_start_top = max_z_top
        lines.extend([
            f"structure {solvent_name}.pdb",
            f"  number {top_solvent_count}",
            f"  inside box 0. 0. {z_start_top * 1 :.2f} {xy:.2f} {xy:.2f} {z_start_top * 1 + solvent_layer_thickness:.2f}",
            "end structure",
            ""
        ])
        # Place water below bottom leaflet
        z_start_bottom = -max_z_bottom - solvent_layer_thickness
        lines.extend([
            f"structure {solvent_name}.pdb",
            f"  number {bottom_solvent_count}",
            f"  inside box 0. 0. {z_start_bottom * 1 :.2f} {xy:.2f} {xy:.2f} {z_start_bottom* 1 + solvent_layer_thickness:.2f}",
            "end structure",
            ""
        ])

    if ion_type and ion_count > 0:
        top_ion_count = ion_count // 2 + ion_count % 2
        bottom_ion_count = ion_count // 2

        if solvent_layer_thickness > 0:
            z_top_min = max_z_top
            z_top_max = max_z_top + solvent_layer_thickness
            z_bottom_min = -max_z_bottom - solvent_layer_thickness
            z_bottom_max = -max_z_bottom
        else:
            padding = 10.0
            z_top_min = max_z_top
            z_top_max = max_z_top + padding
            z_bottom_min = -max_z_bottom - padding
            z_bottom_max = -max_z_bottom
            solvent_layer_thickness = padding

        for count, z_min, z_max in [
            (top_ion_count, z_top_min, z_top_max),
            (bottom_ion_count, z_bottom_min, z_bottom_max),
        ]:
            if count <= 0:
                continue
            lines.extend([
                f"structure {ion_file_stem}.pdb",
                f"  number {count}",
                f"  inside box 0. 0. {z_min:.2f} {xy:.2f} {xy:.2f} {z_max:.2f}",
                "end structure",
                ""
            ])

    # What is our final z dimension for box_dims?

    # Write Packmol input
    packmol_file = cwd / input_name
    with open(packmol_file, "w") as f:
        f.write("\n".join(lines))

    z = max_z_top + max_z_bottom + 2 * solvent_layer_thickness
    box_dims = [xy, xy, z]

    # Save config
    config = init_config_file()
    config["parameters"] = {
        "lipid_names": lipid_names,
        "lipid_counts": lipid_counts,
        "solvent_name": solvent_name,
        "solvent_count": solvent_count,
        "ion_type": ion_type,
        "ion_count": ion_count,
        "box_dimensions": box_dims,
        "input_name": input_name,
        "output_name": output_name,
        "force_field_file": force_field_file,
        "charge_model": charge_model,
        "hmr": hmr
    }
    save_config(cwd / "config.json", config)

    

    return str(packmol_file), box_dims


def run_packmol(
    packmol_input_file: str,
    parameterize: bool,
    force_field_file: str,
    hmr: bool,
    charge_model: str,
    config_path: str = "config.json"
    ):
    """Run Packmol to generate coordinates and parameterize the system with optional HMR."""
    print("Running Packmol...")
    start_time = time.time()

    log_file = "packmol_output.log"
    subprocess.run(f"packmol < {packmol_input_file}", shell=True, check=True,
                   stdout=open(log_file, "w"), stderr=subprocess.STDOUT)

    subprocess.run("tail -n 28 packmol_output.log", shell=True)

    if not parameterize:
        print("System packing complete. Parameterization skipped.")
        return

    print("Starting parameterization...")
    with open(config_path, 'r') as f:
        config = json.load(f)

    box_dims = config['parameters']['box_dimensions']
    lipid_names = config['parameters']['lipid_names']
    lipid_counts = config['parameters']['lipid_counts']
    solvent_name = config['parameters']['solvent_name']
    solvent_count = config['parameters']['solvent_count']
    ion_type = config['parameters'].get('ion_type')
    ion_count = int(config['parameters'].get('ion_count', 0) or 0)
    packmol_output = config['parameters']['output_name']
    # Add charge model and froce field file and HMR 


    # Paths
    lipid_library_path = LIPID_LIBRARY_PATH
    force_field_path = resolve_force_field_path(force_field_file)

    if not force_field_path.exists():
        packaged_force_fields = ", ".join(path.name for path in sorted(FORCEFIELD_DIR.glob("*.offxml")))
        raise FileNotFoundError(
            f"Force field file not found: {force_field_path}. "
            f"Packaged force fields: {packaged_force_fields or 'none'}"
        )

    lipid_library = pd.read_csv(lipid_library_path)
    forcefield = ForceField(force_field_path)

    lipid_molecules = []
    for name in lipid_names:
        try:
            smiles = lipid_library[lipid_library["Name"] == name]["Smiles String"].values[0]
        except IndexError:
            raise ValueError(f"Lipid '{name}' not found in lipid library: {lipid_library_path}")
        
        molecule = Molecule.from_smiles(smiles, allow_undefined_stereo=True)
        if charge_model.endswith(".pt"):
            # Assume a NAGL model
            molecule.assign_partial_charges(
                charge_model,
                toolkit_registry=NAGLToolkitWrapper()
            )
        elif charge_model.lower() == "am1bcc":
            molecule.assign_partial_charges("am1bcc")
        elif charge_model.lower() == "gasteiger":
            molecule.assign_partial_charges("gasteiger")
        elif charge_model.lower() == "am1bccelf10":
            molecule.assign_partial_charges("am1bccelf10")
        else:
            raise ValueError(f"Unsupported charge model: {charge_model}")
        
        molecule.name = name
        for atom in molecule.atoms:
            atom.metadata["residue_name"] = name
        molecule.generate_unique_atom_names()
        lipid_molecules.append(molecule)

    all_lipids = [] # all_lipids is a flat list of Molecule objects representing every single lipid molecule that will go into the system. Used to build the topology. 
    for lipid_index, lipid in enumerate(lipid_molecules):
        top_count = lipid_counts[2*lipid_index]       # top leaflet
        bottom_count = lipid_counts[2*lipid_index+1] # bottom leaflet
        all_lipids.extend([lipid]*top_count)
        all_lipids.extend([lipid]*bottom_count)

    solvent_mol = None
    if solvent_name and solvent_count > 0:
        solvent_mol = Molecule.from_file(f"{solvent_name}.pdb")
        for atom in solvent_mol.atoms:
            atom.metadata["residue_name"] = solvent_name.upper()

    ion_mol = None
    if ion_type and ion_count > 0:
        ion_smiles = {"Na": "[Na+]", "Cl": "[Cl-]"}[ion_type]
        ion_mol = Molecule.from_smiles(ion_smiles)
        ion_mol.name = ion_type
        ion_charge = {"Na": 1.0, "Cl": -1.0}[ion_type]
        ion_mol.partial_charges = np.array([ion_charge]) * unit.elementary_charge
        for atom in ion_mol.atoms:
            atom.name = ion_type
            atom.metadata["residue_name"] = ion_type.upper()

    all_molecules = all_lipids
    if solvent_mol is not None:
        all_molecules.extend([solvent_mol] * solvent_count)
    if ion_mol is not None:
        all_molecules.extend([ion_mol] * ion_count)

    topology = Topology.from_molecules(all_molecules)
    traj = mdtraj.load(packmol_output)
    topology.set_positions(traj.xyz[0] * unit.nanometer)
    topology.box_vectors = np.array(box_dims) * 0.1 * unit.nanometer # 0.1 to convert to nm from Angstrom 

    interchange = Interchange.from_smirnoff(
        force_field=forcefield,
        topology=topology,
        charge_from_molecules=lipid_molecules + ([ion_mol] if ion_mol is not None else []),
    )

    if hmr:
        interchange.to_gromacs(prefix="bilayer", decimal=3, hydrogen_mass=3.0)
        print("System saved with HMR as bilayer.top / bilayer.gro.")
    else:
        interchange.to_gromacs(prefix="bilayer")
        print("System saved without HMR as bilayer.top / bilayer.gro.")

    elapsed = time.time() - start_time
    print(f"Parameterization complete in {elapsed:.2f} seconds.")
