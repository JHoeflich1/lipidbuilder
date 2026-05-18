import json
import time
import shutil
import subprocess
import warnings
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

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
    LIPID_LIBRARY_PATH,
    LIPID_PDB_DIR,
    SOLVENT_PDB_DIR,
    resolve_force_field_path,
)

# functions for seting up packmol input file, running and minimizing 

IONS_MDP = """; ions.mdp - used as input into grompp to generate ions.tpr
; Parameters describing what to do, when to stop and what to save
integrator  = steep         ; Algorithm (steep = steepest descent minimization)
emtol       = 1000.0        ; Stop minimization when the maximum force < 1000.0 kJ/mol/nm
emstep      = 0.01          ; Minimization step size
nsteps      = 50000         ; Maximum number of (minimization) steps to perform

; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
nstlist         = 1         ; Frequency to update the neighbor list and long range forces
cutoff-scheme   = Verlet    ; Buffered neighbor searching
ns_type         = grid      ; Method to determine neighbor list (simple, grid)
rlist           = 1.2       ; Cut-off for making neighbor list (short range forces)
coulombtype     = cutoff    ; Treatment of long range electrostatic interactions
rcoulomb        = 1.2       ; Short-range electrostatic cut-off
rvdw            = 1.2       ; Short-range Van der Waals cut-off
pbc             = xyz       ; Periodic Boundary Conditions in all 3 dimensions
"""


def write_ions_mdp(path: str | Path = "ions.mdp") -> Path:
    """Write the GROMACS MDP file used to generate ions.tpr for genion."""
    path = Path(path)
    path.write_text(IONS_MDP)
    return path

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


def _format_charge(charge: float) -> str:
    if abs(charge - round(charge)) < 1e-6:
        return f"{int(round(charge)):+d}"
    return f"{charge:+.3f}"


def estimate_bilayer_charge(
    lipid_names: List[str],
    lipid_counts: List[int],
    lipid_library_path: Path = LIPID_LIBRARY_PATH,
) -> float:
    """Estimate the bilayer formal charge from lipid SMILES."""
    lipid_library = pd.read_csv(lipid_library_path)
    total_charge = 0.0
    for lipid_index, name in enumerate(lipid_names):
        try:
            smiles = lipid_library[lipid_library["Name"] == name]["Smiles String"].values[0]
        except IndexError:
            raise ValueError(f"Lipid '{name}' not found in lipid library: {lipid_library_path}")

        molecule = Molecule.from_smiles(smiles, allow_undefined_stereo=True)
        lipid_charge = molecule.total_charge.m_as(unit.elementary_charge)
        total_count = lipid_counts[2*lipid_index] + lipid_counts[2*lipid_index+1]
        total_charge += lipid_charge * total_count

    return total_charge


def estimate_neutralizing_ion_counts(
    lipid_names: List[str],
    lipid_counts: List[int],
    lipid_library_path: Path = LIPID_LIBRARY_PATH,
) -> Dict[str, int]:
    """Estimate neutralizing NA/CL counts from lipid formal charges."""
    total_charge = estimate_bilayer_charge(
        lipid_names,
        lipid_counts,
        lipid_library_path=lipid_library_path,
    )
    rounded_charge = int(round(total_charge))
    return {
        "NA": max(-rounded_charge, 0),
        "CL": max(rounded_charge, 0),
    }


def build_packmol_input(
    lipid_names: List[str],
    lipid_counts: List[int],
    solvent_name: str,
    solvent_count: int,
    force_field_file: str,
    charge_model: str,
    hmr: bool,
    neutralize_ions: bool = False,
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
    
    total_lipid_count = sum(lipid_counts)
    target_solvent_count = int(solvent_count)
    estimated_bilayer_charge = 0.0
    estimated_ion_counts = {"NA": 0, "CL": 0}
    if neutralize_ions:
        estimated_bilayer_charge = estimate_bilayer_charge(lipid_names, lipid_counts)
        estimated_ion_counts = estimate_neutralizing_ion_counts(lipid_names, lipid_counts)
    estimated_neutralizing_ions = sum(estimated_ion_counts.values())
    packmol_solvent_count = target_solvent_count + estimated_neutralizing_ions

    if neutralize_ions:
        print(
            "Ion neutralization requested: "
            f"net formal charge of bilayer is {_format_charge(estimated_bilayer_charge)} e."
        )
        print(
            "Estimated neutralizing ions: "
            f"NA={estimated_ion_counts['NA']}, CL={estimated_ion_counts['CL']}."
        )
        print(
            "Hydration target will be preserved by packing extra waters: "
            f"target={target_solvent_count}, extra={estimated_neutralizing_ions}, "
            f"Packmol waters={packmol_solvent_count}."
        )

    if solvent_name and packmol_solvent_count > 0:
        density_water = 1e-21 / 18.01  # mol/nm^3  experimental density of water
        top_solvent_count = packmol_solvent_count // 2 + packmol_solvent_count % 2
        bottom_solvent_count = packmol_solvent_count // 2
        water_per_layer = packmol_solvent_count / 2
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

    # What is our final z dimension for box_dims?

    # Write Packmol input
    packmol_file = cwd / input_name
    with open(packmol_file, "w") as f:
        f.write("\n".join(lines))

    z = max_z_top + max_z_bottom + 2 * solvent_layer_thickness
    box_dims = [xy, xy, z]
    target_hydration_level = (
        target_solvent_count / total_lipid_count
        if total_lipid_count > 0
        else 0.0
    )
    packmol_hydration_level = (
        packmol_solvent_count / total_lipid_count
        if total_lipid_count > 0
        else 0.0
    )

    # Save config
    config = init_config_file()
    config["parameters"] = {
        "lipid_names": lipid_names,
        "lipid_counts": lipid_counts,
        "solvent_name": solvent_name,
        "solvent_count": packmol_solvent_count,
        "target_solvent_count": target_solvent_count,
        "packmol_solvent_count": packmol_solvent_count,
        "target_hydration_level": target_hydration_level,
        "packmol_hydration_level": packmol_hydration_level,
        "hydration_level_basis": "final_target",
        "estimated_bilayer_charge": estimated_bilayer_charge,
        "estimated_neutralizing_ion_counts": estimated_ion_counts,
        "estimated_neutralizing_ions": estimated_neutralizing_ions,
        "neutralize_ions": neutralize_ions,
        "positive_ion_name": "NA",
        "negative_ion_name": "CL",
        "box_dimensions": box_dims,
        "input_name": input_name,
        "output_name": output_name,
        "force_field_file": force_field_file,
        "charge_model": charge_model,
        "hmr": hmr
    }
    save_config(cwd / "config.json", config)

    

    return str(packmol_file), box_dims


def _make_ion_molecule(ion_name: str) -> Molecule:
    smiles_by_name = {"NA": "[Na+]", "Na": "[Na+]", "CL": "[Cl-]", "Cl": "[Cl-]"}
    charge_by_name = {"NA": 1.0, "Na": 1.0, "CL": -1.0, "Cl": -1.0}
    if ion_name not in smiles_by_name:
        raise ValueError(f"Unsupported ion name: {ion_name}")
    molecule = Molecule.from_smiles(smiles_by_name[ion_name])
    molecule.name = ion_name.upper()
    molecule.partial_charges = np.array([charge_by_name[ion_name]]) * unit.elementary_charge
    for atom in molecule.atoms:
        atom.name = ion_name.upper()
        atom.metadata["residue_name"] = ion_name.upper()
    return molecule


def _gro_residue_names(gro_path: str | Path) -> List[str]:
    """Return contiguous residue names from a GRO file in coordinate order."""
    residue_names = []
    previous_residue = None
    lines = Path(gro_path).read_text().splitlines()
    for line in lines[2:-1]:
        residue_id = line[:5].strip()
        residue_name = line[5:10].strip()
        residue_key = (residue_id, residue_name)
        if residue_key != previous_residue:
            residue_names.append(residue_name)
            previous_residue = residue_key
    return residue_names


def _gro_residue_counts(gro_path: str | Path) -> Dict[str, int]:
    counts: Dict[str, int] = {}
    for residue_name in _gro_residue_names(gro_path):
        counts[residue_name] = counts.get(residue_name, 0) + 1
    return counts


def _build_interchange(
    coordinate_file: str,
    output_prefix: str,
    box_dims: List[float],
    lipid_names: List[str],
    lipid_counts: List[int],
    solvent_name: str,
    solvent_count: int,
    forcefield: ForceField,
    lipid_molecules: List[Molecule],
    charge_model: str,
    hmr: bool,
    ion_counts: Optional[Dict[str, int]] = None,
    preserve_coordinate_order: bool = False,
) -> None:
    all_lipids = []
    for lipid_index, lipid in enumerate(lipid_molecules):
        top_count = lipid_counts[2*lipid_index]
        bottom_count = lipid_counts[2*lipid_index+1]
        all_lipids.extend([lipid]*top_count)
        all_lipids.extend([lipid]*bottom_count)

    solvent_mol = None
    if solvent_name and solvent_count > 0:
        solvent_mol = Molecule.from_file(f"{solvent_name}.pdb")
        solvent_mol.name = solvent_name.upper()
        for atom in solvent_mol.atoms:
            atom.metadata["residue_name"] = solvent_name.upper()

    ion_counts = ion_counts or {}
    ion_molecules = {
        ion_name: _make_ion_molecule(ion_name)
        for ion_name, count in ion_counts.items()
        if count > 0
    }

    if preserve_coordinate_order:
        molecule_by_residue_name = {lipid.name[:5]: lipid for lipid in lipid_molecules}
        if solvent_mol is not None:
            molecule_by_residue_name[solvent_name.upper()[:5]] = solvent_mol
        molecule_by_residue_name.update(ion_molecules)
        all_molecules = []
        for residue_name in _gro_residue_names(coordinate_file):
            molecule = molecule_by_residue_name.get(residue_name)
            if molecule is None:
                raise ValueError(
                    f"Could not map residue '{residue_name}' in {coordinate_file} "
                    "to a lipid, solvent, or ion molecule."
                )
            all_molecules.append(molecule)
    else:
        all_molecules = list(all_lipids)
        if solvent_mol is not None:
            all_molecules.extend([solvent_mol] * solvent_count)
        for ion_name, count in ion_counts.items():
            all_molecules.extend([ion_molecules[ion_name]] * count)

    topology = Topology.from_molecules(all_molecules)
    traj = mdtraj.load(coordinate_file)
    topology.set_positions(traj.xyz[0] * unit.nanometer)
    topology.box_vectors = np.array(box_dims) * 0.1 * unit.nanometer

    charge_from_molecules = list(lipid_molecules) + list(ion_molecules.values())
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message=r".*torch\.distributed\.reduce_op.*",
            category=FutureWarning,
        )
        warnings.filterwarnings(
            "ignore",
            message=r"Ambiguous failure while processing constraints.*",
            category=UserWarning,
        )
        interchange = Interchange.from_smirnoff(
            force_field=forcefield,
            topology=topology,
            charge_from_molecules=charge_from_molecules,
        )

        if hmr:
            interchange.to_gromacs(prefix=output_prefix, decimal=3, hydrogen_mass=3.0)
        else:
            interchange.to_gromacs(prefix=output_prefix)


def neutralize_with_gromacs(
    gro_file: str,
    top_file: str,
    solvent_group: str,
    gmx_executable: str = "gmx",
    maxwarn: int = 1,
) -> Tuple[str, Dict[str, int]]:
    """Use gmx genion to replace water with ions and return ion counts."""
    mdp_path = write_ions_mdp("ions.mdp")
    grompp_log = Path("grompp_ions.log")
    genion_log = Path("genion_ions.log")
    print(f"Wrote ion preparation MDP: {mdp_path}.")
    print(f"Running GROMACS grompp to create ions.tpr. Log: {grompp_log}")
    with open(grompp_log, "w") as log:
        try:
            subprocess.run(
                [
                    gmx_executable,
                    "grompp",
                    "-f",
                    "ions.mdp",
                    "-c",
                    gro_file,
                    "-p",
                    top_file,
                    "-o",
                    "ions.tpr",
                    "-maxwarn",
                    str(maxwarn),
                ],
                check=True,
                stdout=log,
                stderr=subprocess.STDOUT,
            )
        except subprocess.CalledProcessError as exc:
            raise RuntimeError(f"GROMACS grompp failed. See {grompp_log}.") from exc
    output_gro = "system_solv.gro"
    print(
        "Running GROMACS genion with "
        f"pname=NA, nname=CL, neutral=True; selecting solvent group '{solvent_group}'. "
        f"Log: {genion_log}"
    )
    with open(genion_log, "w") as log:
        try:
            subprocess.run(
                [
                    gmx_executable,
                    "genion",
                    "-s",
                    "ions.tpr",
                    "-o",
                    output_gro,
                    "-p",
                    top_file,
                    "-pname",
                    "NA",
                    "-nname",
                    "CL",
                    "-neutral",
                ],
                input=f"{solvent_group}\n",
                text=True,
                check=True,
                stdout=log,
                stderr=subprocess.STDOUT,
            )
        except subprocess.CalledProcessError as exc:
            raise RuntimeError(
                f"GROMACS genion failed while selecting solvent group '{solvent_group}'. "
                f"See {genion_log}."
            ) from exc
    residue_counts = _gro_residue_counts(output_gro)
    ion_counts = {
        "NA": residue_counts.get("NA", 0),
        "CL": residue_counts.get("CL", 0),
    }
    return output_gro, ion_counts


def run_packmol(
    packmol_input_file: str,
    parameterize: bool,
    force_field_file: str,
    hmr: bool,
    charge_model: str,
    gmx_executable: Optional[str] = None,
    solvent_group: Optional[str] = None,
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
    target_solvent_count = int(config['parameters'].get('target_solvent_count', solvent_count))
    packmol_solvent_count = int(config['parameters'].get('packmol_solvent_count', solvent_count))
    neutralize_ions = bool(config['parameters'].get('neutralize_ions', False))
    packmol_output = config['parameters']['output_name']
    total_lipid_count = sum(lipid_counts)
    estimated_bilayer_charge = float(config['parameters'].get('estimated_bilayer_charge', 0.0))
    estimated_ion_counts = config['parameters'].get(
        'estimated_neutralizing_ion_counts',
        {"NA": 0, "CL": 0},
    )
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

    if neutralize_ions:
        if not solvent_name or solvent_count <= 0:
            raise ValueError("Ion neutralization requires solvent so genion can replace water molecules.")
        print(
            "Ion neutralization summary before GROMACS: "
            f"net formal charge {_format_charge(estimated_bilayer_charge)} e; "
            f"estimated ions NA={estimated_ion_counts.get('NA', 0)}, "
            f"CL={estimated_ion_counts.get('CL', 0)}."
        )
        print("Writing pre-ionization GROMACS files...")
        _build_interchange(
            coordinate_file=packmol_output,
            output_prefix="pre_ions",
            box_dims=box_dims,
            lipid_names=lipid_names,
            lipid_counts=lipid_counts,
            solvent_name=solvent_name,
            solvent_count=solvent_count,
            forcefield=forcefield,
            lipid_molecules=lipid_molecules,
            charge_model=charge_model,
            hmr=False,
        )
        print("Neutralizing system with GROMACS genion...")
        ionized_gro, ion_counts = neutralize_with_gromacs(
            gro_file="pre_ions.gro",
            top_file="pre_ions.top",
            solvent_group=solvent_group or solvent_name.upper(),
            gmx_executable=gmx_executable or "gmx",
        )
        waters_replaced_by_ions = sum(ion_counts.values())
        solvent_count = solvent_count - waters_replaced_by_ions
        final_hydration_level = (
            solvent_count / total_lipid_count
            if total_lipid_count > 0
            else 0.0
        )
        config["parameters"]["solvent_count_after_genion"] = solvent_count
        config["parameters"]["waters_replaced_by_ions"] = waters_replaced_by_ions
        config["parameters"]["final_hydration_level"] = final_hydration_level
        config["parameters"]["final_hydration_delta"] = (
            solvent_count - target_solvent_count
        )
        config["parameters"]["ion_counts"] = ion_counts
        save_config(Path(config_path), config)
        print(f"GROMACS inserted ions: NA={ion_counts['NA']}, CL={ion_counts['CL']}.")
        print(
            "Hydration after ion replacement: "
            f"{final_hydration_level:.3f} waters/lipid "
            f"({solvent_count} waters; target {target_solvent_count}, "
            f"Packmol placed {packmol_solvent_count})."
        )
        _build_interchange(
            coordinate_file=ionized_gro,
            output_prefix="bilayer",
            box_dims=box_dims,
            lipid_names=lipid_names,
            lipid_counts=lipid_counts,
            solvent_name=solvent_name,
            solvent_count=solvent_count,
            forcefield=forcefield,
            lipid_molecules=lipid_molecules,
            charge_model=charge_model,
            hmr=hmr,
            ion_counts=ion_counts,
            preserve_coordinate_order=True,
        )
    else:
        _build_interchange(
            coordinate_file=packmol_output,
            output_prefix="bilayer",
            box_dims=box_dims,
            lipid_names=lipid_names,
            lipid_counts=lipid_counts,
            solvent_name=solvent_name,
            solvent_count=solvent_count,
            forcefield=forcefield,
            lipid_molecules=lipid_molecules,
            charge_model=charge_model,
            hmr=hmr,
        )

    if hmr:
        print("System saved with HMR as bilayer.top / bilayer.gro.")
    else:
        print("System saved without HMR as bilayer.top / bilayer.gro.")

    elapsed = time.time() - start_time
    print(f"Parameterization complete in {elapsed:.2f} seconds.")
