import os
import pandas as pd
import sys
import numpy as np
import glob
import time
import shutil
from pathlib import Path
import subprocess
from openff.toolkit import ForceField, Molecule, Topology
from openff.interchange import Interchange
from openff.toolkit.utils.nagl_wrapper import NAGLToolkitWrapper
import MDAnalysis as mda
from MDAnalysis.analysis import distances



import subprocess
import os
import glob

from lipidbuilder.utils import SCRATCH_DIR, LIPID_PDB_DIR

def pullNewLipid(Lipid_name, pull_atom):
    """
    Runs the gmx pulling command to straighten out the lipid and cleans up generated files.
    All outputs are dumped into a .log file in ./scratch.
    
    Inputs:
        Lipid_name: shorthand lipid name used for residue names, e.g., POPC
        pull_atom: atom in the lipid headgroup used to straighten the lipid
    
    Returns:
        Path to the final PDB file
    """

    os.makedirs(SCRATCH_DIR, exist_ok=True)

    log_file = f"{SCRATCH_DIR}/{Lipid_name}_pull.log"

    with open(log_file, "w") as log:
        # 1. make_ndx
        command_ndx = f"gmx make_ndx -f {SCRATCH_DIR}/{Lipid_name}.gro -o {SCRATCH_DIR}/{pull_atom}.ndx"
        subprocess.run(
            f"echo 'a {pull_atom}\nname 3 pull_atom\nq\n' | {command_ndx}",
            shell=True, stdout=log, stderr=subprocess.STDOUT
        )

        # 2. grompp
        command_grompp = (
            f"gmx grompp -f scripts/runpull.mdp -c {SCRATCH_DIR}/{Lipid_name}.gro -p {SCRATCH_DIR}/{Lipid_name}.top -n {SCRATCH_DIR}/{pull_atom}.ndx -o {SCRATCH_DIR}/pull.tpr -maxwarn 2")
        subprocess.run(command_grompp, shell=True, stdout=log, stderr=subprocess.STDOUT)
        
        # 3. mdrun
        command_pull = f"gmx mdrun -deffnm {SCRATCH_DIR}/pull"
        subprocess.run(command_pull, shell=True, stdout=log, stderr=subprocess.STDOUT)

        # # move mdout when generated
        # if os.path.exists("pull.mdout"):
        #     shutil.move("pull.mdout", f"{SCRATCH_DIR}/pull.mdout")
        if os.path.exists("mdout.mdp"):
            shutil.move("mdout.mdp", f"{SCRATCH_DIR}/mdout.mdp")

        # 4. editconf to PDB
        command_pdb = f"gmx editconf -f {SCRATCH_DIR}/pull.gro -o {SCRATCH_DIR}/{Lipid_name}.pdb"
        subprocess.run(command_pdb, shell=True, stdout=log, stderr=subprocess.STDOUT)

    # Clean up temporary files in notebook directory
    for file in SCRATCH_DIR.glob('pull*'):
        if os.path.exists(file):
            os.remove(file)

    for file in SCRATCH_DIR.glob('*.ndx'):
        os.remove(file)

    return



class newLipid(object):
    """ Saving lipid parameters into an object"""

    def __init__(self, name, headgroup_atom, headgroup_atom_index, tailgroup_atom,tailgroup_atom_index, distance, experimental_density, smiles):

        self.name = name
        self.headgroup_atom = headgroup_atom
        self.headgroup_atom_index = headgroup_atom_index
        self.tailgroup_atom = tailgroup_atom
        self.tailgroup_atom_index = tailgroup_atom_index
        self.distance = distance
        self.experimental_density = experimental_density #this is in g/cm^3
        self.smiles = smiles


def calcLipidLength(lipid, Lipid_name):
    '''Calculate the distance between headgroup and tailgroup atoms of a lipid,
       and update the lipid object with this information
       
       Inputs:
       lipid: An object representing the lipid with attributes headgroup_atom and tailgroup_atom
       Lipid_name: The name of the lipid used to locate the corresponding pdb
    '''
    pdb_path = SCRATCH_DIR / f"{Lipid_name}.pdb"

    import warnings
    
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=UserWarning)
        u_pdb = mda.Universe(pdb_path)

    hg_atom = lipid.headgroup_atom
    tg_atom = lipid.tailgroup_atom
    head_group = u_pdb.select_atoms(f'name {hg_atom}')
    tail_group = u_pdb.select_atoms(f'name {tg_atom}')
    
    # calculate the distance and save into object
    calc_distance = distances.distance_array(head_group.positions, tail_group.positions)
    lipid.distance = calc_distance[0][0] 
    
    # write in the atom index into object
    lipid.headgroup_atom_index = head_group[0].index 
    lipid.tailgroup_atom_index = tail_group[0].index 

def lipidToDict(lipid):
    return {
        'Name': lipid.name,
        'Headgroup Atom': lipid.headgroup_atom,
        'Headgroup Atom Index': lipid.headgroup_atom_index,
        'Tailgroup Atom': lipid.tailgroup_atom,
        'Tailgroup Atom Index': lipid.tailgroup_atom_index,
        'HG/TG distance': lipid.distance,
        'Experimental Density': lipid.experimental_density,
        'Smiles String': lipid.smiles
    }

def loadExistingData():
    csv_file_path = Path("PulledLipid.csv")

    if csv_file_path.exists():
        df = pd.read_csv(csv_file_path)
    else:
        df = pd.DataFrame(columns=[
            'Name',
            'Headgroup Atom', 'Headgroup Atom Index',
            'Tailgroup Atom', 'Tailgroup Atom Index',
            'HG/TG distance',
            'Experimental Density',
            'Smiles String'
        ])

    return df, csv_file_path

def saveLipidCsv(lipid):
    """Save lipid info to CSV file and place its PDB into lipid_pdbs/lipid/lipid.pdb."""
    #create and add to csv
    df, csv_file_path = loadExistingData()
    lipid_dict = lipidToDict(lipid)

    if lipid_dict['Name'] not in df['Name'].values:
        df = pd.concat([df, pd.DataFrame([lipid_dict])], ignore_index=True)
        df.to_csv(csv_file_path, index=False)
        print(f"Saved lipid '{lipid.name}' to CSV at {csv_file_path}")
    else:
        print(f"Lipid '{lipid.name}' already exists in {csv_file_path}")

    # Expected PDB location (same as your pulling function)
    pdb_source = SCRATCH_DIR / f"{lipid.name}.pdb"

    if not pdb_source.exists():
        print(f"ERROR: Could not find PDB for '{lipid.name}' at {pdb_source}")
        return

    # Folder structure: lipid_pdbs/Lipid/
    lipid_pdb_dir = LIPID_PDB_DIR / lipid.name
    lipid_pdb_dir.mkdir(parents=True, exist_ok=True)

    # Destination: lipid_pdbs/Lipid/Lipid.pdb
    pdb_dest = lipid_pdb_dir / f"{lipid.name}.pdb"

    shutil.copy(pdb_source, pdb_dest)
    print(f"PDB saved to: {pdb_dest}")