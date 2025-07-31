from pathlib import Path

import numpy
import openmm
from openff.toolkit import ForceField, Molecule, Topology
from openff.units import unit


# functions for seting up packmol input file, running and minimizing 
