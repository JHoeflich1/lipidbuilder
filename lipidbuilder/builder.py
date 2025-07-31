from pathlib import Path
from openff.units import unit

class LipidSystemBuilder:
  """
  A class representing a lipid system with a force field, water model, 
  thermodynamic state (pressure, temp, ionic composition).
  """
    def __init__(
      self, 
      result_directory: str,
      force_field_name: str,
      water_model_name: str,
      force_field_file: str, 
      water_model_file: str = None,
      charge_model: str,
      simulation_platform: str = "openmm",
      gmx_executable: str= None,
      lipid_types: list[str] = None,
      lipid_composition: list[int] = None,
      hydration_level: int = 0,
      ion: str = None,
      ion_concentration: list[int] = None
    ):
      """
      Initialize a lipid bilayer system builder.
  
      Parameters
      ----------
      result_directory : str
          The path to the top level directory where results will be stored.
      force_field_name : str
          Identifier for the force field used.
      water_model_name : str
          Identifier for the water model.
      force_field_file : str
          Path to the force field XML or GMX file.
      water_model_file : str, optional
          Path to the water model file, if separate.
      simulation_platform : str, optional
          Simulation engine (default: "openmm").
      gmx_executable : str, optional
          Path to GROMACS executable (if using GROMACS).
      lipid_types : list of str
          Lipid species included in the system.
      lipid_composition : list of int
          Count of each lipid species, ordered as in lipid_types.
      hydration_level : int
          Number of water molecules per lipid.
      ion : str, optional
          Ion type to add (e.g., "Na+", "Cl-").
      """

      self.force_field = force_field_name
      self.water_model = water_model_name
      self.force_field_file = force_field_file
      self.water_model_file = water_model_file
      self.charge_model = charge_model
      self.lipid_types = lipid_types
      self.lipid_composition = lipid_composition
      self.hydration_level = hydration_level
      self.ion = ion
      self.ionic_concentration = ion_concentration


      # Create a directory to store results for system
      #strip lipid types of spaces
      lipid_list = lipid_types.strip().split()
      lipid_str = "_".join(lipid_list) 
      
      self.system_name = f"{lipid_str}-{force_field_name}-{charge_model}"
      self.base_path = Path(result_directory, self.system_name)

      self.setup_dir = Path(self.base_path, "setup")
      self.setup_prefix = Path(self.setup_dir, self.system_name)


    def setup(self):
      """
      Create coordinates for lipid system using Packmol as build engine 
      """
      # Create the setup directory if it doesn't already exist
      self.setup_dir.mkdir(parents=True, exist_ok=True)

      # And then initialize packmol 
      # And then parameterize your system 
      # And then save to output 
