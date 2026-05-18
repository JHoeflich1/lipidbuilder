from pathlib import Path
from typing import List, Optional
import os
import time

from .system_setup import build_packmol_input, run_packmol
from .utils import DEFAULT_RESULTS_DIR

class LipidSystemBuilder:
    """
    A class representing a lipid system with a force field, water model,
    and thermodynamic state (pressure, temperature, ionic composition).
    """

    def __init__(
        self,
        result_directory: str | Path = DEFAULT_RESULTS_DIR,
        force_field_name: str = "openff-2.2.0",
        water_model_name: str = "tip3p",
        force_field_file: str | Path = "openff-2.2.0.offxml",
        charge_model: str = "am1bcc",
        simulation_platform: str = "openmm",
        water_model_file: Optional[str] = None,
        gmx_executable: Optional[str] = None,
        lipid_types: Optional[List[str]] = None,
        lipid_composition: Optional[List[int]] = None,
        hydration_level: int = 0, 
        neutralize_ions: bool = False,
        ion_solvent_group: Optional[str] = None,
        use_hmr: bool = False,
        hmr: Optional[bool] = None,
    ):
        """
        Initialize a lipid bilayer system builder.

        Parameters
        ----------
        result_directory : str
            The path to the top-level directory where results will be stored.
        force_field_name : str
            Identifier for the force field used.
        water_model_name : str
            Identifier for the water model.
        force_field_file : str or Path
            Force-field file name in lipidbuilder/data/forcefields or an absolute path.
        charge_model : str
            Charge model used (e.g., "AM1", "Gasteiger", etc.).
        simulation_platform : str, optional
            Simulation engine (default: "openmm").
        water_model_file : str, optional
            Path to the water model file, if separate.
        gmx_executable : str, optional
            Path to GROMACS executable (if using GROMACS).
        lipid_types : list of str, optional
            Lipid species included in the system.
        lipid_composition : list of int, optional
            Count of each lipid species in each leaflet [top, bottom], ordered as in lipid_types.
        hydration_level : int
            Target final number of water molecules per lipid.
        neutralize_ions : bool
            Replace waters with NA/CL ions using GROMACS genion -neutral after Packmol.
        ion_solvent_group : str
            GROMACS group selected by genion for solvent replacement. Defaults to the solvent residue name.
        use_hmr: default False
            Parametierize your system w/ HMR
        """

        # Normalize lipid_types input
        if isinstance(lipid_types, str):
            self.lipid_types = lipid_types.strip().split()
        else:
            self.lipid_types = lipid_types or []

        # Parse lipid composition
        if isinstance(lipid_composition, str):
            # Convert string like "64 64 15 15" to list of ints
            self.lipid_composition = [int(x) for x in lipid_composition.strip().split()]
        else:
            self.lipid_composition = lipid_composition or []

        # Check that each lipid has a top + bottom leaflet 
        expected_len = 2 * len(self.lipid_types)
        if len(self.lipid_composition) != expected_len:
            raise ValueError(
                f"Lipid composition must have exactly 2 numbers per lipid (top/bottom leaflets). "
                f"Expected {expected_len} numbers, got {len(self.lipid_composition)}."
            )


        self.force_field = force_field_name
        self.water_model = water_model_name
        self.force_field_file = str(force_field_file)
        self.water_model_file = water_model_file
        self.charge_model = charge_model
        self.simulation_platform = simulation_platform
        self.gmx_executable = gmx_executable
        self.hydration_level = hydration_level
        self.neutralize_ions = bool(neutralize_ions)
        self.ion_solvent_group = ion_solvent_group
        self.use_hmr = use_hmr if hmr is None else hmr

        # Define system naming and paths
        clean_charge_model = charge_model.removesuffix(".pt")
        lipid_str = "_".join(self.lipid_types)
        hmr_suffix = "-HMR" if self.use_hmr else ""
        self.system_name = f"{lipid_str}-{force_field_name}-{clean_charge_model}{hmr_suffix}"
        self.base_path = Path(result_directory) / self.system_name
        
        self.setup_dir = self.base_path / "setup"
        self.setup_prefix = self.setup_dir / self.system_name

    def setup(self):
        """
        Create coordinates for lipid system using Packmol as build engine.
        This method should create necessary directories and call packing + parameterization.
        """
        # Create the setup directory if it doesn't already exist
        self.setup_dir.mkdir(parents=True, exist_ok=True)
        original_dir = Path.cwd()
        os.chdir(self.setup_dir)

        try:
            packmol_input_file, box_dims = build_packmol_input(
                lipid_names=self.lipid_types,
                lipid_counts=self.lipid_composition,
                solvent_name=self.water_model,
                solvent_count=self.hydration_level * sum(self.lipid_composition),
                force_field_file=self.force_field_file,
                charge_model=self.charge_model,
                hmr=self.use_hmr,
                neutralize_ions=self.neutralize_ions,
                tolerance=2.0,
            )

            self.packmol_input_file = packmol_input_file
            self.box_dims = box_dims

            start_time = time.time()
            run_packmol(
                packmol_input_file=packmol_input_file,
                parameterize=True,
                force_field_file=self.force_field_file,
                hmr=self.use_hmr,
                charge_model=self.charge_model,
                gmx_executable=self.gmx_executable,
                solvent_group=self.ion_solvent_group,
                config_path="config.json",
            )

            elapsed_time = time.time() - start_time
            print(f"Coordinates created and parameterized. Completed in {elapsed_time:.2f} seconds.")
        finally:
            os.chdir(original_dir)

    def __repr__(self):
        return (
            f"<LipidSystemBuilder system='{self.system_name}', "
            f"{sum(self.lipid_composition)} lipids, "
            f"{self.hydration_level} waters/lipid>"
        )
