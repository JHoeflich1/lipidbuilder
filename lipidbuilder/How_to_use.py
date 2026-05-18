# Run "python How_to_use.py" to build a lipid bilayer system with the specified parameters.


from lipidbuilder.builder import LipidSystemBuilder


builder = LipidSystemBuilder(
    result_directory="results/lipid_system",  # Directory where all generated results (structures, topologies, trajectories) will be saved
    force_field_name="openff-2.2.0", # Name of the force field
    water_model_name="tip3p", # Water model to use for solvation
    force_field_file="openff-2.2.0.offxml", # Path to the specific force field file
    charge_model="openff-gnn-am1bcc-1.0.0.pt", # supported charge methods are openFF am1bcc, or the nagl models ending with .pt. Default charge model openff-gnn-am1bcc-1.0.0.pt
    lipid_types=["POPC"], # List of lipid types to include in the bilayer. ex [lipid1, lipid2]
    lipid_composition= [34, 34],  # Composition of each lipid type, split by leaflet (top and bottom), ex [lipid1_top, lipid1_bottom, lipid2_top, lipid2_bottom,.. ]
    hydration_level=40, # Number of water molecules per lipid (hydration level)
    neutralize_ions=True, # Whether to add ions to neutralize the system, will not neutrlize if the system is already neutral.
    use_hmr=False, # Whether to use Hydrogen Mass Repartitioning (HMR) to allow larger time step
)
    
# Setup the system: builds the lipid bilayer with specified composition, solvates it, 
# and writes the initial structure and topology files to result_directory
builder.setup() 

## TODO: minimize the system in OpenMM
# This step is important to remove bad contacts and relax the system before MD
# builder.minimize()