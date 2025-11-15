from lipidbuilder.builder import LipidSystemBuilder


builder = LipidSystemBuilder(
    result_directory="results/lipid_system",
    force_field_name="openff-2.2.0",
    water_model_name="tip3p",
    force_field_file="openff-2.2.0.offxml",
    charge_model="openff-gnn-am1bcc-0.1.0-rc.3.pt", # supported types are am1bcc, or the nagl models ending with .pt
    lipid_types=["MC3","DSPC", "CHOL", "PEG2000_C_DMG" ],
    lipid_composition= [31, 33, 6, 6, 25, 25, 2, 0],
    hydration_level=40,
    use_hmr=True,
)
    
# If you have a method to build/prep the system files
builder.setup() 
#builder.minimize() TODO minimize in openMM because faster and output to gmx?
