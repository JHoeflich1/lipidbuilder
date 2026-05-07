from lipidbuilder.builder import LipidSystemBuilder


builder = LipidSystemBuilder(
    result_directory="lipidbuilder/results/lipid_system",
    force_field_name="openff-2.2.0",
    water_model_name="tip3p",
    force_field_file="openff-2.2.0.offxml",
    charge_model="am1bccelf10",
    lipid_types=["POPC"],
    lipid_composition=[64, 64],
    hydration_level=40,
    use_hmr=True,
)

builder.setup()
