import pytest

from lipidbuilder.builder import LipidSystemBuilder

def test_lipid_system_builder_initialization(tmp_path):
    builder = LipidSystemBuilder(
        result_directory=str(tmp_path),
        force_field_name="openff-2.2.0",
        water_model_name="tip3p",
        force_field_file="openff-2.2.0.offxml",
        charge_model="AM1BCC",
        lipid_types=["POPC", "POPE"],
        lipid_composition=[64, 64, 16, 16],
        hydration_level=30,
        use_hmr=False,
    )

    assert builder.force_field == "openff-2.2.0"
    assert builder.water_model == "tip3p"
    assert builder.lipid_types == ["POPC", "POPE"]
    assert builder.lipid_composition == [64, 64, 16, 16]
    assert builder.hydration_level == 30
    assert builder.base_path == tmp_path / "POPC_POPE-openff-2.2.0-AM1BCC"
    assert builder.setup_dir.name == "setup"
    assert str(builder).startswith("<LipidSystemBuilder system='POPC_POPE")

def test_lipid_composition_mismatch(tmp_path):
    with pytest.raises(ValueError, match="2 numbers per lipid"):
        LipidSystemBuilder(
            result_directory=str(tmp_path),
            force_field_name="openff-2.2.0",
            water_model_name="tip3p",
            force_field_file="openff-2.2.0.offxml",
            charge_model="AM1BCC",
            lipid_types=["POPC", "POPE"],
            lipid_composition=[64],  # mismatch on purpose
            hydration_level=30,
            use_hmr=False,
        )


def test_hmr_suffix_is_added_to_system_name(tmp_path):
    builder = LipidSystemBuilder(
        result_directory=tmp_path,
        force_field_name="openff-2.2.0",
        water_model_name="tip3p",
        force_field_file="openff-2.2.0.offxml",
        charge_model="am1bcc",
        lipid_types=["POPC"],
        lipid_composition=[1, 1],
        use_hmr=True,
    )

    assert builder.system_name.endswith("-HMR")
