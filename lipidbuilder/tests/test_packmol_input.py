import json
import os
from pathlib import Path

from lipidbuilder.system_setup import (
    build_packmol_input,
    estimate_bilayer_charge,
    estimate_neutralizing_ion_counts,
    write_ions_mdp,
)


def test_build_packmol_input(tmp_path):
    original_dir = Path.cwd()
    test_lipids = ["POPC"]
    test_counts = [2, 2]
    test_solvent = "tip3p"
    test_solvent_count = 10

    os.chdir(tmp_path)
    try:
        inp_path, dims = build_packmol_input(
            lipid_names=test_lipids,
            lipid_counts=test_counts,
            solvent_name=test_solvent,
            solvent_count=test_solvent_count,
            force_field_file="openff-2.2.0.offxml",
            charge_model="am1bcc",
            hmr=False,
        )
    finally:
        os.chdir(original_dir)

    assert Path(inp_path).is_file()
    assert dims[0] > 0
    assert dims[1] > 0
    assert dims[2] > 0
    assert (tmp_path / "config.json").is_file()

    config = json.loads((tmp_path / "config.json").read_text())

    assert config["parameters"]["lipid_names"] == test_lipids
    assert config["parameters"]["solvent_name"] == test_solvent
    assert (tmp_path / "POPC.pdb").is_file()
    assert (tmp_path / "tip3p.pdb").is_file()


def test_build_packmol_input_with_ions(tmp_path):
    original_dir = Path.cwd()

    os.chdir(tmp_path)
    try:
        inp_path, _ = build_packmol_input(
            lipid_names=["POPC"],
            lipid_counts=[2, 2],
            solvent_name="tip3p",
            solvent_count=10,
            force_field_file="openff-2.2.0.offxml",
            charge_model="am1bcc",
            hmr=False,
            neutralize_ions=True,
        )
    finally:
        os.chdir(original_dir)

    packmol_input = Path(inp_path).read_text()
    config = json.loads((tmp_path / "config.json").read_text())

    assert config["parameters"]["neutralize_ions"] is True
    assert config["parameters"]["target_solvent_count"] == 10
    assert config["parameters"]["estimated_bilayer_charge"] == 0.0
    assert config["parameters"]["packmol_solvent_count"] == 10
    assert config["parameters"]["hydration_level_basis"] == "final_target"
    assert config["parameters"]["positive_ion_name"] == "NA"
    assert config["parameters"]["negative_ion_name"] == "CL"
    assert "structure na.pdb" not in packmol_input
    assert "structure cl.pdb" not in packmol_input


def test_build_packmol_input_adds_extra_waters_for_neutralization(tmp_path, monkeypatch):
    monkeypatch.setattr(
        "lipidbuilder.system_setup.estimate_neutralizing_ion_counts",
        lambda lipid_names, lipid_counts: {"NA": 2, "CL": 0},
    )
    original_dir = Path.cwd()

    os.chdir(tmp_path)
    try:
        inp_path, _ = build_packmol_input(
            lipid_names=["POPC"],
            lipid_counts=[1, 1],
            solvent_name="tip3p",
            solvent_count=10,
            force_field_file="openff-2.2.0.offxml",
            charge_model="am1bcc",
            hmr=False,
            neutralize_ions=True,
        )
    finally:
        os.chdir(original_dir)

    packmol_input = Path(inp_path).read_text()
    config = json.loads((tmp_path / "config.json").read_text())

    assert config["parameters"]["target_solvent_count"] == 10
    assert config["parameters"]["estimated_bilayer_charge"] == 0.0
    assert config["parameters"]["estimated_neutralizing_ion_counts"] == {"NA": 2, "CL": 0}
    assert config["parameters"]["estimated_neutralizing_ions"] == 2
    assert config["parameters"]["packmol_solvent_count"] == 12
    assert "  number 6" in packmol_input


def test_estimate_bilayer_charge_and_neutralizing_ion_counts(tmp_path):
    lipid_library = tmp_path / "lipids.csv"
    lipid_library.write_text(
        "Name,Smiles String\n"
        "ANION,CC(=O)[O-]\n"
        "NEUTRAL,O\n"
    )

    assert estimate_bilayer_charge(
        ["ANION"],
        [1, 1],
        lipid_library_path=lipid_library,
    ) == -2.0
    assert estimate_neutralizing_ion_counts(
        ["ANION"],
        [1, 1],
        lipid_library_path=lipid_library,
    ) == {"NA": 2, "CL": 0}
    assert estimate_neutralizing_ion_counts(
        ["NEUTRAL"],
        [1, 1],
        lipid_library_path=lipid_library,
    ) == {"NA": 0, "CL": 0}


def test_write_ions_mdp(tmp_path):
    mdp_path = write_ions_mdp(tmp_path / "ions.mdp")
    mdp = mdp_path.read_text()

    assert "integrator  = steep" in mdp
    assert "coulombtype     = cutoff" in mdp
    assert "pbc             = xyz" in mdp
