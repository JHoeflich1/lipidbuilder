import os
import json
from pathlib import Path
from system_setup import build_packmol_input


def test_build_packmol_input(tmp_path):
    # Arrange: Set working dir and test data
    os.chdir(tmp_path)  # Run test in isolated tmp directory
    test_lipids = ["POPC", "POPE"]
    test_counts = [2, 2]
    test_solvent = "TIP3P"
    test_solvent_count = 10

    # Provide a small, dummy lipid CSV file
    lipid_csv = tmp_path / "PulledLipid.csv"
    lipid_csv.write_text(
        "Name,Lipid Smiles String,HG/TG distance,Headgroup Atom Index,Tailgroup Atom Index,Approximate Volume,Net Charge\n"
        "POPC,C(COC)=O,20,1,10,1000,0\n"
        "POPE,C(COC)=O,20,1,10,1000,0\n"
        "TIP3P,O,0,0,0,30,0\n"
    )

    os.makedirs(tmp_path / "data/available-lipids/lipids_parameterized/POPC")
    os.makedirs(tmp_path / "data/available-lipids/lipids_parameterized/POPE")
    os.makedirs(tmp_path / "data/available-lipids/lipids_parameterized/TIP3P")

    for name in test_lipids + [test_solvent]:
        pdb = tmp_path / f"data/available-lipids/lipids_parameterized/{name}/{name}.pdb"
        top = tmp_path / f"data/available-lipids/lipids_parameterized/{name}/{name}.top"
        pdb.write_text("REMARK dummy pdb\n")
        top.write_text("; dummy top\n")

    # Act: Build Packmol input
    inp_path, out_path, dims = build_packmol_input(
        lipid_names=test_lipids,
        lipid_counts=test_counts,
        solvent_name=test_solvent,
        solvent_count=test_solvent_count,
        lipid_library_csv=str(lipid_csv),
    )

    # Assert: Check output files and values
    assert Path(inp_path).is_file()
    assert dims[2] > 0
    assert Path("config.json").is_file()

    # Load config and check contents
    with open("config.json") as f:
        config = json.load(f)

    assert config["parameters"]["lipid_names"] == test_lipids
    assert config["parameters"]["solvent_name"] == test_solvent
