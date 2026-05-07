import json
import os
from pathlib import Path

from lipidbuilder.system_setup import build_packmol_input


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
