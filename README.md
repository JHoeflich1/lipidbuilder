# lipidbuilder

`lipidbuilder` builds hydrated lipid bilayer starting structures with Packmol and
parameterizes them with OpenFF.

## Layout

- `lipidbuilder/` - Python package source.
- `lipidbuilder/data/available-lipids/` - lipid metadata, lipid PDBs, and solvent PDBs.
- `lipidbuilder/data/forcefields/` - packaged force-field files.
- `lipidbuilder/results/lipid_system/` - default output folder for generated systems.
- `lipidbuilder/How_to_use.py` - runnable example.
- `pixi.toml` - reproducible development and runtime environment.

## Create the Environment

```bash
pixi install
```

Run commands inside the Pixi environment:

```bash
pixi run list-lipids
pixi run test
pixi run build-example
```

## Python Example

```python
from lipidbuilder import LipidSystemBuilder

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
```

`force_field_file` can be either a filename in `lipidbuilder/data/forcefields/`
or an absolute path to another `.offxml` file.

For each lipid, `lipid_composition` takes two values: top leaflet count, then
bottom leaflet count. For example, `lipid_types=["POPC", "POPE"]` expects
`lipid_composition=[POPC_top, POPC_bottom, POPE_top, POPE_bottom]`.

Generated systems are written to:

```text
<result_directory>/<lipids>-<force_field>-<charge_model>[-HMR]/setup/
```
