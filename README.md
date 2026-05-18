# lipidbuilder

`lipidbuilder` builds hydrated lipid bilayer starting structures with Packmol and
parameterizes them with OpenFF.

## Layout

- `lipidbuilder/` - Python package source.
- `lipidbuilder/data/available-molecules/` - lipid metadata, lipid PDBs, and solvent PDBs.
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

cd lipidbuilder
pixi run python How_to_use.py

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
    neutralize_ions=True,
    use_hmr=True,
)

builder.setup()
```

`force_field_file` can be either a filename in `lipidbuilder/data/forcefields/`
or an absolute path to another `.offxml` file.

Set `neutralize_ions=True` to run GROMACS after Packmol. The builder writes
`ions.mdp`, creates a temporary `pre_ions.gro`/`pre_ions.top`, runs `gmx grompp`
and `gmx genion -neutral`, then rebuilds the final OpenFF Interchange system
from `system_solv.gro` with the inserted `NA`/`CL` ions. `hydration_level` is
treated as the target final water/lipid ratio: when neutralization is enabled,
Packmol starts with extra waters equal to the estimated neutralizing ion count,
then the final config reports `ion_counts`, `waters_replaced_by_ions`, and
`final_hydration_level`. By default `genion` selects the solvent group using the
solvent residue name, such as `TIP3P`; override this with `ion_solvent_group`
if your GROMACS groups differ. GROMACS command output is written to
`grompp_ions.log` and `genion_ions.log`.

For each lipid, `lipid_composition` takes two values: top leaflet count, then
bottom leaflet count. For example, `lipid_types=["POPC", "POPE"]` expects
`lipid_composition=[POPC_top, POPC_bottom, POPE_top, POPE_bottom]`.

Generated systems are written to:

```text
<result_directory>/<lipids>-<force_field>-<charge_model>[-HMR]/setup/
```
