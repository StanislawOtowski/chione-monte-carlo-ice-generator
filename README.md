# Monte Carlo of proton configurations in ice Ic and VII

[Wersja polska](README.pl.md)

A Monte Carlo simulation (Metropolis algorithm) of hydrogen arrangements in the ice Ic and ice VII lattices, with gradual cooling. The energy is a penalty for breaking the ice rules: each oxygen should have exactly 2 hydrogens next to it. The output is a set of CIF files with the final configuration.

## Installation

```
pip install -r requirements.txt
```

## Running

```
python main.py
```

Parameters are set in `main.py`:

- in the `if __name__=='__main__':` block, the number of cores (`Pool(4)`) and the list of lattice sizes (`[4]`), which are computed in parallel, one task per size;
- in the arguments of `generate_and_iterate_net`:
  - `ice_type`: `'1c'` or `'7'`,
  - `temperatures`: a list of kT values (from highest to lowest),
  - `iterations_at_one_temperature_at_stable_state`: the number of samples per temperature
    (before sampling, half that many sweeps are run to reach equilibrium),
  - `cif_file_generation`: whether to save a CIF file.

A lattice size of `net_size = n` means n×n×n unit cells: 8n³ oxygens for ice Ic and 16n³ for ice VII (two interpenetrating Ic sublattices).

## Output

- Printed to the console for each temperature: the energy, the per-molecule energy moments (⟨e⟩…⟨e⁴⟩) and the
  per-molecule specific heat cV.
- CIF files go to `Plotting/CIF/`, named `Ice_<type>_<number of oxygens>_oxygens_<date>_<time>.<random>.cif`.
- At the end, a 3D plot of the structure opens (`plt.show()` blocks until the window is closed).

## Structure

| file | role |
|---|---|
| `main.py` | entry point: lattice generation, MC, CIF output, plot |
| `Generating/Ic_ice_net_generator.py` | generates an Ic lattice with a random hydrogen configuration |
| `Ice_Structure_Changes/change_hydrogen_configuration.py` | MC loop (Metropolis) over temperatures |
| `Analyzys/energy_calculations.py` | total energy and local energy change |
| `Plotting/saving_to_CIF_file.py` | writes CIF files from the templates `Ice_1c_template.txt`, `Ice_7_template.txt` |
| `Plotting/Ic_and_7_ice_net_plotter.py` | 3D plot |
