# WAVIConstructor.jl
[![Documentation](https://img.shields.io/badge/Docs%20-github.io%2FWAVIConstructor.jl%2F-red)](https://WAVI-ice-sheet-model.github.io/WAVIConstructor.jl/)
![Work in Progress](https://img.shields.io/badge/status-work--in--progress-orange)
![Experimental](https://img.shields.io/badge/status-experimental-blueviolet)
![License: MIT](https://img.shields.io/badge/license-MIT-green)
![Testing](https://github.com/WAVI-ice-sheet-model/WAVIConstructor.jl/actions/workflows/test.yml/badge.svg)
![Documentation](https://github.com/WAVI-ice-sheet-model/WAVIConstructor.jl/actions/workflows/docs.yml/badge.svg)


Pre-processing library for [WAVI.jl](https://github.com/WAVI-ice-sheet-model/WAVI.jl).

This package is a WIP, but will support the reading and pre-processing of input data, and preparing a WAVI domain for inversion runs.

## Installation

You can install the latest version of `WAVIConstructor` using Julia's in-built package manager:

```julia
using Pkg
Pkg.add(PackageSpec(url="https://github.com/WAVI-ice-sheet-model/WAVIConstructor.jl.git", rev = "main"))
```

## Usage

The main workflow involves three steps:

1. **Initialising grids and loading data** using `init_bedmachine`
2. **Domain selection** using `select_domain_wavi`
3. **Preparing WAVI inversion input files** using `setup_wavi_data`

The `init_bedmachine` and `select_domain_wavi` functions are wrapped within `setup_wavi_data`.

All configurations (data sources, domains, grid spacing, file output) are defined in the constructor params.

### Simple Workflow

```julia
using WAVIConstructor

# Define all data contruction parameters
params = default_constructor_params(
    bedmachine_file = "Data/BedMachineAntarctica-v3.nc",
    dx = 10000.0,  # Grid spacing in metres
    basins = [1, 2, 3],  # Zwally drainage basins to include
    dhdt_data = "Smith",  # Elevation change dataset
    veloc_data = "Measures_2016/17",  # Velocity dataset
    temps = "BISICLES_8km",  # Temperature dataset
    output_path = "wavi_inversion_input"
)

# Complete workflow: initialise, select domain, and write WAVI input files
Gh, Gu, Gv, Gc = setup_wavi_data(params)
```

### Running step-by-step

You can also run the workflow step-by-step for more control and intermediate inspection:

```julia
# Step 1: Initialise grids
Gh, Gu, Gv, Gc = init_bedmachine(params)

# Step 2: Select domain
Gh, Gu, Gv, Gc = select_domain_wavi(Gh, Gu, Gv, Gc, params)

# Step 3: Process and write files
Gh, Gu, Gv, Gc = setup_wavi_data(params)
```

### Loading data only
You can also use the low-level data loading functions directly to inspect specific datasets:

```julia
using WAVIConstructor.DataLoading

# Load specific datasets
bed, x, y, geoid, mask, thickness, surface, firn = get_bedmachine("Data/BedMachineAntarctica-v3.nc")
xx, yy, vx, vy = get_measures_velocities("Data/velocity_file.nc")
```
## Outputs

After running `setup_wavi_data()`, you'll have binary `.bin` files in your output directory containing:

- Grid geometry (coordinates, spacing)
- Ice thickness
- Bed topography
- Surface elevation
- Velocities (u, v components)
- Temperature fields
- Accumulation rates
- Domain masks

These files are in the format expected by WAVI.jl.

## Configuration Parameters

Key parameters you can set:

| Parameter | Description | Default |
|-----------|-------------|---------|
| `bedmachine_file` | Path to BedMachine NetCDF file | Required |
| `dx` | Grid spacing (m) | Required |
| `basins` | Zwally basin IDs to include | `1:27` (all) |
| `dhdt_data` | Elevation change dataset | `"Smith"` |
| `smith_dhdt_dir` | Directory containing Smith dh/dt files | `"Data/Smith_2020_dhdt"` |
| `veloc_data` | Velocity dataset | `"Measures_2016/17"` |
| `temps` | Temperature dataset | `"BISICLES_8km"` |
| `output_path` | Directory for output files | `"outputs"` |
| `clip_edge_padding` | Edge padding when clipping | `3` |
| `density_ice` | Ice density (kg/m³) | `918.0` |
| `density_ocean` | Ocean density (kg/m³) | `1028.0` |
| `min_thick` | Minimum ice thickness (m) | `50.0` |

### Velocity Dataset Options

- `"Measures"`: MEaSUREs from .mat file

### Temperature Dataset Options

- `"Frank"`: Frank's temperature data (.mat file)
- `"BISICLES_8km"`: BISICLES 8km resolution
- `"BISICLES_1km"`: BISICLES 1km resolution

## Obtaining Data
TBC

## Contributing
TBC

## License
This software is licensed under an MIT license.

For more information please see the attached `LICENSE` file.

## Developers
WAVI.jl Ice Sheet Model Team, British Antarctic Survey
