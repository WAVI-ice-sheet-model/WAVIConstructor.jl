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

# Define all data construction parameters
params = default_constructor_params(
    bedmachine_file = "/path/to/data/BedMachineAntarctica-v3.nc",
    albmap_file = "/path/to/data/ALBMAPv1.nc",
    zwally_file = "/path/to/data/DrainageBasins/ZwallyBasins.mat",
    smith_dhdt_dir = "/path/to/data/Smith_2020_dhdt",  # Optional: directory with Smith dhdt data
    measures_velocity_file = "/path/to/data/antarctica_ice_velocity_2016_2017_1km_v01.nc",  # Optional: MEaSUREs velocity
    bisicles_temps_file = "/path/to/data/antarctica-bisicles-xyzT-8km.nc",  # Optional: BISICLES temperature
    dx = 10000.0,  # Grid spacing in metres
    basins = [1, 2, 3],  # Zwally drainage basins to include
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

# Load specific datasets (all file paths should be full paths)
bed, x, y, geoid, mask, thickness, surface, firn = get_bedmachine("/path/to/data/BedMachineAntarctica-v3.nc")
xx, yy, vx, vy = get_measures_velocities("/path/to/data/velocity_file.nc")
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
| `bedmachine_file` | Full path to BedMachine NetCDF file | `"Data/BedMachineAntarctica-v3.nc"` |
| `albmap_file` | Full path to ALBMAP file | `"Data/ALBMAPv1.nc"` |
| `zwally_file` | Full path to Zwally basins file | `"Data/DrainageBasins/ZwallyBasins.mat"` |
| `arthern_file` | Full path to Arthern accumulation file | `"Data/amsr_accumulation_map.txt"` |
| `smith_dhdt_dir` | Full path to directory containing Smith dh/dt files (optional) | `"Data/Smith_2020_dhdt"` |
| `frank_temps_file` | Full path to Frank temperatures file (optional, alternative to bisicles_temps_file) | `"Data/FranksTemps.mat"` |
| `measures_velocity_file` | Full path to MEaSUREs velocity NetCDF file (optional) | `"Data/antarctica_ice_velocity_2016_2017_1km_v01.nc"` |
| `bisicles_temps_file` | Full path to BISICLES temperature file (optional, alternative to frank_temps_file) | `"Data/antarctica-bisicles-xyzT-8km.nc"` |
| `dx` | Grid spacing (m) | `10000.0` |
| `basins` | Zwally basin IDs to include | `1:27` (all) |
| `output_path` | Directory for output files | `"wavi_input"` |
| `clip_edge_padding` | Edge padding when clipping | `3` |
| `density_ice` | Ice density (kg/m³) | `918.0` |
| `density_ocean` | Ocean density (kg/m³) | `1028.0` |
| `min_thick` | Minimum ice thickness (m) | `50.0` |

### Dataset options:
- Temperature data: Choose between `frank_temps_file` (Frank's temperature data) or `bisicles_temps_file` (BISICLES temperature data). Temperature data is required - you must provide one of these.
- Velocity data: Provide `measures_velocity_file` (MEaSUREs NetCDF format). If not provided, zero velocities will be used.
- Elevation change data: Provide `smith_dhdt_dir` to use Smith et al. 2020 elevation change data. If not provided, zeros will be used.

## Obtaining Data
TBC

## Contributing
TBC

## License
This software is licensed under an MIT license.

For more information please see the attached `LICENSE` file.

## Developers
WAVI.jl Ice Sheet Model Team, British Antarctic Survey
