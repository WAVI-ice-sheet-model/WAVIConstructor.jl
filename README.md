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

All configurations (data sources, domains, grid spacing, file output) are defined via the `ConstructorParams` struct. Each geophysical dataset is specified by a **source type** singleton (e.g. `BedMachineV3()`), optionally paired with a custom file path as a tuple. Use `NoData()` for optional categories you want to skip.

### Simple Workflow

```julia
using WAVIConstructor

# Create parameters — each dataset is just a source singleton (default paths built-in)
params = default_constructor_params(
    bed          = BedMachineV3(),
    surface_temp = ALBMAPv1(),
    temperature  = BISICLESTemps(),
    velocity     = MEaSUREs(),
    accumulation = ArthernAccumulation(),
    dhdt         = SmithDhdt(),
    basins       = ZwallyBasins(),
    dx           = 10000.0,          # Grid spacing in metres
    basin_ids    = [1, 2, 3],        # Zwally drainage basins to include
    output_path  = "wavi_inversion_input",
)

# Complete workflow: initialise, select domain, and write WAVI input files
Gh, Gu, Gv, Gc = setup_wavi_data(params)
```

Since all source singletons carry a built-in default path, the above is equivalent to just overriding what you need:

```julia
params = default_constructor_params(
    temperature = BISICLESTemps(),   # uses default path
    basin_ids   = [4],
)
```

To supply a custom file path, pass a `(source, path)` tuple:

```julia
params = default_constructor_params(
    temperature = (BISICLESTemps(), "my/custom/temps.nc"),
    basin_ids   = [4],
)
```

### Skipping optional datasets

Use `NoData()` for categories you don't need. Velocity and dh/dt default to zeros when skipped:

```julia
params = default_constructor_params(
    velocity = NoData(),   # skip velocity — zeros used
    dhdt     = NoData(),   # skip dh/dt   — zeros used
)
```

### Running step-by-step

You can also run the workflow step-by-step for more control and intermediate inspection:

```julia
d = to_dict(params)                           # convert to Dict for the internal API
Gh, Gu, Gv, Gc = init_bedmachine(d)          # Step 1: load & interpolate data
Gh, Gu, Gv, Gc = select_domain_wavi(Gh, Gu, Gv, Gc, d)  # Step 2: mask domain
```

### Loading individual datasets

Use the dispatch-based `load_data` interface to inspect a single dataset:

```julia
using WAVIConstructor

# Each source singleton selects the right loader
bm = load_data(BedMachineV3(), "Data/BedMachineAntarctica-v3.nc")
bm.bed   # bed topography array
bm.x     # x coordinates

vel = load_data(MEaSUREs(), "Data/Antarctica_ice_velocity_2014_2015_1km_v01.nc")
vel.vx   # x-velocity component
```

## Outputs

`setup_wavi_data()` writes all prepared fields to your `output_path` directory.
The output format is controlled by the `output_format` parameter on `ConstructorParams`:

| Value      | Description                                             |
|:-----------|:--------------------------------------------------------|
| `:netcdf`  | **(default)** Single CF-compliant NetCDF-4 file (`wavi_input.nc`) |
| `:bin`     | Raw Float64 binary files (one `.bin` per variable)      |
| `:both`    | Write both formats side by side                         |

### NetCDF output

The default `:netcdf` format produces a single `wavi_input.nc` file containing all
variables with full CF-1.8 metadata:

- **Coordinate variables**: `x`, `y` (H-grid), `xu`/`yu` (U-grid), `xv`/`yv` (V-grid), `sigma` (depth)
- **Per-variable attributes**: `long_name`, `units`, `standard_name`, `_FillValue`
- **Global attributes**: grid spacing, domain origin, projection (EPSG:3031), basin IDs, creation timestamp

H-grid, U-grid (staggered in x), and V-grid (staggered in y) variables are stored
on their own dimension pairs so sizes are always correct.

### Binary output

The `:bin` format writes one file per variable as raw column-major Float64 arrays:

- `thickness.bin`, `surface.bin`, `bed.bin` — ice geometry
- `h_mask.bin` — domain mask
- `udata.bin`, `vdata.bin` — velocity components
- `udata_mask.bin`, `vdata_mask.bin` — velocity data masks
- `u_iszero.bin`, `v_iszero.bin` — velocity zero masks
- `accumulation_data.bin`, `dhdt_data.bin` — forcing fields
- `dhdt_acc_mask.bin` — combined data mask
- `basinID.bin` — basin identifiers
- `temps.bin` — 3-D temperature field (nx × ny × nz)
- `sigma_grid.bin` — sigma coordinate vector

### Choosing the format

```julia
# NetCDF (default)
params = default_constructor_params(basin_ids = [4])

# Binary
params = default_constructor_params(basin_ids = [4], output_format = :bin)

# Both
params = default_constructor_params(basin_ids = [4], output_format = :both)

# Or override at call time:
Gh, Gu, Gv, Gc = setup_wavi_data(params; output_format = :netcdf)
```

All output files are in the format expected by [WAVI.jl](https://github.com/WAVI-ice-sheet-model/WAVI.jl).

## Data Source Configuration

Each data category is specified as a `SourceConfig{S}` pairing a source type with a file path.

### Data source fields on `ConstructorParams`

| Field          | Category              | Required? | Default source          |
|:---------------|:----------------------|:----------|:------------------------|
| `bed`          | Bed topography        | yes       | `BedMachineV3()`        |
| `surface_temp` | Surface temperature   | yes       | `ALBMAPv1()`            |
| `temperature`  | 3-D temperature       | yes       | `FrankTemps()`          |
| `velocity`     | Ice velocity          | no        | `MEaSUREs()`            |
| `accumulation` | Snow accumulation     | no        | `ArthernAccumulation()` |
| `dhdt`         | Elevation change rate | no        | `SmithDhdt()`           |
| `basins`       | Drainage basins       | no        | `ZwallyBasins()`        |

### Scalar / processing parameters

| Parameter          | Description                      | Default      |
|:-------------------|:---------------------------------|:-------------|
| `dx`               | Grid spacing (m)                 | `10000.0`    |
| `basin_ids`        | Zwally basin IDs to include      | `1:27`       |
| `output_path`      | Directory for output files       | `"wavi_input"` |
| `output_format`    | Output format (`:bin`, `:netcdf`, `:both`) | `:netcdf` |
| `clip_edge_padding`| Edge padding when clipping       | `3`          |
| `density_ice`      | Ice density (kg/m³)              | `918.0`      |
| `density_ocean`    | Ocean density (kg/m³)            | `1028.0`     |
| `min_thick`        | Minimum ice thickness (m)        | `50.0`       |
| `sub_samp`         | Velocity subsampling factor      | `8`          |

### Adding a new data source

The dispatch architecture makes it easy to add new datasets:

1. Define a singleton: `struct MyNewSource <: TemperatureSource end`
2. Set its default path: `default_path(::MyNewSource) = "Data/my_file.nc"`
3. Implement `load_data(::MyNewSource, file)` returning the category's standard NamedTuple
4. Implement `interpolate_temperature(::MyNewSource, data, Gh)` (temperature sources only)

No `if/else` changes anywhere — the dispatch system routes automatically.
Users can then write `temperature = MyNewSource()` and the auto-wrapping handles the rest.

## Obtaining Data
TBC

## Contributing
TBC

## License
This software is licensed under an MIT license.

For more information please see the attached `LICENSE` file.

## Developers
WAVI.jl Ice Sheet Model Team, British Antarctic Survey
