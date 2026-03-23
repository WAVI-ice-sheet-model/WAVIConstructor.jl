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

| Parameter          | Description                      | Default    |
|:-------------------|:---------------------------------|:-----------|
| `dx`               | Grid spacing (m)                 | `10000.0`  |
| `basin_ids`        | Zwally basin IDs to include      | `1:27`     |
| `output_path`      | Directory for output files       | `"wavi_input"` |
| `clip_edge_padding`| Edge padding when clipping       | `3`        |
| `density_ice`      | Ice density (kg/m³)              | `918.0`    |
| `density_ocean`    | Ocean density (kg/m³)            | `1028.0`   |
| `min_thick`        | Minimum ice thickness (m)        | `50.0`     |
| `sub_samp`         | Velocity subsampling factor      | `8`        |

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
