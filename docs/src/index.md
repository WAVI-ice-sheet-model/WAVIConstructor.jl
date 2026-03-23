# WAVIConstructor.jl

Pre-processing library for [WAVI.jl](https://github.com/WAVI-ice-sheet-model/WAVI.jl).

WAVIConstructor reads geophysical datasets (bed topography, ice velocity, temperature, accumulation, elevation-change rates, drainage basins), interpolates them onto a common grid, selects a regional domain, and writes files ready for WAVI inversion runs.

## Architecture

The package uses a **singleton-type dispatch** pattern for data source selection.
Each dataset is identified by a lightweight singleton struct (e.g. `BedMachineV3()`, `BISICLESTemps()`), and all loading is routed through a single function:

```julia
data = load_data(BedMachineV3(), "Data/BedMachineAntarctica-v3.nc")
```

A [`SourceConfig`](@ref) wrapper pairs a source singleton with its file path.
The constructor helpers auto-wrap bare singletons and `(source, path)` tuples,
so you rarely need to write `SourceConfig(…)` explicitly:

```julia
# All equivalent for the temperature field:
temperature = BISICLESTemps()                                       # bare singleton (default path)
temperature = (BISICLESTemps(), "my/custom/temps.nc")               # tuple (custom path)
temperature = SourceConfig(BISICLESTemps(), "my/custom/temps.nc")   # explicit (still works)
temperature = NoData()                                              # skip this category
```

### Adding a new data source

1. Define a singleton: `struct MySource <: TemperatureSource end`
2. Set default path: `default_path(::MySource) = "Data/my_file.nc"`
3. Implement `load_data(::MySource, file)` returning the standard NamedTuple for that category
4. (Temperature sources only) Implement `interpolate_temperature(::MySource, data, Gh)`

No if/else chains — dispatch handles routing automatically.
Users just write `temperature = MySource()` and the auto-wrapping handles the rest.

## Quick Start

```julia
using WAVIConstructor

params = default_constructor_params(
    temperature = BISICLESTemps(),
    basin_ids   = [4],
    dx          = 8000.0,
)

Gh, Gu, Gv, Gc = setup_wavi_data(params)
```

See the [API Reference](@ref) for full details on all types and functions.
