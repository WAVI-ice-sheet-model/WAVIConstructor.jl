# API Reference

This page documents the main exported functions and types in WAVIConstructor.jl. Each dataset has its own label type (e.g. `BedMachineV3`, `BISICLESTemps`), and passing that label to `load_data` which calls the right loading code for it.

## Main Workflow Functions

These are the core functions that implement the main workflow for preparing WAVI input data.

```@docs
WAVIConstructor.init_bedmachine
WAVIConstructor.select_domain_wavi
WAVIConstructor.setup_wavi_data
```

## Parameter Helpers

Functions and types for creating and managing the `ConstructorParams` configuration struct.

```@docs
WAVIConstructor.ConstructorParams
WAVIConstructor.default_constructor_params
WAVIConstructor.minimal_constructor_params
WAVIConstructor.to_dict
```

## Data Source Types

The type hierarchy used for dispatch-based data source selection.

```@docs
WAVIConstructor.DataSources.DataSource
WAVIConstructor.DataSources.SourceConfig
WAVIConstructor.DataSources.NoData
WAVIConstructor.DataSources.default_path
```

### Category abstract types

Each data category has its own abstract subtype of `DataSource`:

- `BedSource` — bed topography (e.g. `BedMachineV3`)
- `SurfaceTempSource` — mean annual surface temperature (e.g. `ALBMAPv1`)
- `TemperatureSource` — 3-D temperature fields (e.g. `FrankTemps`, `BISICLESTemps`)
- `VelocitySource` — ice velocity (e.g. `MEaSUREs`)
- `AccumulationSource` — snow accumulation (e.g. `ArthernAccumulation`)
- `DhDtSource` — elevation-change rates (e.g. `SmithDhdt`)
- `BasinSource` — drainage basin delineations (e.g. `ZwallyBasins`)

```@docs
WAVIConstructor.DataSources.BedSource
WAVIConstructor.DataSources.SurfaceTempSource
WAVIConstructor.DataSources.TemperatureSource
WAVIConstructor.DataSources.VelocitySource
WAVIConstructor.DataSources.AccumulationSource
WAVIConstructor.DataSources.DhDtSource
WAVIConstructor.DataSources.BasinSource
```

### Concrete source singletons

```@docs
WAVIConstructor.DataSources.BedMachineV3
WAVIConstructor.DataSources.ALBMAPv1
WAVIConstructor.DataSources.FrankTemps
WAVIConstructor.DataSources.BISICLESTemps
WAVIConstructor.DataSources.MEaSUREs
WAVIConstructor.DataSources.ArthernAccumulation
WAVIConstructor.DataSources.SmithDhdt
WAVIConstructor.DataSources.ZwallyBasins
```

## Data Loading Functions

All dataset loading goes through the dispatch-based `load_data` interface.
Each source singleton selects the appropriate method.

```@docs
WAVIConstructor.DataLoading.load_data
WAVIConstructor.DataLoading.interpolate_to_grid
WAVIConstructor.DataLoading.interpolate_temperature
WAVIConstructor.DataLoading.geotiff_read_axis_only
```

## Output Writing

Functions for writing WAVI model input data in binary and/or NetCDF format.

```@docs
WAVIConstructor.OutputWriting.write_output
WAVIConstructor.OutputWriting.write_binary_files
WAVIConstructor.OutputWriting.write_netcdf_file
```

## Internal Functions

Additional internal functions and helpers are automatically documented below:

```@autodocs
Modules = [WAVIConstructor]
```
