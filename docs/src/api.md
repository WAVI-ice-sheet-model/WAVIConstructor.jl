# API Reference

This page documents the main exported functions and types in WAVIConstructor.jl. The package provides a three-step workflow for preparing WAVI input data: initialize grids, select domain, and write output files.

## Main Workflow Functions

These are the core functions that implement the main workflow for preparing WAVI input data.

```@docs
WAVIConstructor.init_bedmachine
WAVIConstructor.select_domain_wavi
WAVIConstructor.setup_wavi_data
```

## Parameter Helpers

Functions and types for creating and managing parameter dictionaries for WAVI data construction.

```@docs
WAVIConstructor.ConstructorParams
WAVIConstructor.default_constructor_params
WAVIConstructor.minimal_constructor_params
WAVIConstructor.to_dict
```

## Data Loading Functions

Functions for loading various geophysical datasets used in WAVI simulations. These are exported from the `DataLoading` submodule and can be accessed directly or via `WAVIConstructor.DataLoading`.

```@docs
WAVIConstructor.DataLoading.get_albmap
WAVIConstructor.DataLoading.get_bedmachine
WAVIConstructor.DataLoading.get_bisicles_temps
WAVIConstructor.DataLoading.get_measures_velocities
WAVIConstructor.DataLoading.get_smith_dhdt
WAVIConstructor.DataLoading.get_arthern_accumulation
WAVIConstructor.DataLoading.get_zwally_basins
WAVIConstructor.DataLoading.get_frank_temps
WAVIConstructor.DataLoading.geotiff_read_axis_only
```

## Internal Functions

Additional internal functions and helpers are automatically documented below:

```@autodocs
Modules = [WAVIConstructor]
```
