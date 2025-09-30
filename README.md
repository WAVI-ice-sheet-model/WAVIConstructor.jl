# WAVIConstructor.jl
![Work in Progress](https://img.shields.io/badge/status-work--in--progress-orange)
![Experimental](https://img.shields.io/badge/status-experimental-blueviolet)
![License: MIT](https://img.shields.io/badge/license-MIT-green)


Pre-processing library for [WAVI.jl](https://github.com/WAVI-ice-sheet-model/WAVI.jl).

This package is a WIP, but will support the reading and pre-processing of input data, and preparing a WAVI domain.

## Installation

You can install the latest version of `WAVIConstructor` using Julia's in-built package manager:
```julia
using Pkg
Pkg.add(PackageSpec(url="https://github.com/WAVI-ice-sheet-model/WAVIConstructor.jl.git", rev = "main"))
```

## Supported Input Datasets

This package supports loading and processing the following datasets (readers provided in `WAVIConstructor.DataLoading`):

- **ALBMAP** (NetCDF): Antarctic ice sheet surface elevation, bed topography, firn thickness, mean annual temperature, and accumulation rate fields (`get_albmap`).

- **BedMachine Antarctica v3** (NetCDF): Bed topography, ice surface elevation, ice thickness, geoid height, and ice/ocean/land mask (`get_bedmachine`).

- **BISICLES Temperature Data** (NetCDF): 3D temperature fields from BISICLES ice sheet model simulations, including sigma coordinates and spatial grids (`get_bisicles_temps`).

- **MEASUREs Ice Velocity** (NetCDF): NASA MEASUREs program ice velocity products providing annual ice motion vectors (VX, VY)(`get_measures_velocities`).

- **Smith et al. 2020 Elevation Change** (GeoTIFF): Ice surface elevation change (dh/dt) datasets from Smith et al. 2020 study, ice sheet mass change derived from ICESat/ICESat-2(`get_smith_dhdt`).

## Examples

See the `examples/` directory for sample scripts demonstrating how to load and use the supported datasets. The main example (`examples/loading_example.jl`) shows how to read all the data files in the `Data/` directory.
