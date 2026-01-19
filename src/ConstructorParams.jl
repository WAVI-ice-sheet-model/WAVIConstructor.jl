# WAVIConstructor.jl - Parameter Construction Helpers
# Helper functions to create parameter dictionaries for WAVI data construction

module ParamHelpers

export ConstructorParams, default_constructor_params, minimal_constructor_params, to_dict

"""
    ConstructorParams

Struct defining all parameters needed for WAVI data construction.

# Fields

## File Paths
All file paths should be full paths (absolute or relative to current working directory).

- `bedmachine_file::String`: Full path to BedMachine NetCDF file
- `smith_dhdt_dir::String`: Full path to directory with Smith dhdt data
- `zwally_file::String`: Full path to Zwally basins file
- `albmap_file::String`: Full path to ALBMAP file
- `arthern_file::String`: Full path to Arthern accumulation file
- `frank_temps_file::String`: Full path to Frank temperatures file (optional, alternative to bisicles_temps_file)
- `measures_velocity_file::String`: Full path to MEaSUREs velocity NetCDF file (optional)
- `bisicles_temps_file::String`: Full path to BISICLES temperature file (optional, alternative to frank_temps_file)

## Grid Parameters
- `dx::Float64`: Grid spacing in meters
- `basins`: Basin IDs to include (can be range or array)

## Output Parameters
- `output_path::String`: Output directory for binary files
- `clip_edge_padding::Int`: Padding around clipped domain

## Physical Constants
- `density_ice::Float64`: Ice density in kg/m³
- `density_ocean::Float64`: Ocean density in kg/m³
- `min_thick::Float64`: Minimum ice thickness in m

## Processing Parameters
- `sub_samp::Int`: Velocity subsampling factor
"""
Base.@kwdef struct ConstructorParams
    # File paths (all should be full paths, defaults assume Data/ folder)
    bedmachine_file::String = "Data/BedMachineAntarctica-v3.nc"
    smith_dhdt_dir::String = "Data/Smith_2020_dhdt"
    zwally_file::String = "Data/DrainageBasins/ZwallyBasins.mat"
    albmap_file::String = "Data/ALBMAPv1.nc"
    arthern_file::String = "Data/amsr_accumulation_map.txt"
    frank_temps_file::String = "Data/FranksTemps.mat"
    measures_velocity_file::String = "Data/antarctica_ice_velocity_2014_2015_1km_v01.nc"
    bisicles_temps_file::String = "Data/antarctica-bisicles-xyzT-8km.nc"
    
    # Grid parameters
    dx::Float64 = 10000.0
    basins::Union{UnitRange{Int}, Vector{Int}, AbstractVector{Int}} = 1:27
    
    # Output parameters
    output_path::String = "wavi_input"
    clip_edge_padding::Int = 3
    
    # Physical constants
    density_ice::Float64 = 918.0
    density_ocean::Float64 = 1028.0
    min_thick::Float64 = 50.0
    
    # Processing parameters
    sub_samp::Int = 8
end

"""
    to_dict(params::ConstructorParams)

Convert a ConstructorParams struct to a Dict{Symbol, Any} for use with setup_wavi_data.

# Example
```julia
params = ConstructorParams(dx=5000.0)
params_dict = to_dict(params)
Gh, Gu, Gv, Gc = setup_wavi_data(params_dict)
```
"""
function to_dict(params::ConstructorParams)
    return Dict{Symbol, Any}(
        :bedmachine_file => params.bedmachine_file,
        :smith_dhdt_dir => params.smith_dhdt_dir,
        :zwally_file => params.zwally_file,
        :albmap_file => params.albmap_file,
        :arthern_file => params.arthern_file,
        :frank_temps_file => params.frank_temps_file,
        :measures_velocity_file => params.measures_velocity_file,
        :bisicles_temps_file => params.bisicles_temps_file,
        :dx => params.dx,
        :basins => params.basins,
        :output_path => params.output_path,
        :clip_edge_padding => params.clip_edge_padding,
        :density_ice => params.density_ice,
        :density_ocean => params.density_ocean,
        :min_thick => params.min_thick,
        :sub_samp => params.sub_samp
    )
end

"""
    default_constructor_params(; kwargs...)

Create default ConstructorParams for WAVI data construction.

# Keyword Arguments
All fields from `ConstructorParams` struct can be overridden.

# Example
```julia
params = default_constructor_params(dx=5000.0, basins=[1,2,3])
params_dict = to_dict(params)
Gh, Gu, Gv, Gc = setup_wavi_data(params_dict)
```
"""
function default_constructor_params(; kwargs...)
    return ConstructorParams(; kwargs...)
end

"""
    minimal_constructor_params(; kwargs...)

Create minimal ConstructorParams for fast testing.
Uses coarse resolution and skips expensive datasets.

# Default Settings
- Grid spacing: 32 km (very coarse for speed)
- Basins: 0:27 (include basin 0 for missing Zwally data)
- dhdt: none (skipped)
- Velocity: none (skipped)
- Temperature: BISICLES_8km (required)
- Output: wavi_input_test

# Example
```julia
params = minimal_constructor_params()
params_dict = to_dict(params)
Gh, Gu, Gv, Gc = setup_wavi_data(params_dict)
```
"""
function minimal_constructor_params(; kwargs...)
    return ConstructorParams(;
        dx = 32000.0,
        basins = 0:27,
        smith_dhdt_dir = "",  # Empty string means skip dhdt
        measures_velocity_file = "",  # Empty string means skip velocity
        output_path = "wavi_input_test",
        kwargs...
    )
end

end # module ParamHelpers
