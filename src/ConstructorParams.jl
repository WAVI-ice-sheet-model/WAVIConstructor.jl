# WAVIConstructor.jl - Parameter Construction Helpers
# Helper functions to create parameter dictionaries for WAVI data construction

module ParamHelpers

export ConstructorParams, default_constructor_params, minimal_constructor_params

"""
    ConstructorParams

Struct defining all parameters needed for WAVI data construction.

# Fields

## File Paths
- `bedmachine_file::String`: Path to BedMachine NetCDF file
- `smith_dhdt_dir::String`: Directory with Smith dhdt data
- `zwally_file::String`: Path to Zwally basins file
- `albmap_file::String`: Path to ALBMAP file
- `frank_temps_file::String`: Path to Frank temperatures file
- `measures_mat_file::String`: Path to MEaSUREs .mat file

## Grid Parameters
- `dx::Float64`: Grid spacing in meters
- `basins`: Basin IDs to include (can be range or array)

## Dataset Selection
- `dhdt_data::String`: dhdt dataset to use ("Smith" or "none")
- `veloc_data::String`: Velocity dataset ("Measures, or "none")
- `temps::String`: Temperature dataset ("BISICLES_8km", "BISICLES_1km", "Frank")

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
    # File paths
    bedmachine_file::String = "Data/BedMachineAntarctica-v3.nc"
    smith_dhdt_dir::String = "Data/Smith_2020_dhdt"
    zwally_file::String = "Data/DrainageBasins/ZwallyBasins.mat"
    albmap_file::String = "Data/ALBMAPv1.nc"
    frank_temps_file::String = "Data/FranksTemps.mat"
    measures_mat_file::String = "Data/MEaSUREs/MEaSUREsAntVels.mat"
    
    # Grid parameters
    dx::Float64 = 10000.0
    basins::Union{UnitRange{Int}, Vector{Int}, AbstractVector{Int}} = 1:27
    
    # Dataset selection
    dhdt_data::String = "Smith"
    veloc_data::String = "Measures_2016/17"
    temps::String = "BISICLES_8km"
    
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
        :frank_temps_file => params.frank_temps_file,
        :measures_mat_file => params.measures_mat_file,
        :dx => params.dx,
        :basins => params.basins,
        :dhdt_data => params.dhdt_data,
        :veloc_data => params.veloc_data,
        :temps => params.temps,
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
        dhdt_data = "none",
        veloc_data = "none",
        temps = "BISICLES_8km",
        output_path = "wavi_input_test",
        kwargs...
    )
end

end # module ParamHelpers
