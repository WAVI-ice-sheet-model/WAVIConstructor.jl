# WAVIConstructor.jl - Parameter Construction Helpers
# Helper functions to create parameter dictionaries for WAVI data construction

module ParamHelpers

using WAVIConstructor.DataSources

export ConstructorParams, default_constructor_params, minimal_constructor_params, to_dict

"""
    ConstructorParams

Struct defining all parameters needed for WAVI data construction.

Each geophysical dataset is specified as a [`SourceConfig`](@ref) that pairs a
*source type* singleton (e.g. `BISICLESTemps()`) with a file/directory path.
Use `NoData()` for optional categories you want to skip.

# Data-source fields

| Field          | Category              | Required? | Default source          |
|:---------------|:----------------------|:----------|:------------------------|
| `bed`          | Bed topography        | yes       | `BedMachineV3()`        |
| `surface_temp` | Surface temperature   | yes       | `ALBMAPv1()`            |
| `temperature`  | 3-D temperature       | yes       | `FrankTemps()`          |
| `velocity`     | Ice velocity          | no        | `MEaSUREs()`            |
| `accumulation` | Snow accumulation     | no        | `ArthernAccumulation()` |
| `dhdt`         | Elevation change rate | no        | `SmithDhdt()`           |
| `basins`       | Drainage basins       | no        | `ZwallyBasins()`        |

# Other fields

## Grid Parameters
- `dx::Float64`: Grid spacing in metres
- `basin_ids`: Basin IDs to include (can be range or array)

## Output Parameters
- `output_path::String`: Output directory for data files
- `output_format::Symbol`: Output format — `:bin`, `:netcdf`, or `:both`
- `clip_edge_padding::Int`: Padding around clipped domain

## Physical Constants
- `density_ice::Float64`: Ice density in kg/m³
- `density_ocean::Float64`: Ocean density in kg/m³
- `min_thick::Float64`: Minimum ice thickness in m

## Processing Parameters
- `sub_samp::Int`: Velocity subsampling factor
"""
Base.@kwdef struct ConstructorParams
    # ── Data sources (SourceConfig per category) ─────────────────────
    bed::SourceConfig{<:Union{BedSource,NoData}}              = SourceConfig(BedMachineV3())
    surface_temp::SourceConfig{<:Union{SurfaceTempSource,NoData}}    = SourceConfig(ALBMAPv1())
    temperature::SourceConfig{<:Union{TemperatureSource,NoData}} = SourceConfig(FrankTemps())
    velocity::SourceConfig{<:Union{VelocitySource,NoData}}    = SourceConfig(MEaSUREs())
    accumulation::SourceConfig{<:Union{AccumulationSource,NoData}} = SourceConfig(ArthernAccumulation())
    dhdt::SourceConfig{<:Union{DhDtSource,NoData}}            = SourceConfig(SmithDhdt())
    basins::SourceConfig{<:Union{BasinSource,NoData}}         = SourceConfig(ZwallyBasins())

    # ── Grid parameters ──────────────────────────────────────────────
    dx::Float64 = 10000.0
    basin_ids::Union{UnitRange{Int}, Vector{Int}, AbstractVector{Int}} = 1:27

    # ── Output parameters ────────────────────────────────────────────
    output_path::String = "wavi_input"
    output_format::Symbol = :bin          # :bin, :netcdf, or :both
    clip_edge_padding::Int = 3

    # ── Physical constants ───────────────────────────────────────────
    density_ice::Float64 = 918.0
    density_ocean::Float64 = 1028.0
    min_thick::Float64 = 50.0

    # ── Processing parameters ────────────────────────────────────────
    sub_samp::Int = 8
end

"""
    to_dict(params::ConstructorParams)

Convert a `ConstructorParams` struct to a `Dict{Symbol, Any}` for use with
`setup_wavi_data`.

Each `SourceConfig` field is unpacked into two keys:
  `:*_source` (the singleton) and `:*_file` (the path).
The `:basins` key maps to `basin_ids` (the ID vector/range),
while the basin *data source* becomes `:basins_source` / `:basins_file`.

# Example
```julia
params = ConstructorParams(dx = 5000.0)
d = to_dict(params)
Gh, Gu, Gv, Gc = setup_wavi_data(d)
```
"""
function to_dict(params::ConstructorParams)
    return Dict{Symbol, Any}(
        # Data sources (source singleton + path)
        :bed_source          => params.bed.source,
        :bed_file            => params.bed.path,
        :surface_temp_source  => params.surface_temp.source,
        :surface_temp_file    => params.surface_temp.path,
        :temperature_source  => params.temperature.source,
        :temperature_file    => params.temperature.path,
        :velocity_source     => params.velocity.source,
        :velocity_file       => params.velocity.path,
        :accumulation_source => params.accumulation.source,
        :accumulation_file   => params.accumulation.path,
        :dhdt_source         => params.dhdt.source,
        :dhdt_file           => params.dhdt.path,
        :basins_source       => params.basins.source,
        :basins_file         => params.basins.path,

        # Scalar / non-source fields
        :dx              => params.dx,
        :basins          => params.basin_ids,
        :output_path     => params.output_path,
        :output_format   => params.output_format,
        :clip_edge_padding => params.clip_edge_padding,
        :density_ice     => params.density_ice,
        :density_ocean   => params.density_ocean,
        :min_thick       => params.min_thick,
        :sub_samp        => params.sub_samp,
    )
end

# ── Auto-wrapping helper ─────────────────────────────────────────────

# Source-config fields on ConstructorParams that accept auto-wrapping.
const _SOURCE_FIELDS = (:bed, :surface_temp, :temperature, :velocity, :accumulation, :dhdt, :basins)

"""
    _as_source_config(x)

Normalise any accepted shorthand into a `SourceConfig`:
- `SourceConfig(…)`            → pass-through
- bare singleton `BedMachineV3()`     → `SourceConfig(BedMachineV3())`
- tuple `(BedMachineV3(), "path")`    → `SourceConfig(BedMachineV3(), "path")`
"""
_as_source_config(sc::SourceConfig) = sc
_as_source_config(s::DataSource)    = SourceConfig(s)
_as_source_config(t::Tuple{<:DataSource, <:AbstractString}) = SourceConfig(t[1], t[2])

function _wrap_source_kwargs(kwargs)
    processed = Dict{Symbol, Any}(pairs(kwargs))
    for f in _SOURCE_FIELDS
        if haskey(processed, f)
            processed[f] = _as_source_config(processed[f])
        end
    end
    return processed
end

"""
    default_constructor_params(; kwargs...)

Create default `ConstructorParams` for WAVI data construction.

Data-source fields accept any of the following shorthands:
- A bare singleton:  `bed = BedMachineV3()`  (uses its default path)
- A `(source, path)` tuple:  `temperature = (BISICLESTemps(), "my/temps.nc")`
- A full `SourceConfig`:  `bed = SourceConfig(BedMachineV3(), "path")`  (still works)

# Example
```julia
params = default_constructor_params(
    temperature = BISICLESTemps(),          # bare singleton — default path
    velocity    = (MEaSUREs(), "my/vel.nc"), # tuple — custom path
    dx          = 5000.0,
    basin_ids   = [1, 2, 3],
)
```
"""
function default_constructor_params(; kwargs...)
    return ConstructorParams(; _wrap_source_kwargs(kwargs)...)
end

"""
    minimal_constructor_params(; kwargs...)

Create minimal `ConstructorParams` for fast testing.
Uses coarse resolution and skips expensive datasets (velocity, dh/dt).

# Default Settings
- Grid spacing: 32 km (very coarse for speed)
- Basins: 0:27 (include basin 0 for missing Zwally data)
- dh/dt: `NoData` (skipped)
- Velocity: `NoData` (skipped)
- Temperature: `FrankTemps` (default)
- Output: `wavi_input_test`

# Example
```julia
params = minimal_constructor_params()
Gh, Gu, Gv, Gc = setup_wavi_data(params)
```
"""
function minimal_constructor_params(; kwargs...)
    defaults = Dict{Symbol, Any}(
        :dx          => 32000.0,
        :basin_ids   => 0:27,
        :dhdt        => NoData(),
        :velocity    => NoData(),
        :output_path => "wavi_input_test",
    )
    merge!(defaults, Dict{Symbol, Any}(pairs(kwargs)))
    return ConstructorParams(; _wrap_source_kwargs(defaults)...)
end

end # module ParamHelpers
