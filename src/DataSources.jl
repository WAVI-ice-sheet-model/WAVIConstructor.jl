# WAVIConstructor.jl - Data Source Type Hierarchy
# Singleton types for dispatch-based data source selection.
#
# To add a new data source:
#   1. Define `struct MySource <: CategorySource end`
#   2. Define `default_path(::MySource) = "path/to/default/file"`
#   3. Implement `load_data(::MySource, file::String)` in DataLoading.jl
#      returning the category's standard NamedTuple shape.
#   No if/else changes anywhere.

module DataSources

export DataSource, SourceConfig, default_path
# Category abstract types
export BedSource, SurfaceTempSource, TemperatureSource, VelocitySource
export AccumulationSource, DhDtSource, BasinSource
# Concrete source singletons
export BedMachineV3, ALBMAPv1, FrankTemps, BISICLESTemps
export MEaSUREs, ArthernAccumulation, SmithDhdt, ZwallyBasins
# Cross-cutting sentinel
export NoData

# ── Abstract hierarchy ────────────────────────────────────────────────

"""
    DataSource

Root abstract type for all data sources.
Each data *category* (bed topography, temperature, velocity, …) has its own
abstract subtype, and each concrete dataset is a singleton struct under that
subtype.
"""
abstract type DataSource end

"""Bed topography / geometry sources (e.g. BedMachine)."""
abstract type BedSource       <: DataSource end

"""Mean annual surface temperature sources (e.g. ALBMAP)."""
abstract type SurfaceTempSource  <: DataSource end

"""3-D temperature field sources."""
abstract type TemperatureSource <: DataSource end

"""Ice velocity sources."""
abstract type VelocitySource  <: DataSource end

"""Snow/ice accumulation sources."""
abstract type AccumulationSource <: DataSource end

"""Ice elevation-change (dh/dt) sources."""
abstract type DhDtSource      <: DataSource end

"""Drainage-basin delineation sources."""
abstract type BasinSource     <: DataSource end

# ── Cross-cutting sentinel ────────────────────────────────────────────

"""
    NoData <: DataSource

Sentinel singleton indicating that no dataset should be loaded for an
optional data category.  `load_data(::NoData, …)` methods return zeros or
`nothing` as appropriate for the category.
"""
struct NoData <: DataSource end

"""
    default_path(::NoData)

NoData has no file path.
"""
default_path(::NoData) = ""

# ── Concrete singletons ──────────────────────────────────────────────

"""BedMachine Antarctica v3 bed topography."""
struct BedMachineV3 <: BedSource end
default_path(::BedMachineV3) = "Data/BedMachineAntarctica-v3.nc"

"""ALBMAP v1 Antarctic geometry / mean-annual temperature."""
struct ALBMAPv1 <: SurfaceTempSource end
default_path(::ALBMAPv1) = "Data/ALBMAPv1.nc"

"""Frank 3-D temperature field (.mat)."""
struct FrankTemps <: TemperatureSource end
default_path(::FrankTemps) = "Data/FranksTemps.mat"

"""BISICLES 3-D temperature field (NetCDF, 8 km)."""
struct BISICLESTemps <: TemperatureSource end
default_path(::BISICLESTemps) = "Data/antarctica-bisicles-xyzT-8km.nc"

"""MEaSUREs Antarctic ice velocity."""
struct MEaSUREs <: VelocitySource end
default_path(::MEaSUREs) = "Data/Antarctica_ice_velocity_2014_2015_1km_v01.nc"

"""Arthern AMSR accumulation map."""
struct ArthernAccumulation <: AccumulationSource end
default_path(::ArthernAccumulation) = "Data/amsr_accumulation_map.txt"

"""Smith et al. 2020 dh/dt elevation-change rates."""
struct SmithDhdt <: DhDtSource end
default_path(::SmithDhdt) = "Data/Smith_2020_dhdt"

"""Zwally drainage basins."""
struct ZwallyBasins <: BasinSource end
default_path(::ZwallyBasins) = "Data/ZwallyBasins.mat"

# ── SourceConfig wrapper ──────────────────────────────────────────────

"""
    SourceConfig{S<:DataSource}

Pairs a data-source singleton with the file/directory path to the dataset.

# Fields
- `source::S`  – singleton type identifying *which* dataset to load
- `path::String` – file or directory path to the data

# Examples
```julia
SourceConfig(BISICLESTemps())                                      # uses default_path
SourceConfig(BISICLESTemps(), "my/custom/path/temps.nc")           # custom path
SourceConfig(NoData())                                              # skip this category
```
"""
struct SourceConfig{S<:DataSource}
    source::S
    path::String
end

# Convenience: use the default path when only the source is given.
SourceConfig(s::S) where {S<:DataSource} = SourceConfig(s, default_path(s))

# ── Auto-conversion (bare singletons & tuples) ───────────────────────
# Allows writing `bed = BedMachineV3()` or `bed = (BedMachineV3(), "path")`
# instead of the verbose `bed = SourceConfig(BedMachineV3(), "path")`.

"""Auto-convert a bare `DataSource` singleton to a `SourceConfig` (uses `default_path`)."""
Base.convert(::Type{<:SourceConfig}, s::DataSource) = SourceConfig(s)

"""Auto-convert a `(source, path)` tuple to a `SourceConfig`."""
Base.convert(::Type{<:SourceConfig}, t::Tuple{<:DataSource, <:AbstractString}) = SourceConfig(t[1], t[2])

end # module DataSources
