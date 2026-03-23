module WAVIConstructor

include("DataSources.jl")
include("DataLoading.jl")
include("ConstructorParams.jl")
include("InitBedmachine.jl")
include("DomainSelection.jl")
include("SetupData.jl")

# Re-export main functions from submodules
using .InitBedMachine: init_bedmachine
using .DomainSelection: select_domain_wavi
using .SetupData: setup_wavi_data
using .ParamHelpers: ConstructorParams, default_constructor_params, minimal_constructor_params, to_dict
using .DataSources: DataSource, SourceConfig, default_path, NoData,
    BedSource, GeometrySource, TemperatureSource, VelocitySource,
    AccumulationSource, DhDtSource, BasinSource,
    BedMachineV3, ALBMAPv1, FrankTemps, BISICLESTemps,
    MEaSUREs, ArthernAccumulation, SmithDhdt, ZwallyBasins
using .DataLoading: load_data, interpolate_to_grid, interpolate_temperature

export init_bedmachine, select_domain_wavi, setup_wavi_data
export ConstructorParams, default_constructor_params, minimal_constructor_params, to_dict
export DataSource, SourceConfig, default_path, NoData
export BedSource, GeometrySource, TemperatureSource, VelocitySource
export AccumulationSource, DhDtSource, BasinSource
export BedMachineV3, ALBMAPv1, FrankTemps, BISICLESTemps
export MEaSUREs, ArthernAccumulation, SmithDhdt, ZwallyBasins
export load_data, interpolate_to_grid, interpolate_temperature

end # module
