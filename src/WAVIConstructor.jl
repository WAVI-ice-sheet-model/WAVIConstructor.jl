module WAVIConstructor

include("DataLoading.jl")
include("InitBedmachine.jl")
include("DomainSelection.jl")
include("SetupData.jl")
include("ConstructorParams.jl")

# Re-export main functions from submodules
using .InitBedMachine: init_bedmachine
using .DomainSelection: select_domain_wavi
using .SetupData: setup_wavi_data
using .ParamHelpers: ConstructorParams, default_constructor_params, minimal_constructor_params, to_dict

export init_bedmachine, select_domain_wavi, setup_wavi_data
export ConstructorParams, default_constructor_params, minimal_constructor_params, to_dict

end # module
