module WAVIConstructor

include("DataLoading.jl")
include("InitBedmachine.jl")
include("DomainSelection.jl")
include("SetupData.jl")

# Re-export main functions from submodules
using .InitBedMachine: init_bedmachine
using .DomainSelection: select_domain_wavi
using .SetupData: setup_wavi_data

export init_bedmachine, select_domain_wavi, setup_wavi_data

end # module
