module WAVIConstructor

include("DataLoading.jl")
include("InitBedmachine.jl")

# Re-export main functions from submodules
using .InitBedMachine: init_bedmachine
export init_bedmachine

end # module
