# Minimal Test Example for WAVIConstructor.jl
# This is a lightweight test that won't crash your session
# It uses coarse resolution and skips expensive datasets

using WAVIConstructor

println("Minimal WAVIConstructor Test")
println("=" ^ 50)

# Use minimal_constructor_params helper for fast testing
# This uses coarse resolution (32 km) and skips expensive datasets
params = minimal_constructor_params(temps="Frank")

println("Configuration:")
println("  Grid spacing: $(params[:dx]) m")
println("  Basins: $(params[:basins])")
println("  Temperature: $(params[:temps])")
println("  Velocity: $(params[:veloc_data])")
println("  dhdt: $(params[:dhdt_data])")
println()

println("Running setup_wavi_data...")
println("(This should take ~30 seconds on the coarse grid)")
println()

try
    Gh, Gu, Gv, Gc = setup_wavi_data(params)
    
    println()
    println("✓ SUCCESS!")
    println("=" ^ 50)
    println("Output written to: $(params[:output_path])/")
    println("Domain size (H-grid): $(Gh.nx) × $(Gh.ny)")
    println("Number of ice points: $(Gh.n)")
    println()
    println("Files created:")
    if isdir(params[:output_path])
        for file in readdir(params[:output_path])
            println("  - $(file)")
        end
    else
        println("  (Could not list files - directory path issue)")
    end
    println()
    println("Next steps:")
    println("  1. Check the output files in $(params[:output_path])/")
    println("  2. Try examples/basic_usage.jl on a compute node with more resources")
    println("  3. Adjust :dx to balance resolution vs. computation time")
    
catch e
    println("ERROR: $e")
    println()
    println("If this failed, check:")
    println("  1. Does Data/BedMachineAntarctica-v3.nc exist?")
    println("  2. Do you have enough memory? (Try :dx => 64000.0 for even coarser grid)")
    rethrow(e)
end
