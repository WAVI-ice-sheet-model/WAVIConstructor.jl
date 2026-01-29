# Basic Usage Example for WAVIConstructor.jl
# This script demonstrates how to use WAVIConstructor to prepare WAVI input data
#
# ⚠️  WARNING: This script is COMPUTATIONALLY EXPENSIVE! ⚠️
# - It processes the full BedMachine Antarctica dataset
# - It performs multiple large interpolations
# - It may take 5-30 minutes and use several GB of RAM
# - Recommended: Run on a compute node, not interactively
#
# For a quick test, use examples/minimal_test.jl instead!

using WAVIConstructor

# ============================================================================
# Example 1: Complete workflow using setup_wavi_data
# ============================================================================
println("Example 1: Complete workflow")
println("=" ^ 50)

# Use default_constructor_params with specific basins
params = default_constructor_params(
    dx = 10000.0,  # 10 km grid spacing
    basins = [1, 2, 3],  # Pine Island, Thwaites, and Haynes Glaciers region
    output_path = "wavi_input"
)

# This single function call does everything:
# 1. Initializes grids from BedMachine
# 2. Loads additional datasets (velocity, temperature, accumulation, basins)
# 3. Selects the domain based on specified basins
# 4. Clips the domain to remove unused regions
# 5. Replaces NaNs with appropriate values
# 6. Writes binary files for WAVI

# UNCOMMENT TO RUN (WARNING: COMPUTATIONALLY EXPENSIVE!)
# Gh, Gu, Gv, Gc = setup_wavi_data(params)

println("⚠️  Example 1 is commented out - it's very expensive!")
println("   Uncomment the setup_wavi_data call to run it.")
println("   Expected runtime: 5-30 minutes, RAM usage: several GB")

# println("✓ WAVI input files written to: $(params[:output_path])/")
# println("✓ Domain size (H-grid): $(Gh.nx) × $(Gh.ny)")
# println("✓ Number of ice points: $(Gh.n)")
println()

# ============================================================================
# Example 2: Step-by-step workflow with intermediate inspection
# ============================================================================
println("Example 2: Step-by-step workflow")
println("=" ^ 50)

# Step 1: Initialize grids
println("Step 1: Initializing grids from BedMachine...")
params_step = default_constructor_params(dx=10000.0)

# UNCOMMENT TO RUN (WARNING: COMPUTATIONALLY EXPENSIVE!)
# Gh, Gu, Gv, Gc = init_bedmachine(params_step)
# println("✓ Grids initialized: Gh ($(Gh.nx)×$(Gh.ny)), Gu ($(Gu.nx)×$(Gu.ny)), Gv ($(Gv.nx)×$(Gv.ny)), Gc ($(Gc.nx)×$(Gc.ny))")
# println("  Available fields in Gh: $(keys(Gh))")
# println()

# Step 2: Select domain
# println("Step 2: Selecting domain for basins 1-3...")
# params_step[:basins] = [1, 2, 3]
# Gh, Gu, Gv, Gc = select_domain_wavi(Gh, Gu, Gv, Gc, params_step)
# println("✓ Domain selected")
# println("  Ice points in domain (H-grid): $(Gh.n)")
# println("  Ice points in domain (U-grid): $(Gu.n)")
# println("  Ice points in domain (V-grid): $(Gv.n)")
# println("  Ice points in domain (C-grid): $(Gc.n)")
# println()

# At this point you could inspect the grids, plot data, etc.
# For example:
# using Plots
# heatmap(Gh.h', title="Ice Thickness (m)", xlabel="X", ylabel="Y")

# Step 3: Write output files
# println("Step 3: Processing and writing WAVI input files...")
# params_step[:output_path] = "wavi_input_step"
# params_step[:clip_edge_padding] = 3
# Gh, Gu, Gv, Gc = setup_wavi_data(params_step)
# println("✓ Files written to: $(params_step[:output_path])/")
println()

# ============================================================================
# Example 3: Using different datasets
# ============================================================================
println("Example 3: Using different datasets")
println("=" ^ 50)

# Use default_constructor_params for high-resolution setup
params_alt = default_constructor_params(
    dx = 5000.0,  # Higher resolution: 5 km
    basins = [10, 11, 12],  # Different basins
    veloc_data = "Measures_1",  # MEaSUREs from .mat file
    temps = "Frank",  # Frank's temperature data
    output_path = "wavi_input_alt"
)

println("Configuration:")
println("  Grid spacing: $(params_alt[:dx]) m")
println("  Basins: $(params_alt[:basins])")
println("  Velocity: $(params_alt[:veloc_data])")
println("  Temperature: $(params_alt[:temps])")
println()

# Uncomment to run:
# Gh, Gu, Gv, Gc = setup_wavi_data(params_alt)

# ============================================================================
# Example 4: All basins (full Antarctica)
# ============================================================================
println("Example 4: Full Antarctica setup")
println("=" ^ 50)

# Use default_constructor_params for full Antarctica
params_full = default_constructor_params(
    dx = 20000.0,  # Coarser grid: 20 km
    basins = 1:27,  # All Zwally basins
    output_path = "wavi_input_full"
)

println("This would process all of Antarctica at 20 km resolution")
println("Warning: This may take a while and produce large files!")
println()

# Uncomment to run:
# Gh, Gu, Gv, Gc = setup_wavi_data(params_full)

println("=" ^ 50)
println("Examples complete!")
println("=" ^ 50)
