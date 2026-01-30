# Basin 19 Example
#
# This example demonstrates how to set up WAVI model inputs for Basin 19.
#
# Configuration:
#   - Resolution: 8 km (dx = 8000.0 m) - for initial testing
#   - Includes: MEaSUREs velocity data, Smith 2020 dh/dt, BISICLES temperatures
#   - Output: outputs/basin19_brunt_shelf/

using WAVIConstructor

# Get the package root directory (parent of examples/)
const PACKAGE_ROOT = dirname(@__DIR__)
const DATA_DIR = joinpath(PACKAGE_ROOT, "Data")
const OUTPUT_DIR = joinpath(PACKAGE_ROOT, "outputs", "basin19_brunt_shelf")

# Parameter Configuration

params = default_constructor_params(
    # BedMachine v3: bed topography, ice thickness, surface elevation
    bedmachine_file = joinpath(DATA_DIR, "BedMachineAntarctica-v3.nc"),
    
    # Zwally drainage basins (IDs 1-27)
    zwally_file = joinpath(DATA_DIR, "ZwallyBasins.mat"),
    
    # ALBMAP: mean annual temperature
    albmap_file = joinpath(DATA_DIR, "ALBMAPv1.nc"),
    
    # Arthern accumulation data
    arthern_file = joinpath(DATA_DIR, "amsr_accumulation_map.txt"),
    
    # BISICLES 3D temperature field (8 km resolution)
    bisicles_temps_file = joinpath(DATA_DIR, "antarctica-bisicles-xyzT-8km.nc"),
    
    # MEaSUREs ice velocity (1 km resolution)
    measures_velocity_file = joinpath(DATA_DIR, "Antarctica_ice_velocity_2014_2015_1km_v01.nc"),
    
    # Smith 2020 elevation change rates (dh/dt)
    smith_dhdt_dir = joinpath(DATA_DIR, "Smith_2020_dhdt"),
    
    # ----- Grid Settings -----
    # Grid spacing: 8 km for initial testing
    dx = 8000.0,
    
    # Basin selection: Basin 19 = Brunt Ice Shelf region
    basins = [19],
    
    # Velocity subsampling factor (8 for 8 km grid from 1 km source)
    sub_samp = 8,
    
    output_path = OUTPUT_DIR,
    
    # Edge padding for domain clipping (in grid cells)
    clip_edge_padding = 3,
    
    # Physical Parameters
    # Ice density (kg/m³)
    density_ice = 918.0,
    
    # Ocean density (kg/m³)
    density_ocean = 1028.0,
    
    # Minimum ice thickness (m) - thinner ice is masked out
    min_thick = 50.0
)

# Run WAVI Data Setup

# Execute the full workflow:
# 1. init_bedmachine: Load and interpolate BedMachine data to target grid
# 2. select_domain_wavi: Mask domain to selected basins, compute boundaries
# 3. Write binary output files for WAVI model
Gh, Gu, Gv, Gc = setup_wavi_data(to_dict(params); output_path=OUTPUT_DIR)

# Summary
println("Grid Configuration:")
println("  Grid spacing:     $(Gh.dx) m ($(Gh.dx / 1000) km)")
println("  Full grid size:   $(Gh.nx) × $(Gh.ny)")
println("  Clipped grid:     $(Gh.nx_clip) × $(Gh.ny_clip)")
println()
println("Domain Bounds (clipped):")
println("  Origin (x0, y0):  ($(Gh.x0_clip), $(Gh.y0_clip)) m")
println("  Extent:           $(Gh.nx_clip * Gh.dx / 1000) × $(Gh.ny_clip * Gh.dx / 1000) km")
println()
println("Ice Statistics:")
println("  Ice points:       $(Gh.n_clip)")
println()
println("Output: $(OUTPUT_DIR)/")

println("\nSetup complete!")