# Basin 4 Example
#
# This example demonstrates how to set up WAVI model inputs for Basin 4.
#
# Configuration:
#   - Resolution: 8 km (dx = 8000.0 m) - for initial testing
#   - Includes: MEaSUREs velocity data, Smith 2020 dh/dt, BISICLES temperatures
#   - Output: outputs/basin4_brunt_shelf/

using WAVIConstructor

# Parameter Configuration

params = default_constructor_params(
    # BedMachine v3: bed topography, ice thickness, surface elevation
    bed = (BedMachineV3(), "Data/BedMachineAntarctica-v3.nc"),

    # Zwally drainage basins (IDs 1-27)
    basins = (ZwallyBasins(), "Data/ZwallyBasins.mat"),

    # ALBMAP: mean annual surface temperature
    surface_temp = (ALBMAPv1(), "Data/ALBMAPv1.nc"),

    # Arthern accumulation data
    accumulation = (ArthernAccumulation(), "Data/amsr_accumulation_map.txt"),

    # BISICLES 3D temperature field (8 km resolution)
    temperature = (BISICLESTemps(), "Data/antarctica-bisicles-xyzT-8km.nc"),

    # MEaSUREs ice velocity (1 km resolution)
    velocity = (MEaSUREs(), "Data/Antarctica_ice_velocity_2014_2015_1km_v01.nc"),

    # Smith 2020 elevation change rates (dh/dt)
    dhdt = (SmithDhdt(), "Data/Smith_2020_dhdt"),

    # ----- Grid Settings -----
    # Grid spacing: 8 km for initial testing
    dx = 8000.0,

    # Basin selection: Basin 4 = Brunt Ice Shelf region
    basin_ids = [4],

    # Velocity subsampling factor (8 for 8 km grid from 1 km source)
    sub_samp = 8,

    output_path = "outputs/basin4_brunt_shelf",

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
Gh, Gu, Gv, Gc = setup_wavi_data(params; output_path="outputs/basin4_brunt_shelf")

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

println("\nSetup complete!")