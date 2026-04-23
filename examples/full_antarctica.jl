# Full Antarctica Example — matching WAVIConstructor.m/get_setup_data.m
#
# This example replicates the MATLAB workflow in get_setup_data.m with
# identical parameter choices:
#   - BedMachine v3 bed topography
#   - MEaSUREs 2014/15 ice velocities
#   - BISICLES 8 km 3-D temperature field
#   - Smith et al. 2020 dh/dt
#   - Arthern accumulation
#   - ALBMAP surface temperature
#   - Zwally drainage basins (all basins, 1:27)
#   - 8 km grid resolution (subSamp = 8 in MATLAB → dx = 8000 m)
#
# Usage:
#   julia --project=. examples/full_antarctica.jl
#
# Output directory: outputs/full_antarctica/

using WAVIConstructor

# ──────────────────────────────────────────────────────────────────────
# Parameter configuration — mirrors DefaultParamsWAVI.m + get_setup_data.m
# ──────────────────────────────────────────────────────────────────────
params = default_constructor_params(
    # ── Data sources ──────────────────────────────────────────────────
    # BedMachine v3: bed, thickness, surface, geoid, mask
    bed = (BedMachineV3(), "Data/BedMachineAntarctica-v3.nc"),

    # Zwally drainage basins
    basins = (ZwallyBasins(), "Data/ZwallyBasins.mat"),

    # ALBMAP: mean annual surface temperature (Tma)
    surface_temp = (ALBMAPv1(), "Data/ALBMAPv1.nc"),

    # Arthern accumulation (AMSR)
    accumulation = (ArthernAccumulation(), "Data/amsr_accumulation_map.txt"),

    # BISICLES 8 km 3-D temperature field
    # (MATLAB: params.temps = 'BISICLES_8km')
    temperature = (BISICLESTemps(), "Data/antarctica-bisicles-xyzT-8km.nc"),

    # MEaSUREs ice velocity 2014/15 (1 km source resolution)
    # (MATLAB: params.veloc_Data = 'Measures_2014/15')
    velocity = (MEaSUREs(), "Data/Antarctica_ice_velocity_2014_2015_1km_v01.nc"),

    # Smith et al. 2020 elevation change rates (dh/dt)
    # (MATLAB: params.dhdt_data = 'Smith')
    dhdt = (SmithDhdt(), "Data/Smith_2020_dhdt"),

    # ── Grid settings ─────────────────────────────────────────────────
    # 8 km resolution (MATLAB: subSamp = 8 → dx = 500 * 2 * 8 = 8000 m)
    dx = 8000.0,

    # All Zwally basins (MATLAB: params.basins = (1:27))
    basin_ids = collect(1:27),

    # Velocity subsampling factor
    # (MATLAB: params.subSamp = 8, applied to 1 km MEaSUREs data)
    sub_samp = 8,

    # ── Domain clipping ───────────────────────────────────────────────
    # Edge padding around mask bounding box (MATLAB: edge = 3)
    clip_edge_padding = 3,

    # ── Physical constants (MATLAB DefaultParamsWAVI) ─────────────────
    density_ice   = 918.0,   # kg/m³
    density_ocean = 1028.0,  # kg/m³
    min_thick     = 50.0,    # m — ice thinner than this is masked out

    # ── Output ────────────────────────────────────────────────────────
    output_path   = "outputs/full_antarctica",
    output_format = :both,   # write both .bin and .nc files
)

# ──────────────────────────────────────────────────────────────────────
# Run the full workflow
# ──────────────────────────────────────────────────────────────────────
# Steps (mirrors MATLAB get_setup_data.m):
#   1. init_bedmachine  — load BedMachine, accumulation, basins, ALBMAP,
#                          temperatures, velocities, dh/dt; build grids
#   2. select_domain_wavi — apply basin mask, compute halo / rock edges
#   3. Clip domain, replace NaNs, extrapolate temperatures to sigma = 0/1
#   4. Write output files (.bin and/or .nc)
println("Running full Antarctica setup (8 km, basins 1:27)...")
println("This may take a while — loading and interpolating multiple datasets.\n")

Gh, Gu, Gv, Gc = setup_wavi_data(params; output_path="outputs/full_antarctica")

# ──────────────────────────────────────────────────────────────────────
# Summary
# ──────────────────────────────────────────────────────────────────────
println("\n", "=" ^ 60)
println("Full Antarctica Setup — Summary")
println("=" ^ 60)

println("\nGrid Configuration:")
println("  Grid spacing:     $(Gh.dx) m ($(Gh.dx / 1000) km)")
println("  Full grid size:   $(Gh.nx) × $(Gh.ny)")
println("  Clipped grid:     $(Gh.nx_clip) × $(Gh.ny_clip)")

println("\nDomain Bounds (clipped):")
println("  Origin (x0, y0):  ($(Gh.x0_clip), $(Gh.y0_clip)) m")
println("  Extent:           $(Gh.nx_clip * Gh.dx / 1000) × $(Gh.ny_clip * Gh.dx / 1000) km")

println("\nIce Statistics:")
println("  H-grid ice pts:   $(Gh.n_clip)")

println("\nPhysical Parameters (from MATLAB DefaultParamsWAVI):")
density_ice   = 918.0
density_ocean = 1028.0
delta = 1 - density_ice / density_ocean
println("  density_ice:       $density_ice kg/m³")
println("  density_ocean:     $density_ocean kg/m³")
println("  delta:             $delta")
println("  min_thick:         50.0 m")

println("\nOutput written to: outputs/full_antarctica/")
println("Setup complete!")
