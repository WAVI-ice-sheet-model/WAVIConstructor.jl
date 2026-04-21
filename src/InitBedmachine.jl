# WAVIConstructor.jl - BedMachine v3 data initialization
# Julia port of initBedmachineV3.m

module InitBedMachine

using NCDatasets
using ArchGDAL
using Interpolations
using WAVI
using ProgressMeter

using WAVIConstructor.DataSources
using WAVIConstructor.DataLoading: load_data, interpolate_to_grid, interpolate_temperature

"""
    init_bedmachine(params)

Initialise computational grids and load all geophysical data using BedMachine v3.
Julia port of the MATLAB initBedmachineV3 function.

# Arguments
- `params`: NamedTuple with data loading parameters

# Returns
Tuple of grid structures (Gh, Gu, Gv, Gc) with loaded data
"""
function init_bedmachine(params)
    # ── Progress bar: 10 steps covering all major I/O and interpolation ──
    n_steps = 10
    prog = Progress(n_steps; desc="Initialising BedMachine: ", showspeed=true,
                    barglyphs=BarGlyphs("[=> ]"), color=:cyan)

    # Step 1 ─ Load BedMachine
    bed_source = get(params, :bed_source, BedMachineV3())
    bed_file   = get(params, :bed_file, default_path(BedMachineV3()))

    bm = load_data(bed_source, bed_file)
    next!(prog; showvalues = [("Step", "Loading BedMachine data")])
    bed       = bm.bed
    x         = bm.x
    y         = bm.y
    geoid     = bm.geoid
    mask      = bm.mask
    h         = bm.thickness
    s         = bm.surface

    # Flip arrays for correct orientation
    bed   = reverse(bed, dims=1)
    geoid = reverse(geoid, dims=1)
    mask  = reverse(mask, dims=1)
    s     = reverse(s, dims=1)
    h     = reverse(h, dims=1)

    # Create rock mask (mask=1 in BedMachine v3 indicates rock
    rockmask = zeros(size(mask))
    rockmask[mask .== 1] .= 1

    # Grid parameters
    nx_full = 13333
    ny_full = 13333
    x0_full = -3333000 - 250
    y0_full = -3333000 - 250
    dx_full = 500
    dy_full = 500

    x_full = x0_full .+ (0.5:(nx_full-0.5)) .* dx_full
    y_full = y0_full .+ (0.5:(ny_full-0.5)) .* dy_full

    i_pole = floor(Int, -x0_full / dx_full + 0.5)
    j_pole = floor(Int, -y0_full / dy_full + 0.5)


    # Subsampling
    sub_samp = get(params, :sub_samp, 8)
    sub_samp_index_x = get(params, :sub_samp_index_x, 0)
    sub_samp_index_y = get(params, :sub_samp_index_y, 0)

    # Domain size parameters (these may need adjustment for ALBMAP)
    domain_half_width = 2819  # Half-width in grid points
    domain_half_height = 2419  # Half-height in grid points

    isub = i_pole .+ 2*sub_samp .* (-floor(Int, (domain_half_width*2)/(2*sub_samp)):floor(Int, (domain_half_width*2)/(2*sub_samp)))
    jsub = j_pole .+ 2*sub_samp .* (-floor(Int, (domain_half_height*2)/(2*sub_samp)):floor(Int, (domain_half_height*2)/(2*sub_samp)))

    # Subsample all arrays
    bed = bed[isub, jsub]
    h = h[isub, jsub]
    s = s[isub, jsub]
    rockmask = rockmask[isub, jsub]
    geoid = geoid[isub, jsub]
    mask = mask[isub, jsub]

    x = x_full[isub]
    y = y_full[jsub]

    # Step 2 ─ Create H-grid structure
    Gh = create_h_grid(x, y, bed, h, s, geoid, rockmask, mask, params)
    next!(prog; showvalues = [("Step", "Building H-grid")])

    # Steps 3-5 ─ Load additional datasets (accumulation, basins, surface temp)
    Gh = load_additional_datasets(Gh, params, prog)
    
    # Calculate height above floatation and grounded ice mask
    density_ice = get(params, :density_ice, 918.0)
    density_ocean = get(params, :density_ocean, 1028.0)
    min_thick = get(params, :min_thick, 50.0)
    
    # Calculate hAF first
    hAF = Gh.h .+ (density_ocean / density_ice) .* Gh.b
    aground = (hAF .>= 0) .& .!Gh.rockmask .& (mask .!= 4)
    
    Gh = merge(Gh, (
        hAF = hAF,
        aground = aground,
        surfType = zeros(size(mask))
    ))
    
    # Set surfType based on aground mask
    surfType_calc = (Gh.aground .* 1) .+ ((.!Gh.aground .& (Gh.h .> 0.0)) .* 2) .+ (Gh.rockmask .* 3)
    Gh = merge(Gh, (surfType = surfType_calc,))
    
    # Recalculate ice and rock with updated surfType
    Gh = merge(Gh, (
        ice = (Gh.h .> min_thick) .& ((Gh.surfType .== 1) .| (Gh.surfType .== 2)),
        rock = (Gh.surfType .== 3) .| (Gh.aground .& (Gh.h .< min_thick)),
        ok = (Gh.h .> min_thick) .& ((Gh.surfType .== 1) .| (Gh.surfType .== 2))
    ))
    next!(prog; showvalues = [("Step", "Computing surface types & grounding")])
    
    # Step 7 ─ Load dhdt data via dispatch
    dhdt_source = get(params, :dhdt_source, SmithDhdt())
    dhdt_file   = get(params, :dhdt_file, default_path(SmithDhdt()))
    smith_data  = load_data(dhdt_source, dhdt_file)

    if smith_data !== nothing
        # Interpolate grounded data
        fdhdt_grnd = .!isnan.(smith_data.grnd_dhdt)
        dhdt_grnd = interpolate_to_grid(
            smith_data.grnd_xx[fdhdt_grnd][:], 
            smith_data.grnd_yy[fdhdt_grnd][:], 
            smith_data.grnd_dhdt[fdhdt_grnd][:], 
            Gh.xx, Gh.yy
        )
        
        # Interpolate floating data
        fdhdt_flt = .!isnan.(smith_data.flt_dhdt)
        dhdt_flt = interpolate_to_grid(
            smith_data.flt_xx[fdhdt_flt][:], 
            smith_data.flt_yy[fdhdt_flt][:], 
            smith_data.flt_dhdt[fdhdt_flt][:], 
            Gh.xx, Gh.yy
        )
        
        # Combine based on surfType: floating (2) uses flt, grounded (1 or 4) uses grnd
        dhdt = zeros(size(Gh.surfType))
        dhdt[Gh.surfType .== 2] .= dhdt_flt[Gh.surfType .== 2]
        dhdt[(Gh.surfType .== 1) .| (Gh.surfType .== 4)] .= dhdt_grnd[(Gh.surfType .== 1) .| (Gh.surfType .== 4)]
        
        Gh = merge(Gh, (dhdt = dhdt,))
    else
        # NoData or files not found — use zeros
        Gh = merge(Gh, (dhdt = zeros(size(Gh.surfType)),))
    end
    
    # Add aliases for compatibility with SelectDomainWAVI
    Gh = merge(Gh, (
        basinID = Gh.basin_id,
        a = Gh.a_Arthern,
        rock = Gh.rockmask
    ))
    next!(prog; showvalues = [("Step", "Loading dh/dt data")])

    # Step 8 ─ Create U/V/C grids
    Gu = create_u_grid(Gh)
    Gv = create_v_grid(Gh)
    Gc = create_c_grid(Gh)
    next!(prog; showvalues = [("Step", "Building U/V/C grids")])

    # Step 9 ─ Load velocity data
    Gu, Gv = load_velocity_data(Gu, Gv, Gh, params)
    next!(prog; showvalues = [("Step", "Loading velocity data")])

    # Step 10 ─ Load temperature data
    Gh = load_temperature_data(Gh, params)
    next!(prog; showvalues = [("Step", "Loading temperature data")])

    finish!(prog)
    return Gh, Gu, Gv, Gc
end

"""
    create_h_grid(x, y, bed, h, s, geoid, rockmask, mask, params)

Create the H-grid (scalar fields) structure.
"""
function create_h_grid(x, y, bed, h, s, geoid, rockmask, mask, params)
    Gh = (
        dx = x[2] - x[1],
        dy = y[2] - y[1],
        nx = length(x),
        ny = length(y),
        x0 = x[1] - 0.5 * (x[2] - x[1]),
        y0 = y[1] - 0.5 * (y[2] - y[1]),
        xx = [x[i] for i in 1:length(x), j in 1:length(y)],
        yy = [y[j] for i in 1:length(x), j in 1:length(y)],
        geoid = geoid,
        b = bed,
        h = h,
        s = s,
        mask = mask .> 0,  # Convert to boolean
        rockmask = rockmask .> 0
    )

    # Create coordinate grids
    Gh = merge(Gh, (xx = repeat(x, 1, length(y)), yy = repeat(y', length(x), 1)))

    return Gh
end

"""
    load_additional_datasets(Gh, params, prog)

Load accumulation, basin, and surface temperature data via dispatch on source types.
Advances `prog` once per dataset loaded (3 steps).
"""
function load_additional_datasets(Gh, params, prog)
    # ── Step 3 ─ Accumulation ─────────────────────────────────────────
    acc_source = get(params, :accumulation_source, ArthernAccumulation())
    acc_file   = get(params, :accumulation_file, default_path(ArthernAccumulation()))
    acc_data   = load_data(acc_source, acc_file)

    if acc_data !== nothing
        Gh = merge(Gh, (a_Arthern = interpolate_to_grid(acc_data.x, acc_data.y, acc_data.acc, Gh.xx, Gh.yy),))
    else
        @warn "Accumulation data skipped (source: $(typeof(acc_source))). Using zeros."
        Gh = merge(Gh, (a_Arthern = zeros(size(Gh.xx)),))
    end
    next!(prog; showvalues = [("Step", "Loading accumulation data")])

    # ── Step 4 ─ Drainage basins ──────────────────────────────────────
    basins_source = get(params, :basins_source, ZwallyBasins())
    basins_file   = get(params, :basins_file, default_path(ZwallyBasins()))
    basins_data   = load_data(basins_source, basins_file)

    if basins_data !== nothing
        Gh = merge(Gh, (basin_id = interpolate_to_grid(basins_data.xx, basins_data.yy, basins_data.basins, Gh.xx, Gh.yy),))
    else
        @warn "Basin data skipped (source: $(typeof(basins_source))). Using zeros."
        Gh = merge(Gh, (basin_id = zeros(size(Gh.xx)),))
    end
    next!(prog; showvalues = [("Step", "Loading drainage basins")])

    # ── Step 5 ─ Surface temperature / mean-annual temperature (ALBMAP)
    geom_source = get(params, :surface_temp_source, ALBMAPv1())
    geom_file   = get(params, :surface_temp_file, default_path(ALBMAPv1()))

    if geom_source isa NoData
        error("Surface temperature source (e.g. ALBMAP) is required and cannot be NoData.")
    end

    geom_data = load_data(geom_source, geom_file)
    Gh = merge(Gh, (Tma = interpolate_to_grid(geom_data.xx[:], geom_data.yy[:], geom_data.Tma[:], Gh.xx, Gh.yy),))
    next!(prog; showvalues = [("Step", "Loading surface temperature")])

    return Gh
end

"""
    create_u_grid(Gh)

Create U-grid (u-velocity components) structure.
"""
function create_u_grid(Gh)
    Gu = (
        x0 = Gh.x0,
        y0 = Gh.y0,
        dx = Gh.dx,
        dy = Gh.dy,
        nx = Gh.nx + 1,
        ny = Gh.ny,
        xx = [Gh.x0 + (i-1)*Gh.dx for i in 1:(Gh.nx+1), j in 1:Gh.ny],
        yy = [Gh.y0 + (j-0.5)*Gh.dy for i in 1:(Gh.nx+1), j in 1:Gh.ny]
    )
    return Gu
end

"""
    create_v_grid(Gh)

Create V-grid (v-velocity components) structure.
"""
function create_v_grid(Gh)
    Gv = (
        x0 = Gh.x0,
        y0 = Gh.y0,
        dx = Gh.dx,
        dy = Gh.dy,
        nx = Gh.nx,
        ny = Gh.ny + 1,
        xx = [Gh.x0 + (i-0.5)*Gh.dx for i in 1:Gh.nx, j in 1:(Gh.ny+1)],
        yy = [Gh.y0 + (j-1)*Gh.dy for i in 1:Gh.nx, j in 1:(Gh.ny+1)]
    )
    return Gv
end

"""
    create_c_grid(Gh)

Create C-grid (corner points) structure.
"""
function create_c_grid(Gh)
    Gc = (
        x0 = Gh.x0,
        y0 = Gh.y0,
        dx = Gh.dx,
        dy = Gh.dy,
        nx = Gh.nx - 1,
        ny = Gh.ny - 1,
        xx = [Gh.x0 + i*Gh.dx for i in 1:(Gh.nx-1), j in 1:(Gh.ny-1)],
        yy = [Gh.y0 + j*Gh.dy for i in 1:(Gh.nx-1), j in 1:(Gh.ny-1)]
    )
    return Gc
end

"""
    load_velocity_data(Gu, Gv, Gh, params)

Load ice velocity data via dispatch on the velocity source type.
"""
function load_velocity_data(Gu, Gv, Gh, params)
    sub_samp = get(params, :sub_samp, 8)

    vel_source = get(params, :velocity_source, MEaSUREs())
    vel_file   = get(params, :velocity_file, default_path(MEaSUREs()))
    vel_data   = load_data(vel_source, vel_file)

    if vel_data !== nothing
        xx_v = vel_data.xx
        yy_v = vel_data.yy
        vx   = vel_data.vx
        vy   = vel_data.vy
    else
        @warn "Velocity data skipped (source: $(typeof(vel_source))). Using zeros."
        xx_v = Gh.xx
        yy_v = Gh.yy
        vx = zeros(size(Gh.xx))
        vy = zeros(size(Gh.yy))
    end

    # Apply subsampling
    if sub_samp < 1.0
        xx_v = xx_v[1:1:end, 1:1:end]
        yy_v = yy_v[1:1:end, 1:1:end]
        vx = vx[1:1:end, 1:1:end]
        vy = vy[1:1:end, 1:1:end]
    else
        xx_v = xx_v[1:sub_samp:end, 1:sub_samp:end]
        yy_v = yy_v[1:sub_samp:end, 1:sub_samp:end]
        vx = vx[1:sub_samp:end, 1:sub_samp:end]
        vy = vy[1:sub_samp:end, 1:sub_samp:end]
    end

    u_data = interpolate_to_grid(xx_v[:], yy_v[:], vx[:], Gu.xx, Gu.yy)
    v_data = interpolate_to_grid(xx_v[:], yy_v[:], vy[:], Gv.xx, Gv.yy)
    
    Gu = merge(Gu, (
        u_data = u_data,
        u_data_mask = .!isnan.(u_data),
        uData = u_data,  # Alias for compatibility
        uDataMask = .!isnan.(u_data)
    ))

    Gv = merge(Gv, (
        v_data = v_data,
        v_data_mask = .!isnan.(v_data),
        vData = v_data,  # Alias for compatibility
        vDataMask = .!isnan.(v_data)
    ))

    return Gu, Gv
end

"""
    load_temperature_data(Gh, params)

Load and interpolate 3-D temperature data onto the H-grid.

The heavy lifting is fully dispatched:
- `load_data(source, file)`             → raw `(temps, xx, yy, sigma)` NamedTuple
- `interpolate_temperature(source, …)`  → `(temperature, sigmas)` on the H-grid

Adding a new temperature source requires only those two methods —
this orchestrator never needs to change.
"""
function load_temperature_data(Gh, params)
    temp_source = get(params, :temperature_source, FrankTemps())
    temp_file   = get(params, :temperature_file, default_path(FrankTemps()))

    if temp_source isa NoData
        error("Temperature source is required and cannot be NoData.")
    end

    temp_data = load_data(temp_source, temp_file)

    # Source-specific interpolation is fully dispatched — adding a new
    # temperature source only requires a `load_data` and an
    # `interpolate_temperature` method; no changes here.
    temperature, sigmas = interpolate_temperature(temp_source, temp_data, Gh)

    Gh = merge(Gh, (levels = (temperature = temperature, sigmas = sigmas),))
    return Gh
end

end # module InitBedMachine
