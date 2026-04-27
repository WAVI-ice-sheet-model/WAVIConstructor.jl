# WAVIConstructor.jl - BedMachine v3 data initialization
# Julia port of initBedmachineV3.m

module InitBedMachine

using NCDatasets
using ArchGDAL
using Interpolations
using WAVI

using WAVIConstructor.DataSources
using WAVIConstructor.DataLoading: load_data, interpolate_to_grid, interpolate_regular_grid_nearest, interpolate_temperature

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
    bed_source = get(params, :bed_source, BedMachineV3())
    bed_file   = get(params, :bed_file, default_path(BedMachineV3()))

    bm = load_data(bed_source, bed_file)
    bed       = bm.bed
    x         = bm.x
    y         = bm.y
    geoid     = bm.geoid
    mask      = bm.mask
    h         = bm.thickness
    s         = bm.surface

    # Flip arrays to match MATLAB's fliplr (reverse along dimension 2 = columns/y-axis)
    bed   = reverse(bed, dims=2)
    geoid = reverse(geoid, dims=2)
    mask  = reverse(mask, dims=2)
    s     = reverse(s, dims=2)
    h     = reverse(h, dims=2)

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

    # Create H-grid structure
    Gh = create_h_grid(x, y, bed, h, s, geoid, rockmask, mask, params)
    Gh = load_additional_datasets(Gh, params)
    
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
    
    # Load dhdt data via dispatch
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
    # Gh.rock must match MATLAB: (surfType==3) | (grounded & h < min_thick)
    # NOT just rockmask (BedMachine mask==1), which omits thin grounded ice cells.
    Gh = merge(Gh, (
        basinID = Gh.basin_id,
        a = Gh.a_Arthern,
        rock = Gh.rock   # already set above as (surfType==3)|(aground & h<min_thick)
    ))

    # Create U/V/C grids
    Gu = create_u_grid(Gh)
    Gv = create_v_grid(Gh)
    Gc = create_c_grid(Gh)

    Gu, Gv = load_velocity_data(Gu, Gv, Gh, params)
    Gh = load_temperature_data(Gh, params)

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
    load_additional_datasets(Gh, params)

Load accumulation, basin, and surface temperature data via dispatch on source types.
"""
function load_additional_datasets(Gh, params)
    # ── Accumulation ──────────────────────────────────────────────────
    acc_source = get(params, :accumulation_source, ArthernAccumulation())
    acc_file   = get(params, :accumulation_file, default_path(ArthernAccumulation()))
    acc_data   = load_data(acc_source, acc_file)

    if acc_data !== nothing
        Gh = merge(Gh, (a_Arthern = interpolate_to_grid(acc_data.x, acc_data.y, acc_data.acc, Gh.xx, Gh.yy),))
    else
        @warn "Accumulation data skipped (source: $(typeof(acc_source))). Using zeros."
        Gh = merge(Gh, (a_Arthern = zeros(size(Gh.xx)),))
    end

    # ── Drainage basins ───────────────────────────────────────────────
    basins_source = get(params, :basins_source, ZwallyBasins())
    basins_file   = get(params, :basins_file, default_path(ZwallyBasins()))
    basins_data   = load_data(basins_source, basins_file)

    if basins_data !== nothing
        Gh = merge(Gh, (basin_id = interpolate_to_grid(basins_data.xx, basins_data.yy, basins_data.basins, Gh.xx, Gh.yy),))
    else
        @warn "Basin data skipped (source: $(typeof(basins_source))). Using zeros."
        Gh = merge(Gh, (basin_id = zeros(size(Gh.xx)),))
    end

    # ── Surface temperature / mean-annual temperature (ALBMAP) ────────
    geom_source = get(params, :surface_temp_source, ALBMAPv1())
    geom_file   = get(params, :surface_temp_file, default_path(ALBMAPv1()))

    if geom_source isa NoData
        error("Surface temperature source (e.g. ALBMAP) is required and cannot be NoData.")
    end

    geom_data = load_data(geom_source, geom_file)
    Gh = merge(Gh, (Tma = interpolate_to_grid(geom_data.xx[:], geom_data.yy[:], geom_data.Tma[:], Gh.xx, Gh.yy),))

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
        vx   = vel_data.vx
        vy   = vel_data.vy

        # Subsampled 1-D coordinate axes and 2-D value arrays
        step = sub_samp < 1.0 ? 1 : sub_samp
        xs  = vel_data.x[1:step:end]
        ys  = vel_data.y[1:step:end]
        vx_sub = vx[1:step:end, 1:step:end]
        vy_sub = vy[1:step:end, 1:step:end]

        # Use regular-grid nearest-neighbour with MATLAB-style round-half-up
        # tie-breaking (floor(f+0.5)) to replicate TriScatteredInterp behaviour.
        u_data = interpolate_regular_grid_nearest(xs, ys, vx_sub, Gu.xx, Gu.yy)
        v_data = interpolate_regular_grid_nearest(xs, ys, vy_sub, Gv.xx, Gv.yy)
        # Points outside the MEaSUREs bounding box get NaN.  Leave them as NaN
        # here so that replace_nans_in_clipped_velocity can write -9999 for
        # those cells, matching MATLAB's output.  The mask excludes both NaN
        # and zero (MEaSUREs fill value).
    else
        @warn "Velocity data skipped (source: $(typeof(vel_source))). Using zeros."
        u_data = zeros(size(Gu.xx))
        v_data = zeros(size(Gv.xx))
    end

    Gu = merge(Gu, (
        u_data = u_data,
        u_data_mask = .!isnan.(u_data) .& (u_data .!= 0.0),
        uData = u_data,
        uDataMask = .!isnan.(u_data) .& (u_data .!= 0.0)
    ))

    Gv = merge(Gv, (
        v_data = v_data,
        v_data_mask = .!isnan.(v_data) .& (v_data .!= 0.0),
        vData = v_data,
        vDataMask = .!isnan.(v_data) .& (v_data .!= 0.0)
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
