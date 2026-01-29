# WAVIConstructor.jl - BedMachine v3 data initialization
# Julia port of initBedmachineV3.m

module InitBedMachine

using NCDatasets
using ArchGDAL
using Interpolations
using WAVI

using WAVIConstructor.DataLoading: get_arthern_accumulation, get_zwally_basins, get_albmap, 
    get_bisicles_temps, get_frank_temps, get_measures_velocities, 
    get_smith_dhdt, interpolate_to_grid, get_bedmachine

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
    start_data = get(params, :start_data, "BEDMACHINEV3")
    bedmachine_file = get(params, :bedmachine_file, "Data/BedMachineAntarctica-v3.nc")

    if uppercase(start_data) == "BEDMACHINEV3"
        bed, x, y, geoid, mask, h, s = get_bedmachine(bedmachine_file)

        # Flip arrays for correct orientation
        bed = reverse(bed, dims=1)
        geoid = reverse(geoid, dims=1)
        mask = reverse(mask, dims=1)
        s = reverse(s, dims=1)
        h = reverse(h, dims=1)
    end

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
    
    # Load dhdt data if smith_dhdt_dir is provided
    smith_dir = get(params, :smith_dhdt_dir, "Data/Smith_2020_dhdt")
    
    if smith_dir != "" && isdir(smith_dir)
        smith_data = get_smith_dhdt(smith_dir=smith_dir)
        
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
            # Return zeros if files not found
            Gh = merge(Gh, (dhdt = zeros(size(Gh.surfType)),))
        end
    else
        # Return zeros if directory not provided or doesn't exist
        Gh = merge(Gh, (dhdt = zeros(size(Gh.surfType)),))
    end
    
    # Add aliases for compatibility with SelectDomainWAVI
    Gh = merge(Gh, (
        basinID = Gh.basin_id,
        a = Gh.a_Arthern,
        rock = Gh.rockmask
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

Load accumulation, basin, and ALBMAP data.
"""
function load_additional_datasets(Gh, params)
    # Load Arthern accumulation if file exists
    arthern_file = get(params, :arthern_file, "Data/amsr_accumulation_map.txt")
    if isfile(arthern_file)
        aa_lat, aa_lon, aa_x, aa_y, aa_acc, aa_err = get_arthern_accumulation(arthern_file)
        Gh = merge(Gh, (a_Arthern = interpolate_to_grid(aa_x, aa_y, aa_acc, Gh.xx, Gh.yy),))
    else
        @warn "Arthern accumulation file not found: $arthern_file. Using zeros."
        Gh = merge(Gh, (a_Arthern = zeros(size(Gh.xx)),))
    end

    # Load Zwally basins if file exists
    zwally_file = get(params, :zwally_file, "Data/ZwallyBasins.mat")
    if isfile(zwally_file)
        xx_zwally, yy_zwally, zwally_basins = get_zwally_basins(zwally_file)
        Gh = merge(Gh, (basin_id = interpolate_to_grid(xx_zwally, yy_zwally, zwally_basins, Gh.xx, Gh.yy),))
    else
        @warn "Zwally basins file not found: $zwally_file. Using zeros."
        Gh = merge(Gh, (basin_id = zeros(size(Gh.xx)),))
    end

    # Load ALBMAP (required)
    albmap_file = get(params, :albmap_file, "Data/ALBMAPv1.nc")
    if !isfile(albmap_file)
        error("ALBMAP file not found: $albmap_file. Please provide the ALBMAP file or set params[:albmap_file] to the correct path.")
    end
    albmap_data = get_albmap(albmap_file)
    Gh = merge(Gh, (Tma = interpolate_to_grid(albmap_data[:xx][:], albmap_data[:yy][:], albmap_data[:Tma][:], Gh.xx, Gh.yy),))

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

Load MEaSUREs ice velocity data.
"""
function load_velocity_data(Gu, Gv, Gh, params)
    sub_samp = get(params, :sub_samp, 8)

    velocity_loaded = false
    
    # Try measures_velocity_file (if provided and exists)
    measures_velocity_file = get(params, :measures_velocity_file, "Data/antarctica_ice_velocity_2016_2017_1km_v01.nc")
    if measures_velocity_file != "" && isfile(measures_velocity_file)
        xx_v, yy_v, vx, vy = get_measures_velocities(measures_velocity_file)
        velocity_loaded = true
    elseif measures_velocity_file != ""
        @warn "MEaSUREs velocity file not found: $measures_velocity_file"
    end
    
    if !velocity_loaded
        @warn "Using zero velocities"
        xx_v = Gh.xx
        yy_v = Gh.yy
        vx = zeros(size(Gh.xx))
        vy = zeros(size(Gh.yy))
    end

    # Apply subsampling if needed
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

Load 3D temperature data from various sources.
"""
function load_temperature_data(Gh, params)
    temp_loaded = false

    # Try frank_temps_file first (if provided and exists)
    frank_file = get(params, :frank_temps_file, "Data/FranksTemps.mat")
    if frank_file != "" && isfile(frank_file)
        temp_data = get_frank_temps(frank_file)
        temperature = zeros(size(temp_data.FranksTemps, 1), Gh.nx, Gh.ny)
        for i in 1:size(temp_data.FranksTemps, 1)
            # FranksTemps is (sigma_levels, ny, nx), so we need [i, :, :] to get the 2D slice
            temperature[i, :, :] = interpolate_to_grid(
                temp_data.xxTemp[:], temp_data.yyTemp[:], temp_data.FranksTemps[i, :, :][:],
                Gh.xx, Gh.yy
            )
        end
        sigmas = temp_data.sigmaTemp
        temp_loaded = true
    # Else try bisicles_temps_file (if provided and exists)
    else
        bisicles_file = get(params, :bisicles_temps_file, "Data/antarctica-bisicles-xyzT-8km.nc")
        if bisicles_file != "" && isfile(bisicles_file)
            bisicles_sigma, bisicles_x, bisicles_y, bisicles_z, bisicles_temps = get_bisicles_temps(bisicles_file)

            bisicles_yy = repeat(bisicles_y, 1, length(bisicles_x))
            bisicles_xx = repeat(bisicles_x', length(bisicles_y), 1)

            temperature = zeros(length(bisicles_sigma), Gh.nx, Gh.ny)
            bisicles_mask_full = zeros(length(bisicles_sigma), size(bisicles_xx, 1), size(bisicles_xx, 2))

            for i in 1:length(bisicles_sigma)
                bisicles_mask = bisicles_temps[i, :, :] .> 273.1480
                bisicles_mask_full[i, :, :] = bisicles_mask
                bisicles_temps[i, bisicles_mask_full[i, :, :] .== 1] .= NaN

                this_temp = bisicles_temps[i, :, :]
                valid_mask = .~isnan.(this_temp)
                xx_valid = bisicles_xx[valid_mask]
                yy_valid = bisicles_yy[valid_mask]
                temp_valid = this_temp[valid_mask]

                temperature[i, :, :] = interpolate_to_grid(xx_valid, yy_valid, temp_valid, Gh.xx, Gh.yy)
            end
            sigmas = bisicles_sigma
            temp_loaded = true
        else
            @warn "BISICLES temperature file not found: $bisicles_file"
        end
    end
    
    # Temperature data is required
    if !temp_loaded
        error("Temperature data is required but not found. Please provide one of:\n" *
              "  - Frank temperature file: $(get(params, :frank_temps_file, "Data/FranksTemps.mat"))\n" *
              "  - BISICLES temperature file: $(get(params, :bisicles_temps_file, "Data/antarctica-bisicles-xyzT-8km.nc"))\n" *
              "Ensure the file exists and is accessible.")
    end

    Gh = merge(Gh, (levels = (temperature = temperature, sigmas = sigmas),))
    return Gh
end

end # module InitBedMachine
