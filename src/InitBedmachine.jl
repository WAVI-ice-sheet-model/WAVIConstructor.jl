# WAVIConstructor.jl - BedMachine v3 data initialization
# Julia port of initBedmachineV3.m

module InitBedMachine

using NCDatasets
using ArchGDAL
using MAT
using Interpolations
using WAVI

using WAVIConstructor.DataLoading

"""
    init_bedmachine_v3(params)

Initialize computational grids and load all geophysical data using BedMachine v3.
Julia port of the MATLAB initBedmachineV3 function.

# Arguments
- `params`: NamedTuple with data loading parameters

# Returns
Tuple of grid structures (Gh, Gu, Gv, Gc) with loaded data
"""
function init_bedmachine_v3(params)
    start_data = get(params, :start_data, "BEDMACHINEV3")

    # Load BedMachine data
    if uppercase(start_data) == "BEDMACHINEV3"
        bed, x, y, geoid, mask, h, s = get_bedmachine("Data/BedMachineAntarctica-v3.nc")

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

    # BedMachine v3 grid parameters
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

    isub = i_pole + 2*sub_samp * (-floor(Int, (domain_half_width*2)/(2*sub_samp)):floor(Int, (domain_half_width*2)/(2*sub_samp)))
    jsub = j_pole + 2*sub_samp * (-floor(Int, (domain_half_height*2)/(2*sub_samp)):floor(Int, (domain_half_height*2)/(2*sub_samp)))

    # Subsample all arrays
    bed = bed[isub, jsub]
    h = h[isub, jsub]
    s = s[isub, jsub]
    rockmask = rockmask[isub, jsub]
    geoid = geoid[isub, jsub]
    mask = mask[isub, jsub]

    x = bedmachine_x[isub]
    y = bedmachine_y[jsub]

    # Create H-grid structure
    Gh = create_h_grid(x, y, bed, h, s, geoid, rockmask, mask, params)

    # Load additional datasets
    Gh = load_additional_datasets(Gh, params)

    # Create U/V/C grids
    Gu = create_u_grid(Gh)
    Gv = create_v_grid(Gh)
    Gc = create_c_grid(Gh)

    # Load velocity data
    Gu, Gv = load_velocity_data(Gu, Gv, Gh, params)

    # Load temperature data
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
    # Load Arthern accumulation
    aa_lat, aa_lon, aa_x, aa_y, aa_acc, aa_err = load_arthern_accumulation()
    Gh = merge(Gh, (a_Arthern = interpolate_to_grid(aa_x, aa_y, aa_acc, Gh.xx, Gh.yy)))

    # Load Zwally drainage basins
    xx_zwally, yy_zwally, zwally_basins = load_zwally_basins()
    Gh = merge(Gh, (basin_id = interpolate_to_grid(xx_zwally, yy_zwally, zwally_basins, Gh.xx, Gh.yy)))

    # Load ALBMAP data
    albmap_data = get_albmap("Data/ALBMAPv1.nc")
    Gh = merge(Gh, (Tma = interpolate_to_grid(albmap_data[:xx][:], albmap_data[:yy][:], albmap_data[:Tma][:], Gh.xx, Gh.yy)))

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
    veloc_data = get(params, :veloc_data, "Measures_2014/15")
    sub_samp = get(params, :sub_samp, 8)

    if veloc_data == "Measures_1"
        # Load from MAT file
        velocity_data = load_measures_mat()
        xx_v, yy_v, vx, vy = velocity_data.xx_v, velocity_data.yy_v, velocity_data.vx, velocity_data.vy
    elseif veloc_data in ["Measures_2016/17", "Measures_2014/15", "Measures_phase_v1"]
        # Load from NetCDF
        filename = get_measures_filename(veloc_data)
        xx_v, yy_v, vx, vy = get_measures_velocities("Data/$filename")
    else
        error("Velocity dataset '$veloc_data' not recognized")
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

    # Interpolate to grids
    Gu = merge(Gu, (
        u_data = interpolate_to_grid(xx_v[:], yy_v[:], vx[:], Gu.xx, Gu.yy),
        u_data_mask = .~isnan.(interpolate_to_grid(xx_v[:], yy_v[:], vx[:], Gu.xx, Gu.yy))
    ))

    Gv = merge(Gv, (
        v_data = interpolate_to_grid(xx_v[:], yy_v[:], vy[:], Gv.xx, Gv.yy),
        v_data_mask = .~isnan.(interpolate_to_grid(xx_v[:], yy_v[:], vy[:], Gv.xx, Gv.yy))
    ))

    return Gu, Gv
end

"""
    load_temperature_data(Gh, params)

Load 3D temperature data from various sources.
"""
function load_temperature_data(Gh, params)
    temps = get(params, :temps, "BISICLES_8km")

    if temps == "Frank"
        temp_data = load_frank_temps()
        temperature = zeros(size(temp_data.FranksTemps, 1), Gh.nx, Gh.ny)
        for i in 1:size(temp_data.FranksTemps, 1)
            temperature[i, :, :] = interpolate_to_grid(
                temp_data.xxTemp[:], temp_data.yyTemp[:], temp_data.FranksTemps[i, :][:],
                Gh.xx, Gh.yy
            )
        end
        sigmas = temp_data.sigmaTemp

    elseif temps in ["BISICLES_8km", "BISICLES_1km"]
        filename = temps == "BISICLES_8km" ?
            "Data/antarctica-bisicles-xyzT-8km.nc" :
            "Data/ase-bisicles-xyzT-1km_corrected.nc"

        bisicles_sigma, bisicles_x, bisicles_y, bisicles_z, bisicles_temps = get_bisicles_temps(filename)

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
    else
        error("Temperature dataset '$temps' not recognized")
    end

    Gh = merge(Gh, (levels = (temperature = temperature, sigmas = sigmas)))
    return Gh
end

"""
    interpolate_to_grid(x, y, values, xi, yi)

Interpolate scattered data to a regular grid.
"""
function interpolate_to_grid(x, y, values, xi, yi)
    # Use nearest neighbor interpolation for scattered data
    itp = interpolate((x, y), values, Gridded(Linear()))
    return itp(xi, yi)
end

function get_measures_filename(veloc_data)
    filename_map = Dict(
        "Measures_2016/17" => "antarctica_ice_velocity_2016_2017_1km_v01.nc",
        "Measures_2014/15" => "antarctica_ice_velocity_2014_2015_1km_v01.nc",
        "Measures_phase_v1" => "antarctic_ice_vel_phase_map_v01.nc"
    )
    return get(filename_map, veloc_data, "")
end

end # module InitBedMachine
