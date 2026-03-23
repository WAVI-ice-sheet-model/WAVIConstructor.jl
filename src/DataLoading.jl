module DataLoading

using ArchGDAL
using NCDatasets
using MAT
using NearestNeighbors
import DelimitedFiles: readdlm

using WAVIConstructor.DataSources

export load_data, interpolate_to_grid, interpolate_temperature, geotiff_read_axis_only


"""
    load_data(::ALBMAPv1, filename::String)

Load and process the ALBMAP NetCDF file.

# Returns
A `NamedTuple` with grid, surface, and temperature fields:
  `(dx, dy, nx, ny, x0, y0, xx, yy, h, s, b, firn, Tma, a_Arthern, a_vdBerg, surfType)`
"""
function load_data(::ALBMAPv1, filename::String)
    data = Dict()
    ds = NCDataset(filename)
    for var in keys(ds)
        data[var] = ds[var][:]
    end
    close(ds)

    x1, y1 = data["x1"], data["y1"]
    dx, dy = float(x1[2] - x1[1]), float(y1[2] - y1[1])
    nx, ny = length(x1), length(y1)
    x0, y0 = float(x1[1]) - 0.5 * dx, float(y1[1]) - 0.5 * dy
    xx = [x0 + (i-0.5)*dx for i in 1:nx, j in 1:ny]
    yy = [y0 + (j-0.5)*dy for i in 1:nx, j in 1:ny]

    usrf, lsrf, firn = data["usrf"], data["lsrf"], data["firn"]
    h = permutedims(usrf - lsrf - firn)
    s = permutedims(usrf - firn)
    b = permutedims(data["topg"])
    firn_out = permutedims(data["firn"])
    Tma = 273.16 .+ permutedims(data["temp"])
    a_Arthern = permutedims(data["acca"])
    a_vdBerg = permutedims(data["accr"])
    surfType = permutedims(data["mask_plus"])

    return (
        dx = dx,
        dy = dy,
        nx = nx,
        ny = ny,
        x0 = x0,
        y0 = y0,
        xx = xx,
        yy = yy,
        h = h,
        s = s,
        b = b,
        firn = firn_out,
        Tma = Tma,
        a_Arthern = a_Arthern,
        a_vdBerg = a_vdBerg,
        surfType = surfType,
    )
end


"""
    load_data(::BedMachineV3, filename::String)

Load BedMachine v3 NetCDF variables.

# Returns
A `NamedTuple`: `(bed, x, y, geoid, mask, thickness, surface)`
"""
function load_data(::BedMachineV3, filename::String)
    ds = NCDataset(filename)
    try
        bed = Array(ds["bed"])
        x = Array(ds["x"])
        y = Array(ds["y"])
        geoid = Array(ds["geoid"])
        mask = Array(ds["mask"])
        surface = Array(ds["surface"])
        thickness = Array(ds["thickness"])
        return (bed = bed, x = x, y = y, geoid = geoid, mask = mask,
                thickness = thickness, surface = surface)
    finally
        close(ds)
    end
end


"""
    load_data(::BISICLESTemps, filename::String; scale_xy::Real=1)

Load temperature and coordinate data from a BISICLES NetCDF file.

# Returns
A `NamedTuple`: `(temps, xx, yy, sigma)` — same shape as `load_data(::FrankTemps, …)`.
"""
function load_data(::BISICLESTemps, filename::String; scale_xy::Real=1)
    ds = NCDataset(filename)
    try
        bisicles_sigma = vec(Array(ds["sigma"]))
        bisicles_temps = Array(ds["T"])
        bisicles_x = vec(scale_xy .* Array(ds["x"]))
        bisicles_y = vec(scale_xy .* Array(ds["y"]))  

        return (temps = bisicles_temps, xx = bisicles_x, yy = bisicles_y, sigma = bisicles_sigma)
    finally
        close(ds)
    end
end

"""
    load_data(::MEaSUREs, filename::String)

Load ice-velocity data from a MEaSUREs NetCDF file.

# Returns
A `NamedTuple`: `(xx, yy, vx, vy)`
"""
function load_data(::MEaSUREs, filename::String)
    ds = NCDataset(filename)
    try
        Measures_x = Array(ds["x"])
        Measures_y = Array(ds["y"])
        VX = Array(ds["VX"])
        VY = Array(ds["VY"])
        
        # Create coordinate grids (equivalent to MATLAB's ndgrid)
        xx_v = repeat(Measures_x, 1, length(Measures_y))
        yy_v = repeat(Measures_y', length(Measures_x), 1)

        return (xx = xx_v, yy = yy_v, vx = VX, vy = VY)
    finally
        close(ds)
    end
end


"""
    geotiff_read_axis_only(filename::String; pixel_subset=nothing, map_subset=nothing)

Read GeoTIFF file and extract coordinate information without reading the full raster data.

# Arguments
- `filename::String`: Path to the GeoTIFF file
- `pixel_subset`: Optional array [minx, maxx, miny, maxy] for pixel-based subset
- `map_subset`: Optional array [minx, maxx, miny, maxy] for map coordinate-based subset

# Returns
A NamedTuple containing:
- `x`: X coordinates
- `y`: Y coordinates  
- `info`: Metadata information including map_info with dx, dy, mapx, mapy

# Example
```julia
result = geotiff_read_axis_only("myfile.tif")
x_coords = result.x
y_coords = result.y
```
"""
function geotiff_read_axis_only(filename::String; pixel_subset=nothing, map_subset=nothing)
    ArchGDAL.read(filename) do dataset
        # Get basic image info
        width = ArchGDAL.width(dataset)
        height = ArchGDAL.height(dataset)
        bands = ArchGDAL.nraster(dataset)
        
        # Get geotransform information using GDAL standard
        gt = ArchGDAL.getgeotransform(dataset)
        
        geotransform = (
            x_origin = gt[1],
            pixel_width = gt[2], 
            y_origin = gt[4],
            pixel_height = gt[6]
        )
        
        # Calculate pixel width and height (x-direction resolution)
        dx = geotransform.pixel_width    # pixel width (x-direction resolution)
        dy = -geotransform.pixel_height  # pixel height (make positive, GDAL uses negative for north-up)
        mapx = geotransform.x_origin     # x-coordinate of upper-left corner
        mapy = geotransform.y_origin     # y-coordinate of upper-left corner  
        
        x_full = [mapx + (i-1) * dx for i in 1:width]
        y_full = [mapy - (i-1) * dy for i in 1:height]
        
        sub = [1, width, 1, height]  # default: full image
        
        if pixel_subset !== nothing
            sub = pixel_subset
        elseif map_subset !== nothing
            # Convert map coordinates to pixel coordinates
            subx1 = round(Int, (map_subset[1] - mapx) / dx + 1)
            subx2 = round(Int, (map_subset[2] - mapx) / dx + 1)
            suby1 = round(Int, (mapy - map_subset[4]) / dy + 1)
            suby2 = round(Int, (mapy - map_subset[3]) / dy + 1)
            
            # Clamp to valid ranges
            subx1 = clamp(subx1, 1, width)
            subx2 = clamp(subx2, 1, width)
            suby1 = clamp(suby1, 1, height)
            suby2 = clamp(suby2, 1, height)
            
            sub = [subx1, subx2, suby1, suby2]
        end
        
        # Extract subset coordinates
        x_coords = x_full[sub[1]:sub[2]]
        y_coords = y_full[sub[3]:sub[4]]
        
        # Create map_info structure
        map_info = (
            dx = dx,
            dy = dy,
            mapx = mapx,
            mapy = mapy
        )
        
        # Create subset info (if subsetting was applied)
        sub_info = nothing
        if pixel_subset !== nothing || map_subset !== nothing
            sub_info = (
                samples = sub[2] - sub[1] + 1,
                lines = sub[4] - sub[3] + 1,
                mapx = [x_coords[1], x_coords[end]],
                mapy = [y_coords[1], y_coords[end]], 
                pixx = [sub[1], sub[2]],
                pixy = [sub[3], sub[4]]
            )
        end
        
        # Create complete info structure matching MATLAB output
        info = (
            samples = width,
            lines = height,
            bands = bands,
            map_info = map_info,
            sub = sub_info
        )
        
        return (x = x_coords, y = y_coords, info = info)
    end
end

"""
    load_data(::SmithDhdt, smith_dir::String; grnd_file=nothing, flt_file=nothing)

Load dh/dt (elevation/thickness change) data from Smith et al. 2020 GeoTIFF files.

# Returns
A `NamedTuple`: `(grnd_xx, grnd_yy, grnd_dhdt, flt_xx, flt_yy, flt_dhdt)`,
or `nothing` if the files cannot be found.
"""
function load_data(::SmithDhdt, smith_dir::String; grnd_file=nothing, flt_file=nothing)
    # Helper function to load a single file
    function load_single_file(filename)
        dhdt_raw = ArchGDAL.read(filename) do dataset
            ArchGDAL.read(dataset, 1)  # Read band 1
        end

        coord_info = geotiff_read_axis_only(filename)

        # Calculate pixel centre coordinates (PIXEL IS AREA, coordinate is North West corner)
        x = [coord_info.info.map_info.mapx + coord_info.info.map_info.dx * (i - 0.5) for i in 1:coord_info.info.samples]
        y = [coord_info.info.map_info.mapy - coord_info.info.map_info.dy * (i - 0.5) for i in 1:coord_info.info.lines]

        # Create coordinate grids equivalent to MATLAB's ndgrid(flip(y), x)
        yy = repeat(reverse(y), 1, length(x))
        xx = repeat(x', length(y), 1)

        # Transpose dhdt to match coordinate grid shape (height, width)
        # ArchGDAL reads data as (width, height), but we want (height, width) to match xx and yy
        if size(dhdt_raw) == (length(x), length(y))
            dhdt = permutedims(dhdt_raw, (2, 1))
        else
            dhdt = dhdt_raw
        end
        
        # Flip to match MATLAB behavior
        dhdt_flipped = reverse(dhdt, dims=1)

        return xx, yy, dhdt_flipped
    end
    
    # Determine which files to load
    if grnd_file === nothing && flt_file === nothing
        # Load both from default directory
        grnd_file = joinpath(smith_dir, "ais_grounded.tif")
        flt_file = joinpath(smith_dir, "ais_floating.tif")
    end
    
    # Load files
    if grnd_file !== nothing && flt_file !== nothing
        # Load both files
        if !isfile(grnd_file) || !isfile(flt_file)
            return nothing
        end
        
        grnd_xx, grnd_yy, grnd_dhdt = load_single_file(grnd_file)
        flt_xx, flt_yy, flt_dhdt = load_single_file(flt_file)
        
        return (
            grnd_xx = grnd_xx,
            grnd_yy = grnd_yy,
            grnd_dhdt = grnd_dhdt,
            flt_xx = flt_xx,
            flt_yy = flt_yy,
            flt_dhdt = flt_dhdt
        )
    elseif grnd_file !== nothing
        # Load only grounded file — store in grnd fields, zeros for flt
        if !isfile(grnd_file)
            return nothing
        end
        gxx, gyy, gdhdt = load_single_file(grnd_file)
        return (
            grnd_xx = gxx, grnd_yy = gyy, grnd_dhdt = gdhdt,
            flt_xx = gxx, flt_yy = gyy, flt_dhdt = zeros(size(gdhdt)),
        )
    elseif flt_file !== nothing
        # Load only floating file — store in flt fields, zeros for grnd
        if !isfile(flt_file)
            return nothing
        end
        fxx, fyy, fdhdt = load_single_file(flt_file)
        return (
            grnd_xx = fxx, grnd_yy = fyy, grnd_dhdt = zeros(size(fdhdt)),
            flt_xx = fxx, flt_yy = fyy, flt_dhdt = fdhdt,
        )
    end
end

"""
    load_data(::ArthernAccumulation, filename::String)

Load Arthern accumulation data from a text file.

# Returns
A `NamedTuple`: `(x, y, acc)` — only the fields used downstream.
"""
function load_data(::ArthernAccumulation, filename::String)
    # Read the file content, handling potential Windows line endings
    content = read(filename, String)
    content = replace(content, "\r\n" => "\n")  # Convert Windows to Unix line endings
    content = replace(content, "\r" => "\n")    # Handle old Mac line endings
    
    # Split into lines and skip header
    lines = split(content, '\n')
    data_lines = lines[22:end]  # Skip 21 header lines (1-indexed)
    
    # Parse data - handle multiple spaces by splitting on whitespace
    aa_lat = Float64[]
    aa_lon = Float64[]
    aa_x = Float64[]
    aa_y = Float64[]
    aa_acc_raw = Float64[]
    aa_err = Float64[]
    
    for line in data_lines
        stripped = strip(line)
        isempty(stripped) && continue
        
        # Split on whitespace (handles multiple spaces)
        parts = split(stripped)
        length(parts) < 6 && continue
        
        try
            push!(aa_lat, parse(Float64, parts[1]))
            push!(aa_lon, parse(Float64, parts[2]))
            push!(aa_x, parse(Float64, parts[3]))
            push!(aa_y, parse(Float64, parts[4]))
            push!(aa_acc_raw, parse(Float64, parts[5]))
            push!(aa_err, parse(Float64, parts[6]))
        catch e
            # Skip lines that can't be parsed (e.g., malformed data)
            continue
        end
    end

    # Filter out rows where accumulation is NaN
    valid_mask = .!isnan.(aa_acc_raw)

    # Apply mask to all arrays
    aa_lat = aa_lat[valid_mask]
    aa_lon = aa_lon[valid_mask]
    aa_x = aa_x[valid_mask]
    aa_y = aa_y[valid_mask]
    aa_err = aa_err[valid_mask]

    # Convert accumulation from water equivalent to ice equivalent (divide by 917 kg/m³)
    aa_acc = aa_acc_raw[valid_mask] ./ 917.0

    return (x = aa_x, y = aa_y, acc = aa_acc)
end

"""
    load_data(::ZwallyBasins, filename::String)

Load Zwally drainage basins from a MATLAB .mat file.

# Returns
A `NamedTuple`: `(xx, yy, basins)` — filtered to points where basins > 0.
"""
function load_data(::ZwallyBasins, filename::String)
    matopen(filename) do mat_file
        xx_zwally_full = read(mat_file, "xxZwallyBasins")
        yy_zwally_full = read(mat_file, "yyZwallyBasins")
        zwally_basins_full = read(mat_file, "ZwallyBasins")

        # Filter to only include points where ZwallyBasins > 0
        valid_mask = zwally_basins_full .> 0

        # Apply mask - works for both 1D and 2D arrays
        xx_zwally = xx_zwally_full[valid_mask]
        yy_zwally = yy_zwally_full[valid_mask]
        zwally_basins = zwally_basins_full[valid_mask]

        return (xx = xx_zwally, yy = yy_zwally, basins = zwally_basins)
    end
end

"""
    load_data(::FrankTemps, filename::String)

Load Frank temperature data from a MATLAB .mat file.

# Returns
A `NamedTuple`: `(temps, xx, yy, sigma)` — same shape as `load_data(::BISICLESTemps, …)`.
"""
function load_data(::FrankTemps, filename::String)
    matopen(filename) do mat_file
        FranksTemps = read(mat_file, "FranksTemps")
        xxTemp = read(mat_file, "xxTemp")
        yyTemp = read(mat_file, "yyTemp")
        sigmaTemp = read(mat_file, "sigmaTemp")

        return (
            temps = FranksTemps,
            xx = xxTemp,
            yy = yyTemp,
            sigma = vec(sigmaTemp),
        )
    end
end

"""
    interpolate_to_grid(x, y, values, xi, yi)

Interpolate scattered data to a regular grid using nearest neighbor interpolation.
This is equivalent to MATLAB's TriScatteredInterp with 'nearest' method.

# Arguments
- `x, y`: 1D arrays of scattered point coordinates
- `values`: 1D array of values at scattered points
- `xi, yi`: 2D arrays of target grid coordinates

# Returns
- 2D array of interpolated values on the target grid
"""
function interpolate_to_grid(x, y, values, xi, yi)
    # Flatten target grid coordinates
    xi_flat = vec(xi)
    yi_flat = vec(yi)
    
    # Build KDTree for efficient nearest neighbor search
    points = hcat(x, y)
    tree = KDTree(points')
    
    # Find nearest neighbors for each target point
    indices, _ = knn(tree, hcat(xi_flat, yi_flat)', 1, true)
    
    # Extract values at nearest neighbors
    result_flat = [values[idx[1]] for idx in indices]
    
    # Reshape to match target grid
    return reshape(result_flat, size(xi))
end

# ── Temperature interpolation (dispatch per source) ──────────────────

"""
    interpolate_temperature(source::TemperatureSource, temp_data, Gh)

Interpolate raw temperature data onto the H-grid.  Each concrete
`TemperatureSource` subtype implements its own method so that
source-specific pre-processing (masking, grid construction, etc.)
is handled through dispatch rather than if/else chains.

# Returns
`(temperature, sigmas)` where `temperature` is `(nσ, nx, ny)` and
`sigmas` is a 1-D vector of sigma levels.
"""
function interpolate_temperature end   # forward declaration for docstring

"""
    interpolate_temperature(::FrankTemps, temp_data, Gh)

FranksTemps stores `(sigma_levels, ny, nx)`.  `xx` / `yy` are 1-D
vectors that can be passed straight to `interpolate_to_grid`.
"""
function interpolate_temperature(::FrankTemps, temp_data, Gh)
    temps_raw = temp_data.temps
    sigmas    = temp_data.sigma

    temperature = zeros(size(temps_raw, 1), Gh.nx, Gh.ny)
    for i in 1:size(temps_raw, 1)
        temperature[i, :, :] = interpolate_to_grid(
            temp_data.xx[:], temp_data.yy[:], temps_raw[i, :, :][:],
            Gh.xx, Gh.yy
        )
    end
    return temperature, sigmas
end

"""
    interpolate_temperature(::BISICLESTemps, temp_data, Gh)

BISICLES stores `(sigma_levels, ny, nx)` with 1-D coordinate vectors.
Values above 273.148 K are treated as invalid ocean fill and masked
before interpolation.
"""
function interpolate_temperature(::BISICLESTemps, temp_data, Gh)
    temps_raw = copy(temp_data.temps)       # copy so we can NaN-mask in place
    sigmas    = temp_data.sigma

    bisicles_yy = repeat(temp_data.yy, 1, length(temp_data.xx))
    bisicles_xx = repeat(temp_data.xx', length(temp_data.yy), 1)

    temperature = zeros(length(sigmas), Gh.nx, Gh.ny)

    for i in 1:length(sigmas)
        ocean_mask = temps_raw[i, :, :] .> 273.1480
        temps_raw[i, ocean_mask] .= NaN

        this_temp  = temps_raw[i, :, :]
        valid_mask = .!isnan.(this_temp)
        xx_valid   = bisicles_xx[valid_mask]
        yy_valid   = bisicles_yy[valid_mask]
        temp_valid = this_temp[valid_mask]

        temperature[i, :, :] = interpolate_to_grid(xx_valid, yy_valid, temp_valid, Gh.xx, Gh.yy)
    end
    return temperature, sigmas
end

# ── NoData fallback methods ───────────────────────────────────────────

"""
    load_data(::NoData, ::String; kwargs...)

Generic no-op: returns `nothing` when a category is disabled.
"""
load_data(::NoData, ::String; kwargs...) = nothing

end # module