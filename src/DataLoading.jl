module DataLoading

using ArchGDAL
using NCDatasets
using MAT
using NearestNeighbors
import DelimitedFiles: readdlm

export get_albmap, get_bedmachine, get_bisicles_temps, get_measures_velocities, geotiff_read_axis_only, get_smith_dhdt, get_arthern_accumulation, get_zwally_basins, get_frank_temps, get_measures_mat


"""
    get_albmap(filename::String)

Load and process the ALBMAP NetCDF file, returning a dictionary of derived fields.

# Arguments
- `filename::String`: Path to the ALBMAP NetCDF file.

# Returns
- `Dict`: Dictionary with grid, geometry, and accumulation fields.
    - `:dx`, `:dy`: Grid spacing (x, y)
    - `:nx`, `:ny`: Number of grid points (x, y)
    - `:x0`, `:y0`: Grid origin (x, y)
    - `:xx`, `:yy`: 2D coordinate grids
    - `:h`: Ice thickness
    - `:s`: Surface elevation above firn
    - `:b`: Bed elevation
    - `:firn`: Firn thickness
    - `:Tma`: Mean annual temperature (K)
    - `:a_Arthern`, `:a_vdBerg`: Accumulation fields
    - `:surfType`: Surface type mask
"""
function get_albmap(filename::String)
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

    return Dict(
        :dx => dx,
        :dy => dy,
        :nx => nx,
        :ny => ny,
        :x0 => x0,
        :y0 => y0,
        :xx => xx,
        :yy => yy,
        :h => h,
        :s => s,
        :b => b,
        :firn => firn_out,
        :Tma => Tma,
        :a_Arthern => a_Arthern,
        :a_vdBerg => a_vdBerg,
        :surfType => surfType
    )
end


"""
    get_bedmachine(filename::String)

Load BedMachine NetCDF variables and return them.

# Arguments
- `filename::String`: Path to the BedMachine NetCDF file.

# Returns
- Tuple: (bed, x, y, geoid, mask, thickness, surface, firn)
    - `bed`: Bed elevation
    - `x`: x-coordinates
    - `y`: y-coordinates
    - `geoid`: Geoid height
    - `mask`: Mask array
    - `thickness`: Ice thickness
    - `surface`: Surface elevation
    - `firn`: Firn thickness
"""
function get_bedmachine(filename::String)
    ds = NCDataset(filename)
    try
        bed = Array(ds["bed"])
        x = Array(ds["x"])
        y = Array(ds["y"])
        geoid = Array(ds["geoid"])
        mask = Array(ds["mask"])
        surface = Array(ds["surface"])
        thickness = Array(ds["thickness"])
        firn = Array(ds["firn"])
        return bed, x, y, geoid, mask, thickness, surface, firn
    finally
        close(ds)
    end
end


"""
    get_bisicles_temps(fname_start::String; scale_xy::Real=1)

Load temperature and coordinate data from a BISICLES NetCDF file.

# Arguments
- `fname_start::String`: Path to the NetCDF file containing BISICLES temperature data
- `scale_xy::Real=1`: Scaling factor to apply to x and y coordinates (e.g., 1000 to convert km to m)

# Returns
A tuple containing:
- `bisicles_sigma`: Sigma coordinate values
- `bisicles_x`: X coordinate values (scaled by `scale_xy`)
- `bisicles_y`: Y coordinate values (scaled by `scale_xy`)
- `bisicles_z`: Z coordinate values
- `bisicles_temps`: Temperature data

# Example
```julia
sigma, x, y, z, temps = get_bisicles_temps("antarctica-bisicles-xyzT-8km.nc")
sigma, x, y, z, temps = get_bisicles_temps("antarctica-bisicles-xyzT-8km.nc", scale_xy=1000)  # Convert km to m
sigma, x, y, z, temps = get_bisicles_temps("antarctica-bisicles-xyzT-8km.nc", scale_xy=0.001)  # Convert m to km
```
"""
function get_bisicles_temps(fname_start::String; scale_xy::Real=1)
    ds = NCDataset(fname_start)
    try
        bisicles_sigma = vec(Array(ds["sigma"]))  # Ensure it's a vector
        bisicles_temps = Array(ds["T"])
        bisicles_x = vec(scale_xy .* Array(ds["x"]))  # Ensure it's a vector
        bisicles_y = vec(scale_xy .* Array(ds["y"]))  # Ensure it's a vector
        bisicles_z = Array(ds["z"])

        return bisicles_sigma, bisicles_x, bisicles_y, bisicles_z, bisicles_temps
    finally
        close(ds)
    end
end

"""
    get_measures_velocities(filename::String)

Load velocity data from a NetCDF file.

# Arguments
- `filename::String`: Path to the NetCDF file containing velocity data

# Returns
A tuple containing:
- `xx_v`: 2D grid of x coordinates
- `yy_v`: 2D grid of y coordinates  
- `VX`: X-component velocity data
- `VY`: Y-component velocity data

# Example
```julia
xx, yy, vx, vy = get_measures_velocities("Antarctic_ice_velocity_2016_2017_1km_v01.nc")
```
"""
function get_measures_velocities(filename::String)
    ds = NCDataset(filename)
    try
        Measures_x = Array(ds["x"])
        Measures_y = Array(ds["y"])
        VX = Array(ds["VX"])
        VY = Array(ds["VY"])
        
        # Create coordinate grids (equivalent to MATLAB's ndgrid)
        xx_v = repeat(Measures_x, 1, length(Measures_y))
        yy_v = repeat(Measures_y', length(Measures_x), 1)

        return xx_v, yy_v, VX, VY
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
    get_smith_dhdt(; grnd_file=nothing, flt_file=nothing, smith_dir="Data/Smith_2020_dhdt")

Load dhdt (elevation/thickness change) data from Smith et al. 2020 GeoTIFF file(s).

# Arguments
- `grnd_file::Union{String,Nothing}`: Path to grounded ice GeoTIFF file (optional)
- `flt_file::Union{String,Nothing}`: Path to floating ice GeoTIFF file (optional)
- `smith_dir::String`: Directory containing default files if grnd_file/flt_file not specified

# Returns
If both grnd_file and flt_file are provided (or found in smith_dir):
- NamedTuple with (grnd_xx, grnd_yy, grnd_dhdt, flt_xx, flt_yy, flt_dhdt)

If only one file is provided:
- Tuple (xx, yy, dhdt)

# Examples
```julia
# Load single file
xx, yy, dhdt = get_smith_dhdt(grnd_file="ais_grounded.tif")

# Load both files from directory
data = get_smith_dhdt(smith_dir="Data/Smith_2020_dhdt")

# Load both files explicitly
data = get_smith_dhdt(grnd_file="ais_grounded.tif", flt_file="ais_floating.tif")
```
"""
function get_smith_dhdt(; grnd_file=nothing, flt_file=nothing, smith_dir="Data/Smith_2020_dhdt")
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
        # Load only grounded file
        if !isfile(grnd_file)
            return nothing
        end
        return load_single_file(grnd_file)
    elseif flt_file !== nothing
        # Load only floating file
        if !isfile(flt_file)
            return nothing
        end
        return load_single_file(flt_file)
    end
end

"""
    get_arthern_accumulation(filename::String="Data/amsr_accumulation_map.txt")

Load Arthern accumulation data from a text file.

# Arguments
- `filename::String`: Path to the Arthern accumulation text file (default: "Data/amsr_accumulation_map.txt")

# Returns
A tuple containing:
- `aa_lat`: Latitude values
- `aa_lon`: Longitude values
- `aa_x`: X-coordinate values (in meters)
- `aa_y`: Y-coordinate values (in meters)
- `aa_acc`: Accumulation values (converted to ice equivalent, divided by 917)
- `aa_err`: Error values

# Notes
- The file should have 21 header lines
- Data is space-delimited with 6 columns
- Accumulation values are converted from water equivalent to ice equivalent (divided by 917 kg/m³)
- Rows with NaN accumulation values are filtered out
"""
function get_arthern_accumulation(filename::String="Data/amsr_accumulation_map.txt")
    # Read the file, skipping 21 header lines
    data = readdlm(filename, ' ', Float64; skipstart=21)

    aa_lat = data[:, 1]
    aa_lon = data[:, 2]
    aa_x = data[:, 3]
    aa_y = data[:, 4]
    aa_acc_raw = data[:, 5]
    aa_err = data[:, 6]

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

    return aa_lat, aa_lon, aa_x, aa_y, aa_acc, aa_err
end

"""
    get_zwally_basins(filename::String="Data/DrainageBasins/ZwallyBasins.mat")

Load Zwally drainage basins from a MATLAB .mat file.

# Arguments
- `filename::String`: Path to the Zwally basins .mat file (default: "Data/DrainageBasins/ZwallyBasins.mat")

# Returns
A tuple containing:
- `xx_zwally`: X-coordinate values (filtered to only include points where ZwallyBasins > 0)
- `yy_zwally`: Y-coordinate values (filtered to only include points where ZwallyBasins > 0)
- `zwally_basins`: Basin ID values (filtered to only include points where ZwallyBasins > 0)

# Notes
- The .mat file should contain variables: `xxZwallyBasins`, `yyZwallyBasins`, `ZwallyBasins`
- Only points where `ZwallyBasins > 0` are returned
"""
function get_zwally_basins(filename::String="Data/DrainageBasins/ZwallyBasins.mat")
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

        return xx_zwally, yy_zwally, zwally_basins
    end
end

"""
    get_frank_temps(filename::String="Data/FranksTemps.mat")

Load Frank temperature data from a MATLAB .mat file.

# Arguments
- `filename::String`: Path to the FranksTemps .mat file (default: "Data/FranksTemps.mat")

# Returns
A NamedTuple containing:
- `FranksTemps`: 3D temperature data array (levels × points)
- `xxTemp`: X-coordinate values
- `yyTemp`: Y-coordinate values
- `sigmaTemp`: Sigma coordinate values (vertical levels)

# Notes
- The .mat file should contain variables: `FranksTemps`, `xxTemp`, `yyTemp`, `sigmaTemp`
- FranksTemps is a 2D array where rows represent vertical levels and columns represent spatial points
"""
function get_frank_temps(filename::String="Data/FranksTemps.mat")
    matopen(filename) do mat_file
        FranksTemps = read(mat_file, "FranksTemps")
        xxTemp = read(mat_file, "xxTemp")
        yyTemp = read(mat_file, "yyTemp")
        sigmaTemp = read(mat_file, "sigmaTemp")

        return (
            FranksTemps = FranksTemps,
            xxTemp = xxTemp,
            yyTemp = yyTemp,
            sigmaTemp = vec(sigmaTemp)  # Ensure it's a vector, not a matrix
        )
    end
end

"""
    get_measures_mat(filename::String="Data/MEaSUREs/MEaSUREsAntVels.mat")

Load MEaSUREs ice velocity data from a MATLAB .mat file.

# Arguments
- `filename::String`: Path to the MEaSUREs velocity .mat file (default: "Data/MEaSUREs/MEaSUREsAntVels.mat")

# Returns
A NamedTuple containing:
- `xx_v`: 2D grid of x coordinates
- `yy_v`: 2D grid of y coordinates
- `vx`: X-component velocity data
- `vy`: Y-component velocity data

# Notes
- The .mat file should contain variables: `xx_v`, `yy_v`, `vx`, `vy`
- This is the older MEaSUREs dataset format (Measures_1)
"""
function get_measures_mat(filename::String="Data/MEaSUREs/MEaSUREsAntVels.mat")
    matopen(filename) do mat_file
        xx_v = read(mat_file, "xx_v")
        yy_v = read(mat_file, "yy_v")
        vx = read(mat_file, "vx")
        vy = read(mat_file, "vy")

        return (
            xx_v = xx_v,
            yy_v = yy_v,
            vx = vx,
            vy = vy
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

end # module