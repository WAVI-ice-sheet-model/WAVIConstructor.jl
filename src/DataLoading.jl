module DataLoading

using ArchGDAL
using NCDatasets

export get_albmap, get_bedmachine, get_bisicles_temps, get_measures_velocities, geotiff_read_axis_only, get_smith_dhdt


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
    bed = ds["bed"][:]
    x = ds["x"][:]
    y = ds["y"][:]
    geoid = ds["geoid"][:]
    mask = ds["mask"][:]
    surface = ds["surface"][:]
    thickness = ds["thickness"][:]
    firn = ds["firn"][:]
    close(ds)
    return bed, x, y, geoid, mask, thickness, surface, firn
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

# Throws
- `ArgumentError`: If the specified file does not exist

# Example
```julia
sigma, x, y, z, temps = get_bisicles_temps("antarctica-bisicles-xyzT-8km.nc")
sigma, x, y, z, temps = get_bisicles_temps("antarctica-bisicles-xyzT-8km.nc", scale_xy=1000)  # Convert km to m
sigma, x, y, z, temps = get_bisicles_temps("antarctica-bisicles-xyzT-8km.nc", scale_xy=0.001)  # Convert m to km
```
"""
function get_bisicles_temps(fname_start::String; scale_xy::Real=1)
    # Validate input filename
    if !isfile(fname_start)
        throw(ArgumentError("File does not exist: $fname_start"))
    end

    # Use try-finally to ensure file is always closed
    ds = NCDataset(fname_start)
    try
        # Helper function to read and optionally scale data
        read_data(varname, scale_factor=1) = scale_factor .* ds[varname][:]

        # Read data from NetCDF file
        bisicles_sigma = read_data("sigma")
        bisicles_temps = read_data("T")
        bisicles_x = read_data("x", scale_xy)
        bisicles_y = read_data("y", scale_xy)
        bisicles_z = read_data("z")

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
        # Read coordinate and velocity data
        Measures_x = ds["x"][:]
        Measures_y = ds["y"][:]
        VX = ds["VX"][:]
        VY = ds["VY"][:]
        
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
    # Validate file exists
    if !isfile(filename)
        throw(ArgumentError("File does not exist: $filename"))
    end

    ArchGDAL.read(filename) do dataset
        # Get basic image info
        width = ArchGDAL.width(dataset)
        height = ArchGDAL.height(dataset)
        bands = ArchGDAL.nband(dataset)
        
        # Get geotransform information using GDAL standard
        gt = ArchGDAL.getgeotransform(dataset)
        
        # Create a logical interface
        geotransform = (
            x_origin = gt[1],
            pixel_width = gt[2], 
            y_origin = gt[4],
            pixel_height = gt[6]
        )
        
        # Extract the specific values we need
        dx = geotransform.pixel_width    # pixel width (x-direction resolution)
        dy = -geotransform.pixel_height  # pixel height (make positive, GDAL uses negative for north-up)
        mapx = geotransform.x_origin     # x-coordinate of upper-left corner
        mapy = geotransform.y_origin     # y-coordinate of upper-left corner  
        
        # Calculate coordinate arrays
        x_full = [mapx + (i-1) * dx for i in 1:width]
        y_full = [mapy - (i-1) * dy for i in 1:height]
        
        # Handle subsetting
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
    get_smith_dhdt(filename::String)

Load dhdt (elevation/thickness change) data from a single Smith et al. 2020 GeoTIFF file.

# Arguments
- `filename::String`: Path to the GeoTIFF file

# Returns
A tuple containing:
- `xx`: 2D grid of x coordinates
- `yy`: 2D grid of y coordinates  
- `dhdt`: Elevation/thickness change data

# Throws
- `ArgumentError`: If the specified file does not exist

# Example
```julia
xx, yy, dhdt = get_smith_dhdt("ais_grounded.tif")
```
"""
function get_smith_dhdt(filename::String)
    # Validate file exists
    if !isfile(filename)
        throw(ArgumentError("File does not exist: $filename"))
    end

    # Read dhdt data
    dhdt = ArchGDAL.read(filename) do dataset
        ArchGDAL.read(dataset, 1)  # Read band 1
    end

    # Get coordinate information
    coord_info = geotiff_read_axis_only(filename)

    # Calculate pixel center coordinates (PIXEL IS AREA, coordinate is North West corner)
    x = [coord_info.info.map_info.mapx + coord_info.info.map_info.dx * (i - 0.5) for i in 1:coord_info.info.samples]
    y = [coord_info.info.map_info.mapy - coord_info.info.map_info.dy * (i - 0.5) for i in 1:coord_info.info.lines]

    # Create coordinate grids equivalent to MATLAB's ndgrid(flip(y), x)
    yy = repeat(reverse(y), 1, length(x))
    xx = repeat(x', length(y), 1)

    return xx, yy, dhdt
end

end # module