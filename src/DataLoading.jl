module DataLoading

using NCDatasets

export get_albmap, get_bedmachine


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

end # module
