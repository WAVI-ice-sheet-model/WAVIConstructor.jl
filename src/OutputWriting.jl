# WAVIConstructor.jl - Output Writing
# Flexible output writing supporting binary (.bin) and NetCDF (.nc) formats.

module OutputWriting

using NCDatasets
using Dates

export write_output, write_binary_files, write_netcdf_file

# ── Variable metadata ─────────────────────────────────────────────────
# Each entry: (filename / varname, long_name, units, standard_name or nothing)

const _H_GRID_VARS = [
    ("thickness",          "Ice thickness",                       "m",      "land_ice_thickness"),
    ("surface",            "Ice surface elevation",               "m",      "surface_altitude"),
    ("bed",                "Bed topography",                      "m",      "bedrock_altitude"),
    ("h_mask",             "Ice grid mask",                       "1",      nothing),
    ("basinID",            "Drainage basin ID",                   "1",      nothing),
    ("accumulation_data",  "Surface mass balance (accumulation)", "m yr-1", "land_ice_surface_specific_mass_balance"),
    ("dhdt_data",          "Ice thickness change rate",           "m yr-1", "tendency_of_land_ice_thickness"),
    ("dhdt_acc_mask",      "dh/dt and accumulation data mask",    "1",      nothing),
]

const _U_GRID_VARS = [
    ("udata",      "Ice velocity (x-component)",       "m yr-1", "land_ice_x_velocity"),
    ("udata_mask", "Velocity data mask (x-component)", "1",      nothing),
    ("u_iszero",   "Velocity zero mask (x-component)", "1",      nothing),
]

const _V_GRID_VARS = [
    ("vdata",      "Ice velocity (y-component)",       "m yr-1", "land_ice_y_velocity"),
    ("vdata_mask", "Velocity data mask (y-component)", "1",      nothing),
    ("v_iszero",   "Velocity zero mask (y-component)", "1",      nothing),
]

# ── Helpers ───────────────────────────────────────────────────────────

const _H_FIELD_MAP = Dict(
    "thickness"         => :h_clip,
    "surface"           => :s_clip,
    "bed"               => :b_clip,
    "h_mask"            => :mask_clip,
    "basinID"           => :basinID_clip,
    "accumulation_data" => :a_Arthern_clip,
    "dhdt_data"         => :dhdt_clip,
    "dhdt_acc_mask"     => :dhdtAccDataMask_clip,
)

const _U_FIELD_MAP = Dict(
    "udata"      => :uData_clip,
    "udata_mask" => :uDataMaskFull_clip,
    "u_iszero"   => :uiszero_clip,
)

const _V_FIELD_MAP = Dict(
    "vdata"      => :vData_clip,
    "vdata_mask" => :vDataMaskFull_clip,
    "v_iszero"   => :viszero_clip,
)

"""
    _collect_grid_fields(grid, var_meta, field_map)

Return `Vector{Tuple{String, String, String, Union{String,Nothing}, AbstractArray}}`
for every variable whose backing field exists on `grid`.
"""
function _collect_grid_fields(grid, var_meta, field_map)
    out = Vector{Tuple{String, String, String, Union{String,Nothing}, AbstractArray}}()
    for (name, long_name, units, std_name) in var_meta
        field = field_map[name]
        if haskey(grid, field)
            push!(out, (name, long_name, units, std_name, getproperty(grid, field)))
        end
    end
    return out
end

_collect_h_fields(Gh) = _collect_grid_fields(Gh, _H_GRID_VARS, _H_FIELD_MAP)
_collect_u_fields(Gu) = _collect_grid_fields(Gu, _U_GRID_VARS, _U_FIELD_MAP)
_collect_v_fields(Gv) = _collect_grid_fields(Gv, _V_GRID_VARS, _V_FIELD_MAP)

"""Collect all fields across H, U, V grids (used by binary writer)."""
function _collect_all_fields(Gh, Gu, Gv)
    return vcat(_collect_h_fields(Gh), _collect_u_fields(Gu), _collect_v_fields(Gv))
end

# ── Binary writing ────────────────────────────────────────────────────

"""
    write_binary_files(Gh, Gu, Gv, output_path)

Write all clipped data to raw binary files (Float64, column-major).
Each 2-D field is stored in a separate `.bin` file; 3-D temperature is stored
as `temps.bin` and the sigma coordinate as `sigma_grid.bin`.
"""
function write_binary_files(Gh, Gu, Gv, output_path)
    mkpath(output_path)

    field_pairs = _collect_all_fields(Gh, Gu, Gv)

    for (name, _, _, _, data) in field_pairs
        filepath = joinpath(output_path, name * ".bin")
        open(filepath, "w") do io
            write(io, Float64.(data))
        end
    end

    # Temperature & sigma
    if haskey(Gh, :levels) && haskey(Gh.levels, :temperature_clip)
        filepath = joinpath(output_path, "temps.bin")
        open(filepath, "w") do io
            write(io, Float64.(Gh.levels.temperature_clip))
        end
    end
    if haskey(Gh, :levels) && haskey(Gh.levels, :sigma_full)
        filepath = joinpath(output_path, "sigma_grid.bin")
        open(filepath, "w") do io
            write(io, Float64.(Gh.levels.sigma_full))
        end
    end

    @info "Binary output written" path = output_path
end

# ── NetCDF writing ────────────────────────────────────────────────────

"""
    write_netcdf_file(Gh, Gu, Gv, output_path;
                      filename = "wavi_input.nc")

Write all clipped data to a single CF-compliant NetCDF-4 file with full
metadata (units, long names, coordinate variables, global attributes).

# Arguments
- `Gh`, `Gu`, `Gv`: Grid NamedTuples produced by `setup_wavi_data`.
- `output_path`: Directory in which the file will be created.
- `filename`: Name of the NetCDF file (default `"wavi_input.nc"`).
- `overwrite`: If `true`, overwrite an existing file (default `false`).

# Metadata included
- **Global attributes**: Conventions, title, source software, creation date,
  grid spacing, domain origin, basin IDs used.
- **Coordinate variables**: `x`, `y` (from the clipped H-grid), and `sigma`
  (if temperature data is present).
- **Per-variable attributes**: `long_name`, `units`, `standard_name` (where
  applicable), and `_FillValue = -9999.0`.
"""
function write_netcdf_file(Gh, Gu, Gv, output_path;
                           filename::String = "wavi_input.nc")
    mkpath(output_path)
    filepath = joinpath(output_path, filename)



    # ── Determine dimensions ──────────────────────────────────────────
    nx_h, ny_h = Gh.nx_clip, Gh.ny_clip

    # U-grid and V-grid are staggered and may differ in size
    u_fields = _collect_u_fields(Gu)
    v_fields = _collect_v_fields(Gv)
    nx_u, ny_u = isempty(u_fields) ? (nx_h, ny_h) : size(u_fields[1][5])
    nx_v, ny_v = isempty(v_fields) ? (nx_h, ny_h) : size(v_fields[1][5])

    has_temp = haskey(Gh, :levels) &&
               haskey(Gh.levels, :temperature_clip) &&
               haskey(Gh.levels, :sigma_full)
    nz = has_temp ? length(Gh.levels.sigma_full) : 0

    # ── Coordinate vectors ────────────────────────────────────────────
    x_h = Gh.xx_clip[:, 1]   # H-grid x coords
    y_h = Gh.yy_clip[1, :]   # H-grid y coords

    # ── Create file ───────────────────────────────────────────────────
    ds = NCDataset(filepath, "c")

    try
        # — Global attributes —
        ds.attrib["Conventions"]  = "CF-1.8"
        ds.attrib["title"]        = "WAVI model input data"
        ds.attrib["source"]       = "WAVIConstructor.jl"
        ds.attrib["history"]      = "Created $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))"
        ds.attrib["grid_spacing_m"] = Gh.dx
        ds.attrib["domain_origin_x_m"] = Gh.x0_clip
        ds.attrib["domain_origin_y_m"] = Gh.y0_clip
        ds.attrib["nx"] = nx_h
        ds.attrib["ny"] = ny_h
        if haskey(Gh, :basinID_clip)
            basin_ids_present = sort(unique(filter(!isnan, Gh.basinID_clip[Gh.mask_clip])))
            basin_str = join(Int.(basin_ids_present), ", ")
            ds.attrib["basin_ids"] = basin_str
        end
        ds.attrib["projection"] = "Antarctic Polar Stereographic (EPSG:3031)"
        ds.attrib["fill_value_convention"] = "-9999.0 used for masked / missing values"

        # — H-grid dimensions & coordinate variables —
        defDim(ds, "x", nx_h)
        defDim(ds, "y", ny_h)

        xvar = defVar(ds, "x", Float64, ("x",))
        xvar.attrib["long_name"] = "Easting (H-grid)"
        xvar.attrib["units"]     = "m"
        xvar.attrib["standard_name"] = "projection_x_coordinate"
        xvar[:] = x_h

        yvar = defVar(ds, "y", Float64, ("y",))
        yvar.attrib["long_name"] = "Northing (H-grid)"
        yvar.attrib["units"]     = "m"
        yvar.attrib["standard_name"] = "projection_y_coordinate"
        yvar[:] = y_h

        # — U-grid dimensions (staggered in x) —
        if !isempty(u_fields)
            defDim(ds, "xu", nx_u)
            defDim(ds, "yu", ny_u)
            xu_var = defVar(ds, "xu", Float64, ("xu",))
            xu_var.attrib["long_name"] = "Easting (U-grid, staggered in x)"
            xu_var.attrib["units"]     = "m"
            xu_var[:] = haskey(Gu, :xx) ? Gu.xx[1:nx_u, 1] : collect(1.0:nx_u)
            yu_var = defVar(ds, "yu", Float64, ("yu",))
            yu_var.attrib["long_name"] = "Northing (U-grid)"
            yu_var.attrib["units"]     = "m"
            yu_var[:] = haskey(Gu, :yy) ? Gu.yy[1, 1:ny_u] : collect(1.0:ny_u)
        end

        # — V-grid dimensions (staggered in y) —
        if !isempty(v_fields)
            defDim(ds, "xv", nx_v)
            defDim(ds, "yv", ny_v)
            xv_var = defVar(ds, "xv", Float64, ("xv",))
            xv_var.attrib["long_name"] = "Easting (V-grid)"
            xv_var.attrib["units"]     = "m"
            xv_var[:] = haskey(Gv, :xx) ? Gv.xx[1:nx_v, 1] : collect(1.0:nx_v)
            yv_var = defVar(ds, "yv", Float64, ("yv",))
            yv_var.attrib["long_name"] = "Northing (V-grid, staggered in y)"
            yv_var.attrib["units"]     = "m"
            yv_var[:] = haskey(Gv, :yy) ? Gv.yy[1, 1:ny_v] : collect(1.0:ny_v)
        end

        # — Sigma dimension —
        if has_temp
            defDim(ds, "sigma", nz)
            svar = defVar(ds, "sigma", Float64, ("sigma",))
            svar.attrib["long_name"]      = "Normalised depth coordinate (0=surface, 1=base)"
            svar.attrib["units"]          = "1"
            svar.attrib["positive"]       = "down"
            svar.attrib["standard_name"]  = "land_ice_sigma_coordinate"
            svar[:] = Gh.levels.sigma_full
        end

        fill_val = -9999.0

        # — H-grid variables (x, y) —
        for (name, long_name, units, std_name, data) in _collect_h_fields(Gh)
            v = defVar(ds, name, Float64, ("x", "y"), fillvalue = fill_val)
            v.attrib["long_name"] = long_name
            v.attrib["units"]     = units
            v.attrib["grid"]      = "h_grid"
            if std_name !== nothing
                v.attrib["standard_name"] = std_name
            end
            v[:, :] = Float64.(data)
        end

        # — U-grid variables (xu, yu) —
        for (name, long_name, units, std_name, data) in u_fields
            v = defVar(ds, name, Float64, ("xu", "yu"), fillvalue = fill_val)
            v.attrib["long_name"] = long_name
            v.attrib["units"]     = units
            v.attrib["grid"]      = "u_grid"
            if std_name !== nothing
                v.attrib["standard_name"] = std_name
            end
            v[:, :] = Float64.(data)
        end

        # — V-grid variables (xv, yv) —
        for (name, long_name, units, std_name, data) in v_fields
            v = defVar(ds, name, Float64, ("xv", "yv"), fillvalue = fill_val)
            v.attrib["long_name"] = long_name
            v.attrib["units"]     = units
            v.attrib["grid"]      = "v_grid"
            if std_name !== nothing
                v.attrib["standard_name"] = std_name
            end
            v[:, :] = Float64.(data)
        end

        # — 3-D temperature (x, y, sigma) —
        if has_temp
            tvar = defVar(ds, "temperature", Float64, ("x", "y", "sigma"),
                          fillvalue = fill_val)
            tvar.attrib["long_name"]      = "Ice temperature"
            tvar.attrib["units"]          = "K"
            tvar.attrib["standard_name"]  = "land_ice_temperature"
            tvar[:, :, :] = Float64.(Gh.levels.temperature_clip)
        end

    finally
        close(ds)
    end

    @info "NetCDF output written" filepath
    return filepath
end

# ── Unified dispatcher ────────────────────────────────────────────────

"""
    write_output(Gh, Gu, Gv, output_path; format=:bin, nc_filename="wavi_input.nc")

Write WAVI model input data in the requested format.

# Arguments
- `format`: One of `:bin`, `:netcdf`, or `:both`.
  - `:bin`    — raw Float64 binary files (one per variable).
  - `:netcdf` — single CF-compliant NetCDF-4 file.
  - `:both`   — write both formats side by side.
- `nc_filename`: Name of the NetCDF file (default `"wavi_input.nc"`).

"""
function write_output(Gh, Gu, Gv, output_path;
                      format::Symbol = :bin,
                      nc_filename::String = "wavi_input.nc",
)
    if format ∉ (:bin, :netcdf, :both)
        error("Unknown output format :$format — expected :bin, :netcdf, or :both")
    end

    if format in (:bin, :both)
        write_binary_files(Gh, Gu, Gv, output_path)
    end

    if format in (:netcdf, :both)
        write_netcdf_file(Gh, Gu, Gv, output_path;
                          filename = nc_filename)
    end
end

end # module OutputWriting
