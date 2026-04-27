# WAVIConstructor.jl - Setup Data for WAVI
# Julia port of get_setup_data.m

module SetupData

using Interpolations

using WAVIConstructor.InitBedMachine: init_bedmachine
using WAVIConstructor.DomainSelection: select_domain_wavi
using WAVIConstructor.ParamHelpers: ConstructorParams, to_dict
using WAVIConstructor.DataSources: BedMachineV3, NoData
using WAVIConstructor.OutputWriting: write_output

export setup_wavi_data

"""
    setup_wavi_data(params; output_path="outputs", edge=3, output_format=:bin)

Complete setup workflow for WAVI data preparation.
Julia port of MATLAB get_setup_data function.

# Arguments
- `params`: NamedTuple with all WAVI parameters
- `output_path`: Directory to write output files (default: "outputs")
- `edge`: Edge padding for domain clipping (default: 3)
- `output_format`: Output format — `:bin`, `:netcdf`, or `:both` (default: `:netcdf`)
- `diagnostic_plots`: If `true`, save heatmap PNGs of key fields after `init_bedmachine`
  (pre-domain-selection) and after `select_domain_wavi` into `<output_path>/diagnostic_plots/`.
  Useful for diagnosing edge effects introduced by domain selection. (default: `false`)

# Returns
- Modified grid structures (Gh, Gu, Gv, Gc) with clipped data

# Notes
- Calls init_bedmachine to initialize grids
- Calls select_domain_wavi to select domain
- Clips domain to remove unused basins
- Replaces NaNs with -9999 outside masks
- Extrapolates temperature to include surface (sigma=0) and base (sigma=1)
- Writes output files for Julia inversion in the requested format(s)
"""
function setup_wavi_data(params; output_path="outputs", edge=3, output_format=nothing, diagnostic_plots=false)
    # Validate subsampling parameters
    sub_samp = get(params, :sub_samp, 8)
    sub_samp_index_x = get(params, :sub_samp_index_x, 0)
    sub_samp_index_y = get(params, :sub_samp_index_y, 0)
    
    if sub_samp_index_x >= 2 * sub_samp
        error("The subsamp index has to be less than twice the subSamp: exiting")
    end
    
    # Initialise grids using relevant datasets
    bed_source = get(params, :bed_source, BedMachineV3())
    if bed_source isa NoData
        error("Bed topography source is required and cannot be NoData.")
    end
    Gh, Gu, Gv, Gc = init_bedmachine(params)

    if diagnostic_plots
        diag_dir = joinpath(output_path, "diagnostic_plots")
        @info "Saving pre-domain-selection diagnostic plots to $diag_dir"
        _save_diagnostic_plots(Gh, Gu, Gv, diag_dir, "01_post_init_bedmachine")
        Gh_pre, Gu_pre, Gv_pre = Gh, Gu, Gv
        @info """
        ── Diagnostic tip ──────────────────────────────────────────────────────
        Check 01_post_init_bedmachine__speed.png  BEFORE continuing.
        If velocity edge artefacts are already visible there, the cause is the
        MEaSUREs interpolation (NaN / extrapolated values at coverage gaps),
        NOT domain selection.  Domain selection only erodes the ice edge by
        one cell (C-grid ALL-4-ok requirement) which is expected MATLAB behaviour.
        ────────────────────────────────────────────────────────────────────────
        """
    end

    # Select domain using params.basins
    Gh, Gu, Gv, Gc = select_domain_wavi(Gh, Gu, Gv, Gc, params)

    if diagnostic_plots
        diag_dir = joinpath(output_path, "diagnostic_plots")
        @info "Saving post-domain-selection diagnostic plots to $diag_dir"
        _save_diagnostic_plots(Gh, Gu, Gv, diag_dir, "02_post_domain_selection")
        _compare_diagnostic_stages(
            Gh_pre, Gu_pre, Gv_pre,
            Gh,     Gu,     Gv,
            diag_dir,
            "01_post_init_bedmachine",
            "02_post_domain_selection",
        )
    end
    
    # Make mask of where we have accumulation and dhdt data
    Gh = merge(Gh, (
        dhdtAccDataMask = .!isnan.(Gh.dhdt) .& .!isnan.(Gh.a) .& Gh.mask .& Gh.aground,
    ))
    
    # Crop the mask to get rid of the basins we're not using
    # Find bounding box of mask
    I = findall(any(Gh.mask, dims=2)[:])
    J = findall(any(Gh.mask, dims=1)[:])
    
    if isempty(I) || isempty(J)
        @error "No valid mask points found" sum(Gh.mask) sum(Gh.ok) sum(Gh.ice) Gh.n unique(Gh.basin_id)
        error("No valid mask points found. Check that basins are valid and data loaded correctly.")
    end
    
    I_clip_min = max(1, minimum(I) - edge)
    I_clip_max = min(Gh.nx, maximum(I) + edge)
    J_clip_min = max(1, minimum(J) - edge)
    J_clip_max = min(Gh.ny, maximum(J) + edge)
    
    # Clip H-grid fields
    Gh_mask_clip = Gh.mask[I_clip_min:I_clip_max, J_clip_min:J_clip_max]
    
    # Find clipping indices for U, V, C grids
    Iu = findall(any(Gu.mask, dims=2)[:])
    Ju = findall(any(Gu.mask, dims=1)[:])
    Iu_clip_min = max(1, minimum(Iu) - edge)
    Iu_clip_max = min(Gu.nx, maximum(Iu) + edge)
    Ju_clip_min = max(1, minimum(Ju) - edge)
    Ju_clip_max = min(Gu.ny, maximum(Ju) + edge)
    
    Iv = findall(any(Gv.mask, dims=2)[:])
    Jv = findall(any(Gv.mask, dims=1)[:])
    Iv_clip_min = max(1, minimum(Iv) - edge)
    Iv_clip_max = min(Gv.nx, maximum(Iv) + edge)
    Jv_clip_min = max(1, minimum(Jv) - edge)
    Jv_clip_max = min(Gv.ny, maximum(Jv) + edge)
    
    Ic = findall(any(Gc.mask, dims=2)[:])
    Jc = findall(any(Gc.mask, dims=1)[:])
    Ic_clip_min = max(1, minimum(Ic) - edge)
    Ic_clip_max = min(Gc.nx, maximum(Ic) + edge)
    Jc_clip_min = max(1, minimum(Jc) - edge)
    Jc_clip_max = min(Gc.ny, maximum(Jc) + edge)
    
    # Clip all grid fields
    Gh = merge(Gh, (
        mask_clip = Gh_mask_clip,
        xx_clip = Gh.xx[I_clip_min:I_clip_max, J_clip_min:J_clip_max],
        yy_clip = Gh.yy[I_clip_min:I_clip_max, J_clip_min:J_clip_max],
        s_clip = Gh.s[I_clip_min:I_clip_max, J_clip_min:J_clip_max],
        h_clip = Gh.h[I_clip_min:I_clip_max, J_clip_min:J_clip_max],
        b_clip = Gh.b[I_clip_min:I_clip_max, J_clip_min:J_clip_max],
        a_clip = Gh.a[I_clip_min:I_clip_max, J_clip_min:J_clip_max],
        aground_clip = Gh.aground[I_clip_min:I_clip_max, J_clip_min:J_clip_max],
        dhdt_clip = Gh.dhdt[I_clip_min:I_clip_max, J_clip_min:J_clip_max],
        basinID_clip = Gh.basinID[I_clip_min:I_clip_max, J_clip_min:J_clip_max],
        a_Arthern_clip = Gh.a_Arthern[I_clip_min:I_clip_max, J_clip_min:J_clip_max],
        dhdtAccDataMask_clip = Gh.dhdtAccDataMask[I_clip_min:I_clip_max, J_clip_min:J_clip_max],
        nx_clip = size(Gh_mask_clip, 1),
        ny_clip = size(Gh_mask_clip, 2)
    ))
    
    # Calculate new corner point
    Gh = merge(Gh, (
        x0_clip = Gh.xx_clip[1, 1] - Gh.dx / 2,
        y0_clip = Gh.yy_clip[1, 1] - Gh.dy / 2
    ))
    
    # Clip U-grid fields
    Gu_mask_clip = Gu.mask[Iu_clip_min:Iu_clip_max, Ju_clip_min:Ju_clip_max]
    Gu = merge(Gu, (
        mask_clip = Gu_mask_clip,
        uData_clip = Gu.uData[Iu_clip_min:Iu_clip_max, Ju_clip_min:Ju_clip_max],
        uDataMask_clip = Gu.uDataMask[Iu_clip_min:Iu_clip_max, Ju_clip_min:Ju_clip_max]
    ))
    
    # Clip V-grid fields
    Gv_mask_clip = Gv.mask[Iv_clip_min:Iv_clip_max, Jv_clip_min:Jv_clip_max]
    Gv = merge(Gv, (
        mask_clip = Gv_mask_clip,
        vData_clip = Gv.vData[Iv_clip_min:Iv_clip_max, Jv_clip_min:Jv_clip_max],
        vDataMask_clip = Gv.vDataMask[Iv_clip_min:Iv_clip_max, Jv_clip_min:Jv_clip_max]
    ))
    
    # Clip C-grid fields
    Gc_mask_clip = Gc.mask[Ic_clip_min:Ic_clip_max, Jc_clip_min:Jc_clip_max]
    Gc = merge(Gc, (mask_clip = Gc_mask_clip,))
    
    # Find indices to clipped points
    Gh = merge(Gh, (
        f_clip = findall(Gh.mask_clip),
        n_clip = length(findall(Gh.mask_clip))
    ))
    Gu = merge(Gu, (
        f_clip = findall(Gu.mask_clip),
        n_clip = length(findall(Gu.mask_clip))
    ))
    Gv = merge(Gv, (
        f_clip = findall(Gv.mask_clip),
        n_clip = length(findall(Gv.mask_clip))
    ))
    Gc = merge(Gc, (
        f_clip = findall(Gc.mask_clip),
        n_clip = length(findall(Gc.mask_clip))
    ))
    
    # Replace NaNs with -9999 outside masks (Julia doesn't like NaNs)
    Gh = replace_nans_in_clipped_data(Gh)
    Gu = replace_nans_in_clipped_velocity(Gu)
    Gv = replace_nans_in_clipped_velocity(Gv)
    
    # Create uiszero and viszero masks
    Gu_uiszero = (Gu.halo .| Gu.rockEdges) .& Gu.mask
    Gv_viszero = (Gv.halo .| Gv.rockEdges) .& Gv.mask
    
    Gu = merge(Gu, (
        uiszero = Gu_uiszero,
        uiszero_clip = Gu_uiszero[Iu_clip_min:Iu_clip_max, Ju_clip_min:Ju_clip_max]
    ))
    Gv = merge(Gv, (
        viszero = Gv_viszero,
        viszero_clip = Gv_viszero[Iv_clip_min:Iv_clip_max, Jv_clip_min:Jv_clip_max]
    ))
    
    # Create full velocity masks (both data and mask)
    # MATLAB: uDataMask = (uData != 0) & Gu.mask  — velocity data AND ice mask
    Gu = merge(Gu, (
        uDataMaskFull_clip = Gu.uDataMask_clip .& Gu.mask_clip,
    ))
    Gv = merge(Gv, (
        vDataMaskFull_clip = Gv.vDataMask_clip .& Gv.mask_clip,
    ))
    
    # Extrapolate temperature to include surface (sigma=0) and base (sigma=1)
    Gh = extrapolate_temperature(Gh, I_clip_min, I_clip_max, J_clip_min, J_clip_max)
    
    # Determine output format: explicit kwarg > params dict > default (:bin)
    fmt = if output_format !== nothing
        output_format
    elseif isa(params, AbstractDict) && haskey(params, :output_format)
        params[:output_format]
    else
        :netcdf
    end

    # Write output files in the requested format(s)
    write_output(Gh, Gu, Gv, output_path; format=fmt)
    
    return Gh, Gu, Gv, Gc
end

"""
    setup_wavi_data(params::ConstructorParams; kwargs...)

Convenience method that accepts a ConstructorParams struct directly.
Converts it to a Dict and calls the main setup_wavi_data function.

# Example
```julia
params = minimal_constructor_params(temps="Frank")
Gh, Gu, Gv, Gc = setup_wavi_data(params)  # No need for to_dict!
```
"""
function setup_wavi_data(params::ConstructorParams; kwargs...)
    return setup_wavi_data(to_dict(params); kwargs...)
end

"""
    replace_nans_in_clipped_data(Gh)

Replace NaNs with -9999 outside masks in clipped H-grid data.
"""
function replace_nans_in_clipped_data(Gh)
    # Replace NaNs in various fields
    updates = Dict()
    
    if haskey(Gh, :basinID_clip)
        data = Gh.basinID_clip
        nan_mask = isnan.(data) .& .!Gh.mask_clip
        if any(nan_mask)
            data_new = copy(data)
            data_new[nan_mask] .= -9999.0
            updates[:basinID_clip] = data_new
        end
    end
    
    for field in [:h_clip, :s_clip, :b_clip, :a_clip, :dhdt_clip, :a_Arthern_clip]
        if haskey(Gh, field)
            data = getproperty(Gh, field)
            nan_mask = isnan.(data) .& .!Gh.mask_clip
            if any(nan_mask)
                data_new = copy(data)
                data_new[nan_mask] .= -9999.0
                updates[field] = data_new
            end
        end
    end
    
    return isempty(updates) ? Gh : merge(Gh, NamedTuple(updates))
end

"""
    replace_nans_in_clipped_velocity(G)

Replace NaNs with -9999 outside masks in clipped velocity grid data.
"""
function replace_nans_in_clipped_velocity(G)
    updates = Dict()
    
    # Replace NaNs in velocity data
    if haskey(G, :uData_clip)
        data = G.uData_clip
        mask = G.uDataMask_clip
        nan_mask = isnan.(data) .& .!mask
        if any(nan_mask)
            data_new = copy(data)
            data_new[nan_mask] .= -9999.0
            updates[:uData_clip] = data_new
        end
        
        # Also replace NaNs in mask itself
        mask_nan = isnan.(mask) .& .!mask
        if any(mask_nan)
            mask_new = copy(mask)
            mask_new[mask_nan] .= -9999.0
            updates[:uDataMask_clip] = mask_new
        end
    elseif haskey(G, :vData_clip)
        data = G.vData_clip
        mask = G.vDataMask_clip
        nan_mask = isnan.(data) .& .!mask
        if any(nan_mask)
            data_new = copy(data)
            data_new[nan_mask] .= -9999.0
            updates[:vData_clip] = data_new
        end
        
        # Also replace NaNs in mask itself
        mask_nan = isnan.(mask) .& .!mask
        if any(mask_nan)
            mask_new = copy(mask)
            mask_new[mask_nan] .= -9999.0
            updates[:vDataMask_clip] = mask_new
        end
    end
    
    return isempty(updates) ? G : merge(G, NamedTuple(updates))
end

"""
    extrapolate_temperature(Gh, I_min, I_max, J_min, J_max)

Extrapolate temperature to include surface (sigma=0) and base (sigma=1) if not present.
"""
function extrapolate_temperature(Gh, I_min, I_max, J_min, J_max)
    if !haskey(Gh, :levels) || !haskey(Gh.levels, :sigmas) || !haskey(Gh.levels, :temperature)
        return Gh
    end
    
    sigmas = Gh.levels.sigmas
    temperature = Gh.levels.temperature
    
    # Check if surface (0) and base (1) are already in sigmas
    if sigmas[1] == 0.0
        sigma_full = sigmas
        temps = temperature
    else
        # Add surface and base
        sigma_full = [0.0; sigmas; 1.0]
        
        # Extrapolate temperatures to new grid
        nz_full = length(sigma_full)
        temps = zeros(nz_full, Gh.nx, Gh.ny)
        
        for j in 1:Gh.ny
            for i in 1:Gh.nx
                # Create interpolation for this column with linear extrapolation
                # MATLAB's interp1 with 'linear' and 'extrap'
                # Use Interpolations.jl with extrapolation
                itp = LinearInterpolation(sigmas, temperature[:, i, j], 
                                         extrapolation_bc=Line())
                temps[:, i, j] = [itp(s) for s in sigma_full]
            end
        end
    end
    
    # Change array order to be nx, ny, nz (MATLAB permute([2 3 1]))
    temps_for_julia = permutedims(temps, (2, 3, 1))
    
    # Clip temperatures
    temps_for_julia_clip = temps_for_julia[I_min:I_max, J_min:J_max, :]
    
    return merge(Gh, (
        levels = merge(Gh.levels, (
            sigma_full = sigma_full,
            temperature_full = temps,
            temperature_clip = temps_for_julia_clip,
        )),
    ))
end

"""
    _save_diagnostic_plots(Gh, Gu, Gv, plot_dir, stage_label)

Save heatmap PNGs of key geophysical fields for diagnostic purposes.
Called by `setup_wavi_data` when `diagnostic_plots=true`.

- `stage_label`: prefix string used in filenames, e.g. `"01_post_init_bedmachine"`.

Fields plotted (where available):
  H-grid : thickness (`h`), bed (`b`), surface (`s`), accumulation (`a`/`a_Arthern`),
           ice mask (`mask`/`ice`), dhdt
  U-grid : u-velocity (`uData`)
  V-grid : v-velocity (`vData`)
  Derived: horizontal speed = √(u²+v²) interpolated to H-grid centroids
"""
function _save_diagnostic_plots(Gh, Gu, Gv, plot_dir, stage_label)
    # Lazy-load Plots — it is an optional dependency not listed in [deps].
    # Users must have Plots available in their environment (e.g. the test/docs
    # extras, or their own project). If it is missing we bail out with a clear
    # message rather than a hard precompilation error.
    Plots = try
        Base.require(Base.PkgId(Base.UUID("91a5bcdd-55d7-5caf-9e0b-520d859cae80"), "Plots"))
    catch
        @warn "diagnostic_plots=true requires the Plots package, but it is not " *
              "available in the current environment. Skipping diagnostic plots.\n" *
              "Add Plots to your project with `] add Plots` and retry."
        return
    end

    mkpath(plot_dir)

    heatmap_fn = Plots.heatmap
    savefig_fn = Plots.savefig

    function _heatmap_png(data, title_str, fname)
        vals = filter(isfinite, vec(data))
        clims = isempty(vals) ? (0.0, 1.0) : (minimum(vals), maximum(vals))
        p = Base.invokelatest(
            heatmap_fn,
            data';
            c = :viridis,
            clims = clims,
            title = title_str,
            xlabel = "x index",
            ylabel = "y index",
            aspect_ratio = :equal,
            size = (800, 700),
        )
        Base.invokelatest(savefig_fn, p, joinpath(plot_dir, "$(stage_label)__$(fname).png"))
    end

    # ── H-grid scalar fields ─────────────────────────────────────────
    for (field, label, fname) in [
            (:h,         "Ice thickness (m)",       "thickness"),
            (:b,         "Bed elevation (m)",        "bed"),
            (:s,         "Surface elevation (m)",    "surface"),
            (:dhdt,      "dh/dt (m/yr)",             "dhdt"),
            (:basin_id,  "Basin ID",                 "basin_id"),
            (:surfType,  "Surface type (1=grounded, 2=floating, 3=rock)", "surf_type"),
        ]
        if haskey(Gh, field)
            _heatmap_png(Float64.(getproperty(Gh, field)), label, fname)
        end
    end

    # ok mask — the physically meaningful "valid ice" criterion, shown explicitly
    # so it can be directly compared against post-selection Gh.mask
    if haskey(Gh, :ok)
        _heatmap_png(Float64.(Gh.ok), "Valid ice (ok = h>50m & ice)", "ok_mask")
    end

    # accumulation — may be stored as :a or :a_Arthern
    for field in (:a, :a_Arthern)
        if haskey(Gh, field)
            _heatmap_png(Float64.(getproperty(Gh, field)), "Accumulation (m/yr)", "accumulation")
            break
        end
    end

    # boolean masks → Float64 for heatmap
    for (field, label, fname) in [
            (:mask,   "Ice mask",       "mask"),
            (:ice,    "Ice flag",       "ice"),
            (:aground,"Grounded flag",  "aground"),
        ]
        if haskey(Gh, field)
            _heatmap_png(Float64.(getproperty(Gh, field)), label, fname)
        end
    end

    # ── Velocity fields ──────────────────────────────────────────────
    if haskey(Gu, :uData) && haskey(Gv, :vData)
        u = Gu.uData
        v = Gv.vData
        # Interpolate u (nx+1 × ny) and v (nx × ny+1) to H-grid (nx × ny)
        u_h = 0.5 .* (u[1:end-1, :] .+ u[2:end, :])
        v_h = 0.5 .* (v[:, 1:end-1] .+ v[:, 2:end])
        speed = sqrt.(u_h .^ 2 .+ v_h .^ 2)
        _heatmap_png(u_h,    "u-velocity (m/yr)",    "u_velocity")
        _heatmap_png(v_h,    "v-velocity (m/yr)",    "v_velocity")
        _heatmap_png(speed,  "Speed |u| (m/yr)",     "speed")
    elseif haskey(Gu, :uData)
        _heatmap_png(Gu.uData, "u-velocity (m/yr)", "u_velocity")
    elseif haskey(Gv, :vData)
        _heatmap_png(Gv.vData, "v-velocity (m/yr)", "v_velocity")
    end

    @info "  → $(length(readdir(plot_dir; join=false))) files now in $plot_dir"
end

"""
    _compare_diagnostic_stages(plot_dir, stage_before, stage_after)

Load the individual-stage heatmap arrays saved by `_save_diagnostic_plots` and
produce side-by-side comparison PNGs:

    [ before | after | difference ]

One PNG is written per matching field into `<plot_dir>/comparison/`.
Fields that only exist in one stage are skipped.
"""
function _compare_diagnostic_stages(Gh_before, Gu_before, Gv_before,
                                    Gh_after,  Gu_after,  Gv_after,
                                    plot_dir, stage_before, stage_after)
    Plots = try
        Base.require(Base.PkgId(Base.UUID("91a5bcdd-55d7-5caf-9e0b-520d859cae80"), "Plots"))
    catch
        @warn "Plots not available — skipping comparison plots."
        return
    end

    comp_dir = joinpath(plot_dir, "comparison")
    mkpath(comp_dir)

    heatmap_fn   = Plots.heatmap
    histogram_fn = Plots.histogram
    plot_fn      = Plots.plot
    savefig_fn   = Plots.savefig

    # Helper: pad or crop arr2 to match the size of arr1 (fill with NaN)
    function _match_size(arr1, arr2)
        if size(arr1) == size(arr2)
            return arr2
        end
        out = fill(NaN, size(arr1))
        r = min(size(arr1, 1), size(arr2, 1))
        c = min(size(arr1, 2), size(arr2, 2))
        out[1:r, 1:c] .= arr2[1:r, 1:c]
        return out
    end

    function _panel(arr, title_str, colormap, clims)
        Base.invokelatest(
            heatmap_fn, arr';
            c = colormap, clims = clims,
            title = title_str, xlabel = "x", ylabel = "y",
            aspect_ratio = :equal,
        )
    end

    function _compare_pair(before_arr, after_arr, field_name, label, colormap=:viridis)
        after_matched = _match_size(before_arr, Float64.(after_arr))

        all_vals = filter(isfinite, vcat(vec(before_arr), vec(after_matched)))
        isempty(all_vals) && return
        clims = (minimum(all_vals), maximum(all_vals))
        clims = clims[1] == clims[2] ? (clims[1] - 1.0, clims[2] + 1.0) : clims

        diff_arr  = after_matched .- before_arr
        diff_vals = filter(isfinite, vec(diff_arr))
        maxabs    = isempty(diff_vals) ? 1.0 : max(abs(minimum(diff_vals)), abs(maximum(diff_vals)))
        maxabs    = maxabs == 0.0 ? 1.0 : maxabs

        p1 = _panel(before_arr,   "Before\n($stage_before)",   colormap, clims)
        p2 = _panel(after_matched, "After\n($stage_after)",    colormap, clims)
        p3 = _panel(diff_arr, "Difference\n(after − before)", :RdBu, (-maxabs, maxabs))

        if isempty(diff_vals) || all(diff_vals .== 0)
            p4 = Base.invokelatest(plot_fn; title="Diff histogram\n(all zero)", legend=false)
        else
            p4 = Base.invokelatest(
                histogram_fn, diff_vals;
                bins=60, title="Diff histogram",
                xlabel="after − before", ylabel="count", legend=false,
            )
        end

        fig = Base.invokelatest(
            plot_fn, p1, p2, p3, p4;
            layout=(2, 2), size=(1600, 1000),
            plot_title="$label  |  $stage_before  vs  $stage_after",
        )
        out_path = joinpath(comp_dir, "compare__$(field_name).png")
        Base.invokelatest(savefig_fn, fig, out_path)
        @info "    📈  $(out_path)"
    end

    @info "Building stage-comparison plots in $comp_dir"

    # ── H-grid scalars ──────────────────────────────────────────────
    for (field, label, fname) in [
            (:h,        "Ice thickness (m)",     "thickness"),
            (:b,        "Bed elevation (m)",     "bed"),
            (:s,        "Surface elevation (m)", "surface"),
            (:dhdt,     "dh/dt (m/yr)",          "dhdt"),
            (:basin_id, "Basin ID",              "basin_id"),
            (:surfType, "Surface type (1=grounded, 2=floating, 3=rock)", "surf_type"),
        ]
        if haskey(Gh_before, field) && haskey(Gh_after, field)
            _compare_pair(
                Float64.(getproperty(Gh_before, field)),
                Float64.(getproperty(Gh_after,  field)),
                fname, label,
            )
        end
    end

    # accumulation
    acc_field_b = findfirst(f -> haskey(Gh_before, f), (:a, :a_Arthern))
    acc_field_a = findfirst(f -> haskey(Gh_after,  f), (:a, :a_Arthern))
    if acc_field_b !== nothing && acc_field_a !== nothing
        _compare_pair(
            Float64.(getproperty(Gh_before, (:a, :a_Arthern)[acc_field_b])),
            Float64.(getproperty(Gh_after,  (:a, :a_Arthern)[acc_field_a])),
            "accumulation", "Accumulation (m/yr)",
        )
    end

    # boolean masks
    # Before domain selection, Gh.mask is the raw BedMachine mask (ice + rock).
    # After domain selection, Gh.mask is restricted to ok cells (h>50m, ice only).
    # For a like-for-like comparison of "which cells are treated as active ice"
    # we compare Gh.ok (pre) vs Gh.mask (post).
    if haskey(Gh_before, :ok) && haskey(Gh_after, :mask)
        _compare_pair(
            Float64.(Gh_before.ok),
            Float64.(Gh_after.mask),
            "active_ice_mask",
            "Active ice mask  (ok before → mask after domain selection)",
        )
    end
    # Also show the raw mask change (BedMachine mask → domain mask) for reference
    for (field, label, fname) in [
            (:mask,    "Ice mask (BedMachine → domain)", "mask"),
            (:aground, "Grounded flag",                  "aground"),
        ]
        if haskey(Gh_before, field) && haskey(Gh_after, field)
            _compare_pair(
                Float64.(getproperty(Gh_before, field)),
                Float64.(getproperty(Gh_after,  field)),
                fname, label,
            )
        end
    end

    # ── Velocity ─────────────────────────────────────────────────────
    ub = vb = speedb = ua = va = speeda = nothing
    if haskey(Gu_before, :uData) && haskey(Gv_before, :vData)
        u = Gu_before.uData;  v = Gv_before.vData
        ub     = 0.5 .* (u[1:end-1, :] .+ u[2:end, :])
        vb     = 0.5 .* (v[:, 1:end-1] .+ v[:, 2:end])
        speedb = sqrt.(ub .^ 2 .+ vb .^ 2)
    end
    if haskey(Gu_after, :uData) && haskey(Gv_after, :vData)
        u = Gu_after.uData;  v = Gv_after.vData
        ua     = 0.5 .* (u[1:end-1, :] .+ u[2:end, :])
        va     = 0.5 .* (v[:, 1:end-1] .+ v[:, 2:end])
        speeda = sqrt.(ua .^ 2 .+ va .^ 2)
    end
    if speedb !== nothing && speeda !== nothing
        _compare_pair(Float64.(ub),     Float64.(ua),     "u_velocity", "u-velocity (m/yr)")
        _compare_pair(Float64.(vb),     Float64.(va),     "v_velocity", "v-velocity (m/yr)")
        _compare_pair(Float64.(speedb), Float64.(speeda), "speed",      "Speed |u| (m/yr)")
    end

    @info "Comparison plots written to $comp_dir"
end

end # module SetupData
