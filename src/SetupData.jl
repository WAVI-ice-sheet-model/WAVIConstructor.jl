# WAVIConstructor.jl - Setup Data for WAVI
# Julia port of get_setup_data.m

module SetupData

using Interpolations

using WAVIConstructor.InitBedMachine: init_bedmachine
using WAVIConstructor.DomainSelection: select_domain_wavi

export setup_wavi_data

"""
    setup_wavi_data(params; output_path="outputs", edge=3)

Complete setup workflow for WAVI data preparation.
Julia port of MATLAB get_setup_data function.

# Arguments
- `params`: NamedTuple with all WAVI parameters
- `output_path`: Directory to write binary output files (default: "outputs")
- `edge`: Edge padding for domain clipping (default: 3)

# Returns
- Modified grid structures (Gh, Gu, Gv, Gc) with clipped data

# Notes
- Calls init_bedmachine to initialize grids
- Calls select_domain_wavi to select domain
- Clips domain to remove unused basins
- Replaces NaNs with -9999 outside masks
- Extrapolates temperature to include surface (sigma=0) and base (sigma=1)
- Writes binary files for Julia inversion
"""
function setup_wavi_data(params; output_path="outputs", edge=3)
    # Validate subsampling parameters
    sub_samp = get(params, :sub_samp, 8)
    sub_samp_index_x = get(params, :sub_samp_index_x, 0)
    sub_samp_index_y = get(params, :sub_samp_index_y, 0)
    
    if sub_samp_index_x >= 2 * sub_samp
        error("The subsamp index has to be less than twice the subSamp: exiting")
    end
    
    # Initialise grids using relevant datasets
    start_data = get(params, :start_data, "BEDMACHINEV3")
    if uppercase(start_data) == "BEDMACHINEV3"
        Gh, Gu, Gv, Gc = init_bedmachine(params)
    else
        error("Starting Dataset is not defined")
    end
    
    # Select domain using params.basins
    Gh, Gu, Gv, Gc = select_domain_wavi(Gh, Gu, Gv, Gc, params)
    
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
    Gu = merge(Gu, (
        uDataMaskFull_clip = Gu.uDataMask_clip .& Gu.mask_clip,
    ))
    Gv = merge(Gv, (
        vDataMaskFull_clip = Gv.vDataMask_clip .& Gv.mask_clip,
    ))
    
    # Extrapolate temperature to include surface (sigma=0) and base (sigma=1)
    Gh = extrapolate_temperature(Gh, I_clip_min, I_clip_max, J_clip_min, J_clip_max)
    
    # Write binary files
    write_binary_files(Gh, Gu, Gv, output_path)
    
    return Gh, Gu, Gv, Gc
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
    write_binary_files(Gh, Gu, Gv, output_path)

Write all clipped data to binary files for Julia inversion.
"""
function write_binary_files(Gh, Gu, Gv, output_path)
    # Create output directory if it doesn't exist
    if !isdir(output_path)
        mkpath(output_path)
    end
    
    # List of files to write (using Vector{Tuple{String, Any}} to accommodate 2D and 3D arrays)
    files_to_write = Tuple{String, Any}[
        ("thickness.bin", Gh.h_clip),
        ("surface.bin", Gh.s_clip),
        ("bed.bin", Gh.b_clip),
        ("h_mask.bin", Float64.(Gh.mask_clip)),
        ("u_iszero.bin", Float64.(Gu.uiszero_clip)),
        ("v_iszero.bin", Float64.(Gv.viszero_clip)),
        ("basinID.bin", Gh.basinID_clip),
        ("accumulation_data.bin", Gh.a_Arthern_clip),
        ("dhdt_data.bin", Gh.dhdt_clip),
        ("dhdt_acc_mask.bin", Float64.(Gh.dhdtAccDataMask_clip)),
        ("udata.bin", Gu.uData_clip),
        ("vdata.bin", Gv.vData_clip),
        ("udata_mask.bin", Float64.(Gu.uDataMaskFull_clip)),
        ("vdata_mask.bin", Float64.(Gv.vDataMaskFull_clip))
    ]
    
    # Add temperature and sigma files if available
    if haskey(Gh, :levels) && haskey(Gh.levels, :temperature_clip)
        push!(files_to_write, ("temps.bin", Gh.levels.temperature_clip))
    end
    if haskey(Gh, :levels) && haskey(Gh.levels, :sigma_full)
        push!(files_to_write, ("sigma_grid.bin", Gh.levels.sigma_full))
    end
    
    for (filename, data) in files_to_write
        filepath = joinpath(output_path, filename)
        if !isfile(filepath)
            open(filepath, "w") do io
                write(io, Float64.(data))
            end
        end
    end
end

end # module SetupData
