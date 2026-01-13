# WAVIConstructor.jl - Domain Selection for WAVI
# Julia port of SelectDomainWAVI.m

module DomainSelection

using ImageFiltering

export select_domain_wavi

"""
    select_domain_wavi(Gh, Gu, Gv, Gc, params)

Select domain and create masks for WAVI simulation based on basin selection.
Julia port of MATLAB SelectDomainWAVI function.

# Arguments
- `Gh`: H-grid structure with ice data
- `Gu`: U-grid structure for u-velocity
- `Gv`: V-grid structure for v-velocity  
- `Gc`: C-grid structure for corner points
- `params`: NamedTuple with parameters including `:basins` (array of basin IDs)

# Returns
- Modified grid structures (Gh, Gu, Gv, Gc) with masks and indices added

# Notes
- Creates global masks (all ice, not basin-restricted)
- Creates basin-specific masks based on params.basins
- Creates rock edge masks for velocity boundaries
- Creates halo masks for velocity exchange between basins
"""
function select_domain_wavi(Gh, Gu, Gv, Gc, params)
    basins = get(params, :basins, 1:27)
    
    # Find masks for points to include as ice on each grid
    # These global masks are not restricted to basins of interest
    Gc_globalMask = (Gh.ok[1:end-1, 1:end-1] .&
                     Gh.ok[1:end-1, 2:end] .&
                     Gh.ok[2:end, 1:end-1] .&
                     Gh.ok[2:end, 2:end])
    
    # Expand C-grid global mask to H-grid
    # MATLAB conv2([1 1]', [1 1], Gc_globalMask) expands from (nx-1, ny-1) to (nx, ny)
    # We need to OR adjacent cells to expand the mask
    gh_target_size = size(Gh.ok)  # Target size for H-grid mask
    
    # Expand by OR-ing: each H-grid cell is true if any adjacent C-grid cell is true
    Gh_globalMask = falses(gh_target_size)
    Gh_globalMask[1:end-1, 1:end-1] .|= Gc_globalMask
    Gh_globalMask[1:end-1, 2:end] .|= Gc_globalMask
    Gh_globalMask[2:end, 1:end-1] .|= Gc_globalMask
    Gh_globalMask[2:end, 2:end] .|= Gc_globalMask
    
    # Create mask for locations on the h-grid to include in model
    # MATLAB's ismember equivalent
    Gh_basinOK = [basin_id in basins for basin_id in Gh.basinID]
    
    # Find masks for points to include as ice on each grid (restricted to basins)
    Gc_mask = (Gh.ok[1:end-1, 1:end-1] .&
               Gh.ok[1:end-1, 2:end] .&
               Gh.ok[2:end, 1:end-1] .&
               Gh.ok[2:end, 2:end]) .&
              (Gh_basinOK[1:end-1, 1:end-1] .|
               Gh_basinOK[1:end-1, 2:end] .|
               Gh_basinOK[2:end, 1:end-1] .|
               Gh_basinOK[2:end, 2:end])
    
    # Expand C-grid mask to H-grid
    # MATLAB conv2([1 1]', [1 1], Gc_mask) expands from (nx-1, ny-1) to (nx, ny)
    # We need to OR adjacent cells to expand the mask
    gc_size = size(Gc_mask)
    gh_target_size = size(Gh.ok)  # Target size for H-grid mask
    
    # Expand by OR-ing: each H-grid cell is true if any adjacent C-grid cell is true
    Gh_mask = falses(gh_target_size)
    Gh_mask[1:end-1, 1:end-1] .|= Gc_mask
    Gh_mask[1:end-1, 2:end] .|= Gc_mask
    Gh_mask[2:end, 1:end-1] .|= Gc_mask
    Gh_mask[2:end, 2:end] .|= Gc_mask
    
    # Expand H-grid mask to U-grid and V-grid to match actual grid sizes
    # U-grid and V-grid sizes are determined by Gu.xx and Gv.xx
    gu_size = size(Gu.xx)
    gv_size = size(Gv.xx)
    gh_mask_size = size(Gh_mask)
    
    # For U-grid: should match size of Gu.xx 
    # MATLAB conv2([1 1]', 1, A) expands first dimension by 1
    if gu_size[1] == gh_mask_size[1] + 1 && gu_size[2] == gh_mask_size[2]
        # Standard case: U-grid has one more row than H-grid
        Gu_mask = [Gh_mask[1:end-1, :] .| Gh_mask[2:end, :]; Gh_mask[end:end, :]]
    else
        # Fall back to resizing
        Gu_mask = falses(gu_size)
        Gu_mask[1:min(gu_size[1], gh_mask_size[1]), 1:min(gu_size[2], gh_mask_size[2])] .= Gh_mask[1:min(gu_size[1], gh_mask_size[1]), 1:min(gu_size[2], gh_mask_size[2])]
    end
    
    # For V-grid: should match size of Gv.xx
    # MATLAB conv2(1, [1 1], A) expands second dimension by 1
    if gv_size[1] == gh_mask_size[1] && gv_size[2] == gh_mask_size[2] + 1
        # Standard case: V-grid has one more column than H-grid
        Gv_mask = [Gh_mask[:, 1:end-1] .| Gh_mask[:, 2:end] Gh_mask[:, end:end]]
    else
        # Fall back to resizing
        Gv_mask = falses(gv_size)
        Gv_mask[1:min(gv_size[1], gh_mask_size[1]), 1:min(gv_size[2], gh_mask_size[2])] .= Gh_mask[1:min(gv_size[1], gh_mask_size[1]), 1:min(gv_size[2], gh_mask_size[2])]
    end
    
    # Now compute rock edges on u- and v-grids (needs Gu_mask and Gv_mask to be defined)
    # Logical masks on u- and v-grids corresponding to edges of rock mask
    # MATLAB: diff([zeros(1,size(Gh.rock,2));Gh.rock;zeros(1,size(Gh.rock,2))],1,1)
    # This pads with zeros, then diff finds where values change
    rock_size = size(Gh.rock)
    rock_padded_u = [zeros(1, rock_size[2]); Gh.rock; zeros(1, rock_size[2])]
    Gu_rockEdges_temp = diff(rock_padded_u, dims=1) .!= 0
    # Resize to match Gu_mask size exactly
    Gu_rockEdges = falses(size(Gu_mask))
    Gu_rockEdges[1:min(size(Gu_rockEdges, 1), size(Gu_rockEdges_temp, 1)), 1:min(size(Gu_rockEdges, 2), size(Gu_rockEdges_temp, 2))] .= Gu_rockEdges_temp[1:min(size(Gu_rockEdges, 1), size(Gu_rockEdges_temp, 1)), 1:min(size(Gu_rockEdges, 2), size(Gu_rockEdges_temp, 2))]
    
    # MATLAB: diff([zeros(size(Gh.rock,1),1) Gh.rock zeros(size(Gh.rock,1),1)],1,2)
    rock_padded_v = [zeros(rock_size[1], 1) Gh.rock zeros(rock_size[1], 1)]
    Gv_rockEdges_temp = diff(rock_padded_v, dims=2) .!= 0
    # Resize to match Gv_mask size exactly
    Gv_rockEdges = falses(size(Gv_mask))
    Gv_rockEdges[1:min(size(Gv_rockEdges, 1), size(Gv_rockEdges_temp, 1)), 1:min(size(Gv_rockEdges, 2), size(Gv_rockEdges_temp, 2))] .= Gv_rockEdges_temp[1:min(size(Gv_rockEdges, 1), size(Gv_rockEdges_temp, 1)), 1:min(size(Gv_rockEdges, 2), size(Gv_rockEdges_temp, 2))]
    
    # Find indices to those points (MATLAB's find returns linear indices)
    Gh_f = LinearIndices(Gh_mask)[Gh_mask]
    Gu_f = LinearIndices(Gu_mask)[Gu_mask]
    Gv_f = LinearIndices(Gv_mask)[Gv_mask]
    Gc_f = LinearIndices(Gc_mask)[Gc_mask]
    
    # Number of points included in each mask
    Gh_n = length(Gh_f)
    Gu_n = length(Gu_f)
    Gv_n = length(Gv_f)
    Gc_n = length(Gc_f)
    
    # Logical mask corresponding to other basins, but not the current basin
    Gh_otherBasins = Gh_globalMask .& .!Gh_mask
    
    # Logical masks on u- and v-grids corresponding to halo
    # The halo corresponds to velocities exchanged with other basins
    # MATLAB: diff([zeros(1,size(Gh.otherBasins,2));Gh.otherBasins;zeros(1,size(Gh.otherBasins,2))],1,1)
    gh_size = size(Gh_otherBasins)
    otherBasins_padded_u = [zeros(1, gh_size[2]); Gh_otherBasins; zeros(1, gh_size[2])]
    otherBasins_diff_u_temp = diff(otherBasins_padded_u, dims=1) .!= 0
    # Resize to match Gu_mask size
    otherBasins_diff_u = falses(size(Gu_mask))
    otherBasins_diff_u[1:min(size(otherBasins_diff_u, 1), size(otherBasins_diff_u_temp, 1)), 1:min(size(otherBasins_diff_u, 2), size(otherBasins_diff_u_temp, 2))] .= otherBasins_diff_u_temp[1:min(size(otherBasins_diff_u, 1), size(otherBasins_diff_u_temp, 1)), 1:min(size(otherBasins_diff_u, 2), size(otherBasins_diff_u_temp, 2))]
    Gu_halo = Gu_mask .& .!Gu_rockEdges .& otherBasins_diff_u
    
    # MATLAB: diff([zeros(size(Gh.otherBasins,1),1) Gh.otherBasins zeros(size(Gh.otherBasins,1),1)],1,2)
    otherBasins_padded_v = [zeros(gh_size[1], 1) Gh_otherBasins zeros(gh_size[1], 1)]
    otherBasins_diff_v_temp = diff(otherBasins_padded_v, dims=2) .!= 0
    # Resize to match Gv_mask size
    otherBasins_diff_v = falses(size(Gv_mask))
    otherBasins_diff_v[1:min(size(otherBasins_diff_v, 1), size(otherBasins_diff_v_temp, 1)), 1:min(size(otherBasins_diff_v, 2), size(otherBasins_diff_v_temp, 2))] .= otherBasins_diff_v_temp[1:min(size(otherBasins_diff_v, 1), size(otherBasins_diff_v_temp, 1)), 1:min(size(otherBasins_diff_v, 2), size(otherBasins_diff_v_temp, 2))]
    Gv_halo = Gv_mask .& .!Gv_rockEdges .& otherBasins_diff_v
    
    # Update grid structures with new fields
    Gh = merge(Gh, (
        globalMask = Gh_globalMask,
        basinOK = Gh_basinOK,
        mask = Gh_mask,
        f = Gh_f,
        n = Gh_n,
        otherBasins = Gh_otherBasins
    ))
    
    Gu = merge(Gu, (
        rockEdges = Gu_rockEdges,
        mask = Gu_mask,
        f = Gu_f,
        n = Gu_n,
        halo = Gu_halo
    ))
    
    Gv = merge(Gv, (
        rockEdges = Gv_rockEdges,
        mask = Gv_mask,
        f = Gv_f,
        n = Gv_n,
        halo = Gv_halo
    ))
    
    Gc = merge(Gc, (
        globalMask = Gc_globalMask,
        mask = Gc_mask,
        f = Gc_f,
        n = Gc_n
    ))
    
    return Gh, Gu, Gv, Gc
end

end # module DomainSelection
