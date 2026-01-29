# WAVIConstructor.jl Quick Reference Guide

# ============================================================================
# INSTALLATION
# ============================================================================

# From Julia REPL:
# using Pkg
# Pkg.add(PackageSpec(url="https://github.com/WAVI-ice-sheet-model/WAVIConstructor.jl.git", rev="main"))

# ============================================================================
# BASIC USAGE
# ============================================================================

using WAVIConstructor

# Minimal example (uses all defaults):
params = Dict(
    :bedmachine_file => "Data/BedMachineAntarctica-v3.nc",
    :dx => 10000.0
)
Gh, Gu, Gv, Gc = setup_wavi_data(params)

# ============================================================================
# MAIN FUNCTIONS
# ============================================================================

# 1. init_bedmachine(params)
#    - Loads BedMachine data and creates computational grids
#    - Loads additional datasets (velocity, temperature, accumulation)
#    - Returns: Gh, Gu, Gv, Gc (NamedTuples for each grid)

# 2. select_domain_wavi(Gh, Gu, Gv, Gc, params)
#    - Selects domain based on specified drainage basins
#    - Creates masks and identifies domain boundaries
#    - Returns: Updated Gh, Gu, Gv, Gc with mask information

# 3. setup_wavi_data(params)
#    - Complete workflow: init → select → clip → write files
#    - Calls init_bedmachine and select_domain_wavi internally
#    - Writes binary files for WAVI
#    - Returns: Processed Gh, Gu, Gv, Gc

# ============================================================================
# PARAMETER OPTIONS
# ============================================================================

params = Dict(
    # Required parameters:
    :bedmachine_file => "Data/BedMachineAntarctica-v3.nc",  # Path to BedMachine file
    :dx => 10000.0,                                         # Grid spacing (m)
    
    # Domain selection:
    :basins => [1, 2, 3],           # Zwally basin IDs (default: 1:27 for all)
    
    # Dataset selection:
    :dhdt_data => "Smith",          # Elevation change: "Smith" or nothing
    :smith_dhdt_dir => "Data/Smith_2020_dhdt",  # Smith data directory
    
    :veloc_data => "Measures_2016/17",  # Options:
                                        # "Measures_1" (from .mat)
                                        # "Measures_2016/17"
                                        # "Measures_2014/15"
                                        # "Measures_phase_v1"
    
    :temps => "BISICLES_8km",      # Options:
                                   # "Frank"
                                   # "BISICLES_8km"
                                   # "BISICLES_1km"
    
    # Output:
    :output_path => "wavi_input",   # Output directory for binary files
    
    # Processing options:
    :clip_edge_padding => 3,        # Grid cells to pad edges when clipping
    
    # Physical parameters:
    :density_ice => 918.0,          # Ice density (kg/m³)
    :density_ocean => 1028.0,       # Ocean density (kg/m³)
    :min_thick => 50.0              # Minimum ice thickness (m)
)

# ============================================================================
# GRID STRUCTURE (NamedTuples returned by functions)
# ============================================================================

# Gh (H-grid) - Main grid with scalar fields
# Key fields:
#   - h: Ice thickness (m)
#   - b: Bed elevation (m)
#   - s: Surface elevation (m)
#   - a_Arthern: Accumulation rate (ice equivalent)
#   - basin_id: Zwally basin IDs
#   - dhdt: Elevation change rate (m/yr)
#   - temperature: 3D temperature field
#   - mask: Boolean mask of domain points
#   - f: Linear indices of domain points
#   - n: Number of domain points
#   - xx, yy: Coordinate grids

# Gu (U-grid) - X-velocity grid
# Key fields:
#   - u: X-component of velocity (m/yr)
#   - uDataMask: Mask where velocity data is valid
#   - mask: Boolean mask of domain points
#   - rockEdges: Edges at rock boundaries
#   - halo: Boundary/halo points

# Gv (V-grid) - Y-velocity grid
# Key fields:
#   - v: Y-component of velocity (m/yr)
#   - vDataMask: Mask where velocity data is valid
#   - mask: Boolean mask of domain points
#   - rockEdges: Edges at rock boundaries
#   - halo: Boundary/halo points

# Gc (C-grid) - Corner/cell-centre grid
# Key fields:
#   - mask: Boolean mask of domain points
#   - f: Linear indices of domain points
#   - n: Number of domain points

# ============================================================================
# COMMON WORKFLOWS
# ============================================================================

# Workflow 1: Quick setup for a specific region
params = Dict(
    :bedmachine_file => "Data/BedMachineAntarctica-v3.nc",
    :dx => 10000.0,
    :basins => [1, 2, 3],  # Pine Island region
    :output_path => "my_wavi_run"
)
Gh, Gu, Gv, Gc = setup_wavi_data(params)

# Workflow 2: Initialize and inspect before writing files
Gh, Gu, Gv, Gc = init_bedmachine(params)
# ... inspect/plot data ...
Gh, Gu, Gv, Gc = select_domain_wavi(Gh, Gu, Gv, Gc, params)
# ... more inspection ...
Gh, Gu, Gv, Gc = setup_wavi_data(params)

# Workflow 3: Use data in Julia without writing files
Gh, Gu, Gv, Gc = init_bedmachine(params)
Gh, Gu, Gv, Gc = select_domain_wavi(Gh, Gu, Gv, Gc, params)
# Now use Gh, Gu, Gv, Gc directly in your analysis

# ============================================================================
# DATA LOADING (Low-level access)
# ============================================================================

# If you need direct access to data loading functions:
using WAVIConstructor.DataLoading

# Load individual datasets:
bed, x, y, geoid, mask, thickness, surface, firn = get_bedmachine("Data/BedMachineAntarctica-v3.nc")
xx, yy, vx, vy = get_measures_velocities("Data/Antarctic_ice_velocity_2016_2017_1km_v01.nc")
xx, yy, dhdt = get_smith_dhdt(grnd_file="Data/Smith_2020_dhdt/ais_grounded.tif")
sigma, x, y, z, temps = get_bisicles_temps("Data/antarctica-bisicles-xyzT-8km.nc")

# ============================================================================
# ZWALLY DRAINAGE BASIN IDs
# ============================================================================

# Basin 1: Pine Island Glacier
# Basin 2: Thwaites Glacier
# Basin 3: Haynes Glacier
# Basin 4-6: Marie Byrd Land
# Basin 7-12: Ross Ice Shelf drainage
# Basin 13-17: Transantarctic Mountains/Victoria Land
# Basin 18-20: Wilkes Land
# Basin 21-23: Aurora/Vincennes Bays
# Basin 24-27: East Antarctic interior

# Use basins 1:27 for full Antarctica

# ============================================================================
# TROUBLESHOOTING
# ============================================================================

# Issue: Package not found
# Solution: Make sure you've installed with the full URL:
#   Pkg.add(PackageSpec(url="https://github.com/WAVI-ice-sheet-model/WAVIConstructor.jl.git", rev="main"))

# Issue: Data files not found
# Solution: Check that paths in params Dict point to correct locations
#   Use absolute paths if needed

# Issue: Out of memory
# Solution: Use coarser grid spacing (:dx => 20000.0 or higher)
#   Or select fewer basins

# Issue: Slow performance
# Solution: Start with a small domain (1-3 basins) and coarser grid
#   Process full Antarctica only when necessary
