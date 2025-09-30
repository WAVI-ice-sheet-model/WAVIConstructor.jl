#!/usr/bin/env julia
"""
Example script demonstrating how to use WAVIConstructor.DataLoading
to read various Antarctic ice sheet datasets.

This script loads data from the Data/ directory and prints basic information
about each dataset.
"""

using WAVIConstructor.DataLoading

# Define data directory
data_dir = joinpath(@__DIR__, "..", "Data")

println("Loading Antarctic ice sheet datasets using WAVIConstructor.jl")
println("=" ^ 60)

# 1. Load ALBMAP data
println("\n1. Loading ALBMAP data...")
albmap_file = joinpath(data_dir, "ALBMAPv1.nc")
if isfile(albmap_file)
    Gh = get_albmap(albmap_file)
    println("   ALBMAP loaded successfully!")
    println("   Grid size: $(Gh[:nx]) x $(Gh[:ny])")
    println("   Ice thickness range: $(minimum(Gh[:h])) to $(maximum(Gh[:h])) m")
    println("   Surface elevation range: $(minimum(Gh[:s])) to $(maximum(Gh[:s])) m")
else
    println("   ALBMAP file not found: $albmap_file")
end

# 2. Load BedMachine data
println("\n2. Loading BedMachine Antarctica v3 data...")
bedmachine_file = joinpath(data_dir, "BedMachineAntarctica-v3.nc")
if isfile(bedmachine_file)
    bed, x, y, geoid, mask, thickness, surface, firn = get_bedmachine(bedmachine_file)
    println("   BedMachine loaded successfully!")
    println("   Grid dimensions: $(size(bed))")
    println("   Bed elevation range: $(minimum(bed)) to $(maximum(bed)) m")
    println("   Ice thickness range: $(minimum(thickness)) to $(maximum(thickness)) m")
else
    println("   BedMachine file not found: $bedmachine_file")
end

# 3. Load BISICLES temperature data
println("\n3. Loading BISICLES temperature data...")
bisicles_file = joinpath(data_dir, "antarctica-bisicles-xyzT-8km.nc")
if isfile(bisicles_file)
    sigma, x, y, z, temps = get_bisicles_temps(bisicles_file)
    println("   BISICLES temperature data loaded successfully!")
    println("   Temperature grid: $(size(temps))")
    println("   Temperature range: $(minimum(temps)) to $(maximum(temps)) K")
    println("   Vertical levels: $(length(sigma))")
else
    println("   BISICLES file not found: $bisicles_file")
end

# 4. Load MEASURES ice velocity data
println("\n4. Loading MEASURES ice velocity data...")
measures_file = joinpath(data_dir, "Antarctica_ice_velocity_2014_2015_1km_v01.nc")
if isfile(measures_file)
    xx_v, yy_v, VX, VY = get_measures_velocities(measures_file)
    println("   MEASURES velocity data loaded successfully!")
    println("   Velocity grid: $(size(VX))")
    println("   Velocity magnitude range: $(minimum(sqrt.(VX.^2 .+ VY.^2))) to $(maximum(sqrt.(VX.^2 .+ VY.^2))) m/yr")
else
    println("   MEASURES file not found: $measures_file")
end

# 5. Load Smith et al. dh/dt data (example with one file)
println("\n5. Loading Smith et al. dh/dt data...")
smith_dir = joinpath(data_dir, "Smith_2020_dhdt")
smith_files = filter(f -> endswith(f, ".tif"), readdir(smith_dir; join=true))
if !isempty(smith_files)
    # Load the first .tif file as an example
    example_file = smith_files[1]
    xx, yy, dhdt = get_smith_dhdt(example_file)
    println("   Smith dh/dt data loaded successfully from: $(basename(example_file))")
    println("   dh/dt grid: $(size(dhdt))")
    println("   Elevation change range: $(minimum(dhdt)) to $(maximum(dhdt)) m/yr")
else
    println("   No Smith dh/dt .tif files found in: $smith_dir")
end

println("\n" * "=" ^ 60)
println("Example completed! All datasets loaded successfully.")
