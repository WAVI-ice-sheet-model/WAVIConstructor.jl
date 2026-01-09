#!/usr/bin/env julia
"""
Test script for the newly implemented data loading functions
"""

using WAVIConstructor.DataLoading

println("Testing newly implemented data loading functions...")
println("=" ^ 60)

# Test 1: get_arthern_accumulation
println("\n1. Testing get_arthern_accumulation()...")
if isfile("Data/amsr_accumulation_map.txt")
    try
        aa_lat, aa_lon, aa_x, aa_y, aa_acc, aa_err = get_arthern_accumulation("Data/amsr_accumulation_map.txt")
        println("   ✓ Success! Loaded $(length(aa_acc)) data points")
        println("   ✓ Accumulation range: $(minimum(aa_acc)) to $(maximum(aa_acc)) m/yr (ice equivalent)")
        println("   ✓ X range: $(minimum(aa_x)) to $(maximum(aa_x)) m")
        println("   ✓ Y range: $(minimum(aa_y)) to $(maximum(aa_y)) m")
    catch e
        println("   ✗ Error: ", e)
    end
else
    println("   ⚠ Skipping: Data file not found (Data/amsr_accumulation_map.txt)")
end

# Test 2: get_zwally_basins
println("\n2. Testing get_zwally_basins()...")
if isfile("Data/DrainageBasins/ZwallyBasins.mat")
    try
        xx_zwally, yy_zwally, zwally_basins = get_zwally_basins("Data/DrainageBasins/ZwallyBasins.mat")
        println("   ✓ Success! Loaded $(length(zwally_basins)) basin points")
        println("   ✓ Basin ID range: $(minimum(zwally_basins)) to $(maximum(zwally_basins))")
        println("   ✓ X range: $(minimum(xx_zwally)) to $(maximum(xx_zwally)) m")
        println("   ✓ Y range: $(minimum(yy_zwally)) to $(maximum(yy_zwally)) m")
    catch e
        println("   ✗ Error: ", e)
    end
else
    println("   ⚠ Skipping: Data file not found (Data/DrainageBasins/ZwallyBasins.mat)")
end

# Test 3: get_frank_temps
println("\n3. Testing get_frank_temps()...")
try
    temp_data = get_frank_temps("Data/FranksTemps.mat")
    println("   ✓ Success! Loaded Frank temperature data")
    println("   ✓ FranksTemps shape: $(size(temp_data.FranksTemps))")
    println("   ✓ Number of levels: $(size(temp_data.FranksTemps, 1))")
    println("   ✓ Number of spatial points: $(size(temp_data.FranksTemps, 2))")
    println("   ✓ Temperature range: $(minimum(temp_data.FranksTemps)) to $(maximum(temp_data.FranksTemps)) K")
    println("   ✓ Sigma levels: $(length(temp_data.sigmaTemp))")
catch e
    println("   ✗ Error: ", e)
end

# Test 4: get_measures_mat
println("\n4. Testing get_measures_mat()...")
if isfile("Data/MEaSUREs/MEaSUREsAntVels.mat")
    try
        velocity_data = get_measures_mat("Data/MEaSUREs/MEaSUREsAntVels.mat")
        println("   ✓ Success! Loaded MEaSUREs velocity data")
        println("   ✓ xx_v shape: $(size(velocity_data.xx_v))")
        println("   ✓ yy_v shape: $(size(velocity_data.yy_v))")
        println("   ✓ vx shape: $(size(velocity_data.vx))")
        println("   ✓ vy shape: $(size(velocity_data.vy))")
        println("   ✓ Velocity magnitude range: $(minimum(sqrt.(velocity_data.vx.^2 .+ velocity_data.vy.^2))) to $(maximum(sqrt.(velocity_data.vx.^2 .+ velocity_data.vy.^2))) m/yr")
    catch e
        println("   ✗ Error: ", e)
    end
else
    println("   ⚠ Skipping: Data file not found (Data/MEaSUREs/MEaSUREsAntVels.mat)")
end

println("\n" * "=" ^ 60)
println("Testing complete!")
