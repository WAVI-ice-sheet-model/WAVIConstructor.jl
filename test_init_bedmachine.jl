#!/usr/bin/env julia
"""
Test script for init_bedmachine_v3 function
"""

using WAVIConstructor

println("Testing init_bedmachine function...")
println("=" ^ 60)

# Create a minimal params NamedTuple
# The function uses get() with defaults, so we can use an empty tuple or minimal params
params = (
    start_data = "BEDMACHINEV3",
    sub_samp = 8,
    veloc_data = "Measures_2014/15",
    temps = "BISICLES_8km"
)

println("\nParameters:")
println("  start_data: $(params.start_data)")
println("  sub_samp: $(params.sub_samp)")
println("  veloc_data: $(params.veloc_data)")
println("  temps: $(params.temps)")

println("\nInitializing BedMachine v3 grids...")
try
    Gh, Gu, Gv, Gc = init_bedmachine(params)
    
    println("✓ Success! Grids initialized")
    println("\nGrid Information:")
    println("  H-grid (Gh):")
    println("    Size: $(Gh.nx) x $(Gh.ny)")
    println("    dx: $(Gh.dx) m, dy: $(Gh.dy) m")
    println("    x0: $(Gh.x0) m, y0: $(Gh.y0) m")
    println("    Bed elevation range: $(minimum(Gh.b)) to $(maximum(Gh.b)) m")
    println("    Ice thickness range: $(minimum(Gh.h)) to $(maximum(Gh.h)) m")
    println("    Surface elevation range: $(minimum(Gh.s)) to $(maximum(Gh.s)) m")
    
    if haskey(Gh, :a_Arthern)
        println("    Arthern accumulation loaded: $(!isempty(Gh.a_Arthern))")
    end
    if haskey(Gh, :basin_id)
        println("    Zwally basins loaded: $(!isempty(Gh.basin_id))")
    end
    if haskey(Gh, :Tma)
        println("    ALBMAP Tma loaded: $(!isempty(Gh.Tma))")
    end
    
    println("\n  U-grid (Gu):")
    println("    Size: $(Gu.nx) x $(Gu.ny)")
    if haskey(Gu, :u_data)
        println("    Velocity data loaded: $(!isempty(Gu.u_data))")
        println("    Velocity range: $(minimum(filter(!isnan, Gu.u_data))) to $(maximum(filter(!isnan, Gu.u_data))) m/yr")
    end
    
    println("\n  V-grid (Gv):")
    println("    Size: $(Gv.nx) x $(Gv.ny)")
    if haskey(Gv, :v_data)
        println("    Velocity data loaded: $(!isempty(Gv.v_data))")
        println("    Velocity range: $(minimum(filter(!isnan, Gv.v_data))) to $(maximum(filter(!isnan, Gv.v_data))) m/yr")
    end
    
    println("\n  C-grid (Gc):")
    println("    Size: $(Gc.nx) x $(Gc.ny)")
    
    if haskey(Gh, :levels)
        println("\n  Temperature data:")
        println("    Number of levels: $(length(Gh.levels.sigmas))")
        println("    Temperature shape: $(size(Gh.levels.temperature))")
    end
    
    println("\n" * "=" ^ 60)
    println("Test completed successfully!")
    
catch e
    println("✗ Error occurred:")
    println("  ", e)
    if isa(e, LoadError)
        println("  Original error: ", e.error)
    end
    rethrow(e)
end
