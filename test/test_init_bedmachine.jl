using WAVIConstructor
using WAVIConstructor.InitBedMachine
using Test

# Test the grid creation helper functions
# Test create_h_grid
@testset "Grid creation functions" begin
    nx, ny = 10, 8
    x = collect(1000.0:500.0:(1000.0 + (nx-1)*500.0))
    y = collect(2000.0:500.0:(2000.0 + (ny-1)*500.0))
    bed = zeros(nx, ny)
    h = ones(nx, ny) .* 100.0
    s = ones(nx, ny) .* 200.0
    geoid = zeros(nx, ny)
    rockmask = zeros(nx, ny)
    mask = zeros(nx, ny)
    
    params = NamedTuple()
    Gh = InitBedMachine.create_h_grid(x, y, bed, h, s, geoid, rockmask, mask, params)
    
    @test Gh.nx == nx
    @test Gh.ny == ny
    @test Gh.dx == 500.0
    @test Gh.dy == 500.0
    @test size(Gh.xx) == (nx, ny)
    @test size(Gh.yy) == (nx, ny)
    @test size(Gh.b) == (nx, ny)
    @test size(Gh.h) == (nx, ny)
    @test size(Gh.s) == (nx, ny)
    @test haskey(Gh, :geoid)
    @test haskey(Gh, :mask)
    @test haskey(Gh, :rockmask)
    
    # Test U-grid
    Gu = InitBedMachine.create_u_grid(Gh)
    @test Gu.nx == Gh.nx + 1
    @test Gu.ny == Gh.ny
    @test size(Gu.xx) == (Gu.nx, Gu.ny)
    @test size(Gu.yy) == (Gu.nx, Gu.ny)
    
    # Test V-grid
    Gv = InitBedMachine.create_v_grid(Gh)
    @test Gv.nx == Gh.nx
    @test Gv.ny == Gh.ny + 1
    @test size(Gv.xx) == (Gv.nx, Gv.ny)
    @test size(Gv.yy) == (Gv.nx, Gv.ny)
    
    # Test C-grid
    Gc = InitBedMachine.create_c_grid(Gh)
    @test Gc.nx == Gh.nx - 1
    @test Gc.ny == Gh.ny - 1
    @test size(Gc.xx) == (Gc.nx, Gc.ny)
    @test size(Gc.yy) == (Gc.nx, Gc.ny)
    
    # Test coordinate consistency
    @test all(Gh.xx[:, 1] .== Gh.xx[:, 1])  # x-coords constant along y
    @test all(Gh.yy[1, :] .== Gh.yy[1, :])  # y-coords constant along x
end