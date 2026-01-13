using Test
using WAVIConstructor.DomainSelection

@testset "Domain Selection" begin
    # Create mock grid structures
    nx, ny = 10, 10
    
    # Create H-grid
    Gh = (
        nx = nx,
        ny = ny,
        ok = trues(nx, ny),  # All ice points
        rock = falses(nx, ny),
        basinID = reshape(collect(1:100), nx, ny),
        mask = trues(nx, ny)
    )
    
    # Create U-grid
    Gu = (
        nx = nx + 1,
        ny = ny,
        mask = trues(nx + 1, ny)
    )
    
    # Create V-grid
    Gv = (
        nx = nx,
        ny = ny + 1,
        mask = trues(nx, ny + 1)
    )
    
    # Create C-grid
    Gc = (
        nx = nx - 1,
        ny = ny - 1,
        mask = trues(nx - 1, ny - 1)
    )
    
    # Test with single basin
    params = (basins = [1],)
    Gh_test, Gu_test, Gv_test, Gc_test = select_domain_wavi(Gh, Gu, Gv, Gc, params)
    
    @test haskey(Gh_test, :globalMask)
    @test haskey(Gh_test, :basinOK)
    @test haskey(Gh_test, :mask)
    @test haskey(Gh_test, :f)
    @test haskey(Gh_test, :n)
    @test haskey(Gh_test, :otherBasins)
    
    @test haskey(Gu_test, :rockEdges)
    @test haskey(Gu_test, :mask)
    @test haskey(Gu_test, :f)
    @test haskey(Gu_test, :n)
    @test haskey(Gu_test, :halo)
    
    @test haskey(Gv_test, :rockEdges)
    @test haskey(Gv_test, :mask)
    @test haskey(Gv_test, :f)
    @test haskey(Gv_test, :n)
    @test haskey(Gv_test, :halo)
    
    @test haskey(Gc_test, :globalMask)
    @test haskey(Gc_test, :mask)
    @test haskey(Gc_test, :f)
    @test haskey(Gc_test, :n)
    
    # Test that masks are boolean arrays
    @test eltype(Gh_test.globalMask) == Bool
    @test eltype(Gh_test.mask) == Bool
    @test eltype(Gu_test.mask) == Bool
    @test eltype(Gv_test.mask) == Bool
    @test eltype(Gc_test.mask) == Bool
    
    # Test that indices are vectors
    @test Gh_test.f isa Vector
    @test Gu_test.f isa Vector
    @test Gv_test.f isa Vector
    @test Gc_test.f isa Vector
    
    # Test that counts are integers
    @test Gh_test.n isa Integer
    @test Gu_test.n isa Integer
    @test Gv_test.n isa Integer
    @test Gc_test.n isa Integer
end
