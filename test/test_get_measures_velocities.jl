using WAVIConstructor.DataLoading
using NCDatasets
using Test

# Create a mock NetCDF file for testing
mktempdir() do tmpdir
    mockfile = joinpath(tmpdir, "mock_measures.nc")
    nx, ny = 4, 3
    ds = NCDataset(mockfile, "c")
    defDim(ds, "x", nx)
    defDim(ds, "y", ny)

    vx = defVar(ds, "x", Float32, ("x",))
    vy = defVar(ds, "y", Float32, ("y",))
    vVX = defVar(ds, "VX", Float32, ("x", "y"))
    vVY = defVar(ds, "VY", Float32, ("x", "y"))

    vx[:] .= collect(1.0f0:nx)
    vy[:] .= collect(1.0f0:ny)
    vVX[:, :] .= reshape(collect(1.0f0:(nx*ny)), nx, ny)
    vVY[:, :] .= reshape(collect(1.0f0:(nx*ny)), nx, ny)

    close(ds)

    println("Testing get_measures_velocities on: ", mockfile)
    xx, yy, vx, vy = get_measures_velocities(mockfile)

    @test size(xx) == (nx, ny)
    @test size(yy) == (nx, ny)
    @test size(vx) == (nx, ny)
    @test size(vy) == (nx, ny)

    @test all(xx .== repeat(collect(1.0f0:nx), 1, ny))
    @test all(yy .== repeat(collect(1.0f0:ny)', nx, 1))
    @test all(vx .== reshape(collect(1.0f0:(nx*ny)), nx, ny))
    @test all(vy .== reshape(collect(1.0f0:(nx*ny)), nx, ny))
end
