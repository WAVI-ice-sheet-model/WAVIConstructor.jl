using WAVIConstructor.DataLoading
using NCDatasets
using Test

# Create a mock BISICLES NetCDF file for testing
mktempdir() do tmpdir
    mockfile = joinpath(tmpdir, "mock_bisicles.nc")
    nx, ny, nz = 4, 3, 2
    ds = NCDataset(mockfile, "c")
    defDim(ds, "x", nx)
    defDim(ds, "y", ny)
    defDim(ds, "z", nz)
    defDim(ds, "sigma", nz)

    vx = defVar(ds, "x", Float32, ("x",))
    vy = defVar(ds, "y", Float32, ("y",))
    vz = defVar(ds, "z", Float32, ("z",))
    vsigma = defVar(ds, "sigma", Float32, ("sigma",))
    # BISICLES temperature is stored as (z, y, x)
    vtemp = defVar(ds, "T", Float32, ("z", "y", "x"))

    vx[:] .= collect(1.0f0:nx)
    vy[:] .= collect(1.0f0:ny)
    vz[:] .= collect(1.0f0:nz)
    vsigma[:] .= collect(0.1f0:0.1f0:0.2f0)
    vtemp[:, :, :] .= reshape(collect(1.0f0:(nx*ny*nz)), nz, ny, nx)

    close(ds)

    println("Testing get_bisicles_temps on: ", mockfile)
    sigma, x, y, z, temps = get_bisicles_temps(mockfile, scale_xy=1000)

    @test size(sigma) == (nz,)
    @test size(x) == (nx,)
    @test size(y) == (ny,)
    @test size(z) == (nz,)
    # The code uses temps[i, :, :], so it expects (nz, ny, nx)
    @test size(temps) == (nz, ny, nx)

    @test all(x .== 1000 .* collect(1.0f0:nx))
    @test all(y .== 1000 .* collect(1.0f0:ny))
    @test all(z .== collect(1.0f0:nz))
    @test all(sigma .== collect(0.1f0:0.1f0:0.2f0))
    # Check that temps contains the expected values
    expected_temps = reshape(collect(1.0f0:(nx*ny*nz)), nz, ny, nx)
    @test all(temps .== expected_temps)
end
