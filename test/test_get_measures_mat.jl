using WAVIConstructor.DataLoading
using MAT
using Test

# Create a mock MEaSUREs .mat file for testing
mktempdir() do tmpdir
    mockfile = joinpath(tmpdir, "mock_MEaSUREsAntVels.mat")
    
    # Create test data
    nx, ny = 5, 4
    xx_v = repeat(collect(-3000000.0:100000.0:-2600000.0), 1, ny)
    yy_v = repeat(collect(-3000000.0:100000.0:-2700000.0)', nx, 1)
    vx = reshape(collect(10.0:1.0:(10.0 + nx*ny - 1)), nx, ny)
    vy = reshape(collect(20.0:1.0:(20.0 + nx*ny - 1)), nx, ny)
    
    # Write to .mat file
    matopen(mockfile, "w") do mat_file
        write(mat_file, "xx_v", xx_v)
        write(mat_file, "yy_v", yy_v)
        write(mat_file, "vx", vx)
        write(mat_file, "vy", vy)
    end

    println("Testing get_measures_mat on: ", mockfile)
    velocity_data = get_measures_mat(mockfile)

    @test size(velocity_data.xx_v) == (nx, ny)
    @test size(velocity_data.yy_v) == (nx, ny)
    @test size(velocity_data.vx) == (nx, ny)
    @test size(velocity_data.vy) == (nx, ny)

    # Check that data matches
    @test all(velocity_data.xx_v .== xx_v)
    @test all(velocity_data.yy_v .== yy_v)
    @test all(velocity_data.vx .== vx)
    @test all(velocity_data.vy .== vy)

    # Check velocity values are reasonable
    @test minimum(velocity_data.vx) >= 10.0
    @test maximum(velocity_data.vx) <= 30.0
    @test minimum(velocity_data.vy) >= 20.0
    @test maximum(velocity_data.vy) <= 40.0
end
