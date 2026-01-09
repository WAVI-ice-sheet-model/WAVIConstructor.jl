using WAVIConstructor.DataLoading
using MAT
using Test

# Create a mock Zwally basins .mat file for testing
mktempdir() do tmpdir
    mockfile = joinpath(tmpdir, "mock_ZwallyBasins.mat")
    
    # Create test data
    nx, ny = 5, 4
    xx_zwally_full = repeat(collect(-3000000.0:100000.0:-2600000.0), 1, ny)
    yy_zwally_full = repeat(collect(-3000000.0:100000.0:-2700000.0)', nx, 1)
    zwally_basins_full = reshape(collect(0:1:(nx*ny-1)), nx, ny)  # Some zeros, some positive IDs
    
    # Write to .mat file
    matopen(mockfile, "w") do mat_file
        write(mat_file, "xxZwallyBasins", xx_zwally_full)
        write(mat_file, "yyZwallyBasins", yy_zwally_full)
        write(mat_file, "ZwallyBasins", zwally_basins_full)
    end

    println("Testing get_zwally_basins on: ", mockfile)
    xx_zwally, yy_zwally, zwally_basins = get_zwally_basins(mockfile)

    # Should have filtered out zeros (1 out of 20)
    expected_count = nx * ny - 1  # One zero value
    @test length(xx_zwally) == expected_count
    @test length(yy_zwally) == expected_count
    @test length(zwally_basins) == expected_count

    # Check that all arrays have the same length
    @test length(xx_zwally) == length(yy_zwally) == length(zwally_basins)

    # Check that all basin IDs are > 0
    @test all(zwally_basins .> 0)

    # Check that values are reasonable
    @test minimum(zwally_basins) > 0
    @test maximum(zwally_basins) <= (nx * ny - 1)
end
