using WAVIConstructor.DataLoading
using MAT
using Test

# Create a mock Arthern accumulation text file for testing
mktempdir() do tmpdir
    mockfile = joinpath(tmpdir, "mock_amsr_accumulation_map.txt")
    
    # Write 21 header lines
    open(mockfile, "w") do f
        for i in 1:21
            println(f, "# Header line $i")
        end
        
        # Write data rows: lat, lon, x, y, acc, err
        # Create some test data with a few NaN values
        for i in 1:10
            lat = -80.0 + i * 0.1
            lon = -100.0 + i * 0.1
            x = -3000000.0 + i * 1000.0
            y = -3000000.0 + i * 1000.0
            acc = i % 3 == 0 ? NaN : 0.1 * i * 917.0  # Some NaNs, others in water equivalent
            err = 0.01 * i
            println(f, "$lat $lon $x $y $acc $err")
        end
    end

    println("Testing get_arthern_accumulation on: ", mockfile)
    aa_lat, aa_lon, aa_x, aa_y, aa_acc, aa_err = get_arthern_accumulation(mockfile)

    # Should have filtered out NaN values (3 out of 10)
    @test length(aa_lat) == 7
    @test length(aa_lon) == 7
    @test length(aa_x) == 7
    @test length(aa_y) == 7
    @test length(aa_acc) == 7
    @test length(aa_err) == 7

    # Check that all arrays have the same length
    @test length(aa_lat) == length(aa_lon) == length(aa_x) == length(aa_y) == length(aa_acc) == length(aa_err)

    # Check that accumulation values are converted to ice equivalent (divided by 917)
    @test all(aa_acc .> 0)  # All should be positive after conversion
    @test !any(isnan.(aa_acc))  # No NaNs should remain

    # Check that values are reasonable
    @test minimum(aa_acc) > 0
    @test maximum(aa_acc) < 10  # Should be reasonable accumulation values
end
