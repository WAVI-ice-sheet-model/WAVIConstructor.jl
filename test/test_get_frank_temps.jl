using WAVIConstructor.DataLoading
using MAT
using Test

# Create a mock FranksTemps .mat file for testing
mktempdir() do tmpdir
    mockfile = joinpath(tmpdir, "mock_FranksTemps.mat")
    
    # Create test data
    nlevels = 3
    npoints = 10
    FranksTemps = reshape(collect(200.0:1.0:(200.0 + nlevels*npoints - 1)), nlevels, npoints)
    xxTemp = collect(-3000000.0:100000.0:(-3000000.0 + (npoints-1)*100000.0))
    yyTemp = collect(-3000000.0:100000.0:(-3000000.0 + (npoints-1)*100000.0))
    sigmaTemp = collect(0.1:0.1:0.3)
    
    # Write to .mat file
    matopen(mockfile, "w") do mat_file
        write(mat_file, "FranksTemps", FranksTemps)
        write(mat_file, "xxTemp", xxTemp)
        write(mat_file, "yyTemp", yyTemp)
        write(mat_file, "sigmaTemp", sigmaTemp)
    end

    println("Testing get_frank_temps on: ", mockfile)
    temp_data = get_frank_temps(mockfile)

    @test size(temp_data.FranksTemps) == (nlevels, npoints)
    @test length(temp_data.xxTemp) == npoints
    @test length(temp_data.yyTemp) == npoints
    @test length(temp_data.sigmaTemp) == nlevels

    # Check that data matches
    @test all(temp_data.FranksTemps .== FranksTemps)
    @test all(temp_data.xxTemp .== xxTemp)
    @test all(temp_data.yyTemp .== yyTemp)
    @test all(temp_data.sigmaTemp .== sigmaTemp)

    # Check temperature values are reasonable (in Kelvin)
    @test minimum(temp_data.FranksTemps) >= 200.0
    @test maximum(temp_data.FranksTemps) <= 300.0
end
