using WAVIConstructor.DataLoading
using WAVIConstructor.DataSources
using NCDatasets
using ArchGDAL
using MAT
using Test


# Helper function to create a GeoTIFF file
function create_mock_geotiff(tmpdir, filename, width, height, dx, dy, mapx, mapy, data)
    filepath = joinpath(tmpdir, filename)
    driver = ArchGDAL.getdriver("GTiff")
    ArchGDAL.create(filepath; driver=driver, width=width, height=height, nbands=1, dtype=Float32) do dataset
        ArchGDAL.setgeotransform!(dataset, [mapx, dx, 0.0, mapy, 0.0, dy])
        ArchGDAL.write!(ArchGDAL.getband(dataset, 1), data)
    end
    return filepath
end

# Helper function to create a MATLAB .mat file
function create_mock_matfile(tmpdir, filename, vars)
    filepath = joinpath(tmpdir, filename)
    matopen(filepath, "w") do mat_file
        for (name, data) in vars
            write(mat_file, name, data)
        end
    end
    return filepath
end

@testset "Data Loading Functions" begin
    # Test load_data(ALBMAPv1(), ...)
    mktempdir() do tmpdir
        mockfile = joinpath(tmpdir, "mock_albmap.nc")
        nx, ny = 4, 3
        ds = NCDataset(mockfile, "c")
        defDim(ds, "x1", nx)
        defDim(ds, "y1", ny)
        vx = defVar(ds, "x1", Float32, ("x1",))
        vy = defVar(ds, "y1", Float32, ("y1",))
        vx[:] .= collect(1.0f0:nx)
        vy[:] .= collect(1.0f0:ny)
        for var in ["usrf", "lsrf", "topg", "acca", "accr", "firn", "temp", "mask", "mask_plus"]
            v = defVar(ds, var, Float32, ("x1", "y1"))
            v[:, :] .= ones(Float32, nx, ny)
        end
        close(ds)

        Gh = load_data(ALBMAPv1(), mockfile)
        @test size(Gh.h) == (1, 12)
        @test all(Gh.h .== -1.0f0)
    end

    # Test load_data(ArthernAccumulation(), ...)
    mktempdir() do tmpdir
        mockfile = joinpath(tmpdir, "mock_amsr_accumulation_map.txt")
        
        open(mockfile, "w") do f
            # Write 21 header lines (realistic header format)
            println(f, "# Arthern, R. J., D. P. Winebrenner, and D. G. Vaughan (2006)")
            println(f, "# Antarctic snow accumulation mapped using polarization of 4.3-cm wavelength microwave emission")
            println(f, "# J. Geophys. Res., 111, D06107, doi:10.1029/2004JD005667")
            println(f, "#")
            println(f, "# Effective resolution approximately 100 km")
            println(f, "#")
            println(f, "# WARNING: values for locations subject to melt may be unreliable")
            println(f, "#")
            println(f, "# Columns:")
            println(f, "# (1) Latitude (degrees)")
            println(f, "# (2) Longitude (degrees, 0-360 convention)")
            println(f, "# (3) X(m) - Polar Stereographic (latitude of true scale 71S, ellipsoid WGS84)")
            println(f, "# (4) Y(m) - Polar Stereographic (latitude of true scale 71S, ellipsoid WGS84)")
            println(f, "# (5) Accumulation Rate (kg/m^2/a water equivalent)")
            println(f, "# (6) Estimated RMS error at 100 km scale (% of column 5)")
            println(f, "#")
            for i in 17:21
                println(f, "# Header line $i")
            end
            
            # Write realistic data with Antarctic coordinates
            # Realistic accumulation values: 100-1500 kg/m²/a (water equivalent)
            # Realistic error: 15-25%
            # Longitude in 0-360 convention
            # Scatter NaN values throughout (not just at end)
            test_data = [
                # lat, lon (0-360), x (polar stereo), y (polar stereo), acc, err
                (-65.0, 120.0, -1200000.0, 1500000.0, NaN, NaN),       # NaN - ocean
                (-70.5, 150.0, -800000.0, 2000000.0, 450.0, 20.0),     # Valid
                (-75.0, 180.0, 0.0, 2200000.0, 180.0, 18.0),           # Valid - low acc
                (-80.0, 200.0, 500000.0, 1800000.0, NaN, NaN),         # NaN - gap
                (-72.0, 250.0, 1200000.0, 800000.0, 850.0, 22.0),      # Valid
                (-78.0, 280.0, 1500000.0, -200000.0, 1200.0, 19.0),    # Valid - high acc
                (-82.0, 300.0, 1000000.0, -800000.0, 350.0, 24.0),     # Valid
                (-85.0, 0.0, -500000.0, -1500000.0, NaN, NaN),         # NaN - interior
                (-68.0, 45.0, -2000000.0, 500000.0, 680.0, 21.0),      # Valid
                (-73.5, 90.0, -1500000.0, 1200000.0, 520.0, 17.0),     # Valid
                (-76.0, 330.0, -800000.0, -500000.0, 1050.0, 23.0),    # Valid
                (-88.0, 180.0, 100000.0, -100000.0, 50.0, 25.0),       # Valid - very low (interior)
            ]
            
            for (lat, lon, x, y, acc, err) in test_data
                println(f, "$lat $lon $x $y $acc $err")
            end
        end

        acc_data = load_data(ArthernAccumulation(), mockfile)

        # Should have 9 valid rows (3 NaN rows filtered out from 12 total)
        @test length(acc_data.x) == 9
        @test length(acc_data.y) == 9
        @test length(acc_data.acc) == 9
        @test length(acc_data.x) == length(acc_data.y) == length(acc_data.acc)
        
        # No NaN values in filtered output
        @test !any(isnan.(acc_data.acc))
        
        # Accumulation converted from water equivalent to ice equivalent (divided by 917)
        # Original values: 450, 180, 850, 1200, 350, 680, 520, 1050, 50
        # Min ice equivalent: 50/917 ≈ 0.055
        # Max ice equivalent: 1200/917 ≈ 1.31
        @test all(acc_data.acc .> 0)
        @test minimum(acc_data.acc) > 0.05
        @test maximum(acc_data.acc) < 1.5
        
        # Verify coordinate ranges are realistic (Antarctic polar stereographic)
        @test all(abs.(acc_data.x) .< 3000000)  # Reasonable polar stereographic range
        @test all(abs.(acc_data.y) .< 3000000)
    end

    # Test load_data(BISICLESTemps(), ...)
    mktempdir() do tmpdir
        nx, ny, nz = 4, 3, 2
        mockfile = joinpath(tmpdir, "mock_bisicles.nc")
        ds = NCDataset(mockfile, "c")
        defDim(ds, "x", nx)
        defDim(ds, "y", ny)
        defDim(ds, "z", nz)
        defDim(ds, "sigma", nz)

        vx = defVar(ds, "x", Float32, ("x",))
        vy = defVar(ds, "y", Float32, ("y",))
        vz = defVar(ds, "z", Float32, ("z",))
        vsigma = defVar(ds, "sigma", Float32, ("sigma",))
        vtemp = defVar(ds, "T", Float32, ("z", "y", "x"))

        vx[:] .= collect(1.0f0:nx)
        vy[:] .= collect(1.0f0:ny)
        vz[:] .= collect(1.0f0:nz)
        vsigma[:] .= collect(0.1f0:0.1f0:0.2f0)
        vtemp[:, :, :] .= reshape(collect(1.0f0:(nx*ny*nz)), nz, ny, nx)
        close(ds)

        data = load_data(BISICLESTemps(), mockfile, scale_xy=1000)

        @test size(data.sigma) == (nz,)
        @test size(data.xx) == (nx,)
        @test size(data.yy) == (ny,)
        @test size(data.temps) == (nz, ny, nx)
        @test all(data.xx .== 1000 .* collect(1.0f0:nx))
        @test all(data.yy .== 1000 .* collect(1.0f0:ny))
        @test all(data.sigma .== collect(0.1f0:0.1f0:0.2f0))
        expected_temps = reshape(collect(1.0f0:(nx*ny*nz)), nz, ny, nx)
        @test all(data.temps .== expected_temps)
    end

    # Test load_data(FrankTemps(), ...)
    mktempdir() do tmpdir
        nlevels = 3
        npoints = 10
        FranksTemps = reshape(collect(200.0:1.0:(200.0 + nlevels*npoints - 1)), nlevels, npoints)
        xxTemp = collect(-3000000.0:100000.0:(-3000000.0 + (npoints-1)*100000.0))
        yyTemp = collect(-3000000.0:100000.0:(-3000000.0 + (npoints-1)*100000.0))
        sigmaTemp = collect(0.1:0.1:0.3)
        
        mockfile = create_mock_matfile(tmpdir, "mock_FranksTemps.mat",
            [("FranksTemps", FranksTemps), ("xxTemp", xxTemp), ("yyTemp", yyTemp), ("sigmaTemp", sigmaTemp)]
        )

        temp_data = load_data(FrankTemps(), mockfile)

        @test size(temp_data.temps) == (nlevels, npoints)
        @test length(temp_data.xx) == npoints
        @test length(temp_data.yy) == npoints
        @test length(temp_data.sigma) == nlevels
        @test all(temp_data.temps .== FranksTemps)
        @test all(temp_data.xx .== xxTemp)
        @test all(temp_data.yy .== yyTemp)
        @test all(temp_data.sigma .== sigmaTemp)
        @test minimum(temp_data.temps) >= 200.0
        @test maximum(temp_data.temps) <= 300.0
    end

    # Test load_data(MEaSUREs(), ...)
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

        data = load_data(MEaSUREs(), mockfile)

        @test size(data.xx) == (nx, ny)
        @test size(data.yy) == (nx, ny)
        @test size(data.vx) == (nx, ny)
        @test size(data.vy) == (nx, ny)
        @test all(data.xx .== repeat(collect(1.0f0:nx), 1, ny))
        @test all(data.yy .== repeat(collect(1.0f0:ny)', nx, 1))
        @test all(data.vx .== reshape(collect(1.0f0:(nx*ny)), nx, ny))
        @test all(data.vy .== reshape(collect(1.0f0:(nx*ny)), nx, ny))
    end

    # Test load_data(SmithDhdt(), ...)
    mktempdir() do tmpdir
        width, height = 4, 3
        dx, dy = 10.0, -10.0
        mapx, mapy = 100.0, 200.0
        mock_dhdt = Float32[1.5 2.5 3.5 4.5; -1.0 0.0 1.0 2.0; -2.5 -1.5 -0.5 0.5]

        grnd_file = create_mock_geotiff(tmpdir, "mock_smith_grnd.tif", width, height, dx, dy, mapx, mapy,
            permutedims(mock_dhdt, (2, 1)))
        flt_file = create_mock_geotiff(tmpdir, "mock_smith_flt.tif", width, height, dx, dy, mapx, mapy,
            permutedims(mock_dhdt .* 0.5, (2, 1)))

        result = load_data(SmithDhdt(), "", grnd_file=grnd_file, flt_file=flt_file)
        
        # The function flips data with reverse(dhdt, dims=1), so expected data should be flipped
        expected_grnd = reverse(mock_dhdt, dims=1)
        expected_flt = reverse(mock_dhdt .* 0.5, dims=1)
        
        # Test grounded data
        @test size(result.grnd_xx) == (height, width)
        @test size(result.grnd_yy) == (height, width) 
        @test size(result.grnd_dhdt) == (height, width)
        expected_x = [mapx + dx * (i - 0.5) for i in 1:width]
        expected_y = [mapy - abs(dy) * (i - 0.5) for i in 1:height]
        @test all(result.grnd_xx[1, :] .≈ expected_x)
        @test all(result.grnd_yy[:, 1] .≈ reverse(expected_y))
        @test all(result.grnd_dhdt .≈ expected_grnd)
        
        # Test floating data
        @test size(result.flt_xx) == (height, width)
        @test size(result.flt_yy) == (height, width)
        @test size(result.flt_dhdt) == (height, width)
        @test all(result.flt_dhdt .≈ expected_flt)
    end

    # Test load_data(ZwallyBasins(), ...)
    mktempdir() do tmpdir
        # Create a more realistic mock: 20x20 grid with realistic Antarctic coordinates
        # and realistic basin IDs (1-27 range) with contiguous regions
        nx, ny = 20, 20
        
        # Antarctic polar stereographic coordinates (approx range)
        x_coords = collect(range(-2500000.0, 500000.0, length=nx))
        y_coords = collect(range(-2500000.0, 500000.0, length=ny))
        
        xx_zwally_full = repeat(x_coords, 1, ny)
        yy_zwally_full = repeat(y_coords', nx, 1)
        
        # Create basin IDs with contiguous regions (realistic spatial clustering)
        # Basin IDs should be 0 (ocean/ice-free) or 1-27
        zwally_basins_full = zeros(Int, nx, ny)
        
        # Create 9 basin regions (3x3 quadrants) with basin IDs 1-9
        # Add some zeros around edges (ocean/ice-free areas)
        for i in 1:nx
            for j in 1:ny
                # Set edges to 0 (ocean)
                if i <= 2 || i >= nx-1 || j <= 2 || j >= ny-1
                    zwally_basins_full[i, j] = 0
                else
                    # Divide interior into 9 regions with basin IDs 1-9
                    region_i = div(i - 3, 5) + 1
                    region_j = div(j - 3, 5) + 1
                    zwally_basins_full[i, j] = min((region_i - 1) * 3 + region_j, 27)
                end
            end
        end
        
        # Convert to Float64 to match actual file format
        xx_zwally_full = Float64.(xx_zwally_full)
        yy_zwally_full = Float64.(yy_zwally_full)
        zwally_basins_full = Float64.(zwally_basins_full)
        
        mockfile = create_mock_matfile(tmpdir, "mock_ZwallyBasins.mat",
            [("xxZwallyBasins", xx_zwally_full), ("yyZwallyBasins", yy_zwally_full), ("ZwallyBasins", zwally_basins_full)]
        )

        data = load_data(ZwallyBasins(), mockfile)

        # Count expected valid points (basin > 0)
        expected_count = count(x -> x > 0, zwally_basins_full)
        
        @test length(data.xx) == expected_count
        @test length(data.yy) == expected_count
        @test length(data.basins) == expected_count
        @test length(data.xx) == length(data.yy) == length(data.basins)
        @test all(data.basins .> 0)
        @test minimum(data.basins) >= 1
        @test maximum(data.basins) <= 27
        
        # Verify coordinate ranges are realistic (Antarctic extent)
        @test minimum(data.xx) >= -3000000.0
        @test maximum(data.xx) <= 3000000.0
        @test minimum(data.yy) >= -3000000.0
        @test maximum(data.yy) <= 3000000.0
    end

    # Test geotiff_read_axis_only
    mktempdir() do tmpdir
        width, height = 4, 3
        dx, dy = 10.0, -10.0
        mapx, mapy = 100.0, 200.0

        mockfile = create_mock_geotiff(tmpdir, "mock_geotiff.tif", width, height, dx, dy, mapx, mapy,
            ones(Float32, width, height))

        result = geotiff_read_axis_only(mockfile)

        @test length(result.x) == width
        @test length(result.y) == height
        @test result.x[1] == mapx
        @test result.x[end] == mapx + (width - 1) * dx
        @test result.y[1] == mapy
        @test result.y[end] == mapy - (height - 1) * abs(dy)
        @test result.info.map_info.dx == dx
        @test result.info.map_info.dy == abs(dy)
        @test result.info.map_info.mapx == mapx
        @test result.info.map_info.mapy == mapy
    end
end

# Integration tests with real data files (skipped if files not present)
@testset "Integration: Real Data Files" begin
    
    # Test loading actual Arthern accumulation file
    arthern_file = "Data/amsr_accumulation_map.txt"
    if isfile(arthern_file)
        @testset "Real Arthern Accumulation File" begin
            acc_data = load_data(ArthernAccumulation(), arthern_file)
            
            # File should contain valid data points
            @test length(acc_data.x) > 0
            @test length(acc_data.x) == length(acc_data.y) == length(acc_data.acc)
            
            # No NaN values in filtered output
            @test !any(isnan.(acc_data.acc))
            @test !any(isnan.(acc_data.x))
            
            # X and Y should be in reasonable polar stereographic range (within ~3000 km of pole)
            @test all(abs.(acc_data.x) .< 3500000)
            @test all(abs.(acc_data.y) .< 3500000)
            
            # Accumulation (ice equivalent) should be positive and reasonable
            # Original water equiv is ~50-2000 kg/m²/a, so ice equiv is ~0.05-2.2 m/a
            @test all(acc_data.acc .> 0)
            @test minimum(acc_data.acc) > 0.01  # At least 10 mm/a
            @test maximum(acc_data.acc) < 5.0   # Less than 5 m/a
            
            @info "Loaded $(length(acc_data.acc)) valid accumulation data points from real file"
        end
    else
        @warn "Skipping Arthern accumulation integration test: file not found at $arthern_file"
    end
    
    # Test loading actual Zwally basins file
    zwally_file = "Data/ZwallyBasins.mat"
    if isfile(zwally_file)
        @testset "Real Zwally Basins File" begin
            data = load_data(ZwallyBasins(), zwally_file)
            
            # File should contain valid data points
            @test length(data.xx) > 0
            @test length(data.xx) == length(data.yy) == length(data.basins)
            
            # All basin IDs should be positive (filtered out zeros)
            @test all(data.basins .> 0)
            
            # Basin IDs should be in valid range (1-27 for Zwally drainage basins)
            @test minimum(data.basins) >= 1
            @test maximum(data.basins) <= 27
            
            # Coordinates should be in reasonable polar stereographic range
            @test all(abs.(data.xx) .< 4000000)
            @test all(abs.(data.yy) .< 4000000)
            
            # Check that multiple basins are represented
            unique_basins = unique(data.basins)
            @test length(unique_basins) > 1  # Should have more than 1 basin
            
            @info "Loaded $(length(data.basins)) valid basin data points from real file"
            @info "Unique basin IDs: $(sort(unique(data.basins)))"
        end
    else
        @warn "Skipping Zwally basins integration test: file not found at $zwally_file"
    end
end

