using WAVIConstructor.DataLoading
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
    # Test get_albmap
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

        Gh = get_albmap(mockfile)
        @test size(Gh[:h]) == (1, 12)
        @test all(Gh[:h] .== -1.0f0)
    end

    # Test get_arthern_accumulation
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

        aa_lat, aa_lon, aa_x, aa_y, aa_acc, aa_err = get_arthern_accumulation(mockfile)

        # Should have 9 valid rows (3 NaN rows filtered out from 12 total)
        @test length(aa_lat) == 9
        @test length(aa_lon) == 9
        @test length(aa_x) == 9
        @test length(aa_y) == 9
        @test length(aa_acc) == 9
        @test length(aa_err) == 9
        @test length(aa_lat) == length(aa_lon) == length(aa_x) == length(aa_y) == length(aa_acc) == length(aa_err)
        
        # No NaN values in filtered output
        @test !any(isnan.(aa_acc))
        @test !any(isnan.(aa_err))
        
        # Accumulation converted from water equivalent to ice equivalent (divided by 917)
        # Original values: 450, 180, 850, 1200, 350, 680, 520, 1050, 50
        # Min ice equivalent: 50/917 ≈ 0.055
        # Max ice equivalent: 1200/917 ≈ 1.31
        @test all(aa_acc .> 0)
        @test minimum(aa_acc) > 0.05
        @test maximum(aa_acc) < 1.5
        
        # Error values should be realistic percentages (15-25%)
        @test minimum(aa_err) >= 15
        @test maximum(aa_err) <= 25
        
        # Verify coordinate ranges are realistic (Antarctic polar stereographic)
        @test all(aa_lat .< -60)  # All Antarctic latitudes
        @test all(aa_lon .>= 0) && all(aa_lon .<= 360)  # 0-360 longitude convention
        @test all(abs.(aa_x) .< 3000000)  # Reasonable polar stereographic range
        @test all(abs.(aa_y) .< 3000000)
    end

    # Test get_bisicles_temps
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

        sigma, x, y, z, temps = get_bisicles_temps(mockfile, scale_xy=1000)

        @test size(sigma) == (nz,)
        @test size(x) == (nx,)
        @test size(y) == (ny,)
        @test size(z) == (nz,)
        @test size(temps) == (nz, ny, nx)
        @test all(x .== 1000 .* collect(1.0f0:nx))
        @test all(y .== 1000 .* collect(1.0f0:ny))
        @test all(z .== collect(1.0f0:nz))
        @test all(sigma .== collect(0.1f0:0.1f0:0.2f0))
        expected_temps = reshape(collect(1.0f0:(nx*ny*nz)), nz, ny, nx)
        @test all(temps .== expected_temps)
    end

    # Test get_frank_temps
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

        temp_data = get_frank_temps(mockfile)

        @test size(temp_data.FranksTemps) == (nlevels, npoints)
        @test length(temp_data.xxTemp) == npoints
        @test length(temp_data.yyTemp) == npoints
        @test length(temp_data.sigmaTemp) == nlevels
        @test all(temp_data.FranksTemps .== FranksTemps)
        @test all(temp_data.xxTemp .== xxTemp)
        @test all(temp_data.yyTemp .== yyTemp)
        @test all(temp_data.sigmaTemp .== sigmaTemp)
        @test minimum(temp_data.FranksTemps) >= 200.0
        @test maximum(temp_data.FranksTemps) <= 300.0
    end

    # Test get_measures_velocities
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

    # Test get_smith_dhdt
    mktempdir() do tmpdir
        width, height = 4, 3
        dx, dy = 10.0, -10.0
        mapx, mapy = 100.0, 200.0
        mock_dhdt = Float32[1.5 2.5 3.5 4.5; -1.0 0.0 1.0 2.0; -2.5 -1.5 -0.5 0.5]

        grnd_file = create_mock_geotiff(tmpdir, "mock_smith_grnd.tif", width, height, dx, dy, mapx, mapy,
            permutedims(mock_dhdt, (2, 1)))
        flt_file = create_mock_geotiff(tmpdir, "mock_smith_flt.tif", width, height, dx, dy, mapx, mapy,
            permutedims(mock_dhdt .* 0.5, (2, 1)))

        result = get_smith_dhdt(grnd_file=grnd_file, flt_file=flt_file)
        
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

    # Test get_zwally_basins
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

        xx_zwally, yy_zwally, zwally_basins = get_zwally_basins(mockfile)

        # Count expected valid points (basin > 0)
        expected_count = count(x -> x > 0, zwally_basins_full)
        
        @test length(xx_zwally) == expected_count
        @test length(yy_zwally) == expected_count
        @test length(zwally_basins) == expected_count
        @test length(xx_zwally) == length(yy_zwally) == length(zwally_basins)
        @test all(zwally_basins .> 0)
        @test minimum(zwally_basins) >= 1
        @test maximum(zwally_basins) <= 27
        
        # Verify coordinate ranges are realistic (Antarctic extent)
        @test minimum(xx_zwally) >= -3000000.0
        @test maximum(xx_zwally) <= 3000000.0
        @test minimum(yy_zwally) >= -3000000.0
        @test maximum(yy_zwally) <= 3000000.0
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
            aa_lat, aa_lon, aa_x, aa_y, aa_acc, aa_err = get_arthern_accumulation(arthern_file)
            
            # File should contain valid data points
            @test length(aa_lat) > 0
            @test length(aa_lat) == length(aa_lon) == length(aa_x) == length(aa_y) == length(aa_acc) == length(aa_err)
            
            # No NaN values in filtered output
            @test !any(isnan.(aa_acc))
            @test !any(isnan.(aa_lat))
            @test !any(isnan.(aa_x))
            
            # Latitude should be Antarctic (south of 60°S)
            @test all(aa_lat .< -60)
            @test all(aa_lat .>= -90)
            
            # Longitude should be in 0-360 convention
            @test all(aa_lon .>= 0)
            @test all(aa_lon .<= 360)
            
            # X and Y should be in reasonable polar stereographic range (within ~3000 km of pole)
            @test all(abs.(aa_x) .< 3500000)
            @test all(abs.(aa_y) .< 3500000)
            
            # Accumulation (ice equivalent) should be positive and reasonable
            # Original water equiv is ~50-2000 kg/m²/a, so ice equiv is ~0.05-2.2 m/a
            @test all(aa_acc .> 0)
            @test minimum(aa_acc) > 0.01  # At least 10 mm/a
            @test maximum(aa_acc) < 5.0   # Less than 5 m/a
            
            # Error percentages should be reasonable (typically 5-30%)
            @test all(aa_err .> 0)
            @test minimum(aa_err) >= 5    # At least 5%
            @test maximum(aa_err) < 50    # Less than 50%
            
            @info "Loaded $(length(aa_acc)) valid accumulation data points from real file"
        end
    else
        @warn "Skipping Arthern accumulation integration test: file not found at $arthern_file"
    end
    
    # Test loading actual Zwally basins file
    zwally_file = "Data/ZwallyBasins.mat"
    if isfile(zwally_file)
        @testset "Real Zwally Basins File" begin
            xx_zwally, yy_zwally, zwally_basins = get_zwally_basins(zwally_file)
            
            # File should contain valid data points
            @test length(xx_zwally) > 0
            @test length(xx_zwally) == length(yy_zwally) == length(zwally_basins)
            
            # All basin IDs should be positive (filtered out zeros)
            @test all(zwally_basins .> 0)
            
            # Basin IDs should be in valid range (1-27 for Zwally drainage basins)
            @test minimum(zwally_basins) >= 1
            @test maximum(zwally_basins) <= 27
            
            # Coordinates should be in reasonable polar stereographic range
            @test all(abs.(xx_zwally) .< 4000000)
            @test all(abs.(yy_zwally) .< 4000000)
            
            # Check that multiple basins are represented
            unique_basins = unique(zwally_basins)
            @test length(unique_basins) > 1  # Should have more than 1 basin
            
            @info "Loaded $(length(zwally_basins)) valid basin data points from real file"
            @info "Unique basin IDs: $(sort(unique(zwally_basins)))"
        end
    else
        @warn "Skipping Zwally basins integration test: file not found at $zwally_file"
    end
end

