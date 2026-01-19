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
            for i in 1:21
                println(f, "# Header line $i")
            end
            
            for i in 1:10
                lat = -80.0 + i * 0.1
                lon = -100.0 + i * 0.1
                x = -3000000.0 + i * 1000.0
                y = -3000000.0 + i * 1000.0
                acc = i % 3 == 0 ? NaN : 0.1 * i * 917.0
                err = 0.01 * i
                println(f, "$lat $lon $x $y $acc $err")
            end
        end

        aa_lat, aa_lon, aa_x, aa_y, aa_acc, aa_err = get_arthern_accumulation(mockfile)

        @test length(aa_lat) == 7
        @test length(aa_lon) == 7
        @test length(aa_x) == 7
        @test length(aa_y) == 7
        @test length(aa_acc) == 7
        @test length(aa_err) == 7
        @test length(aa_lat) == length(aa_lon) == length(aa_x) == length(aa_y) == length(aa_acc) == length(aa_err)
        @test all(aa_acc .> 0)
        @test !any(isnan.(aa_acc))
        @test minimum(aa_acc) > 0
        @test maximum(aa_acc) < 10
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
        nx, ny = 5, 4
        xx_zwally_full = repeat(collect(-3000000.0:100000.0:-2600000.0), 1, ny)
        yy_zwally_full = repeat(collect(-3000000.0:100000.0:-2700000.0)', nx, 1)
        zwally_basins_full = reshape(collect(0:1:(nx*ny-1)), nx, ny)
        
        mockfile = create_mock_matfile(tmpdir, "mock_ZwallyBasins.mat",
            [("xxZwallyBasins", xx_zwally_full), ("yyZwallyBasins", yy_zwally_full), ("ZwallyBasins", zwally_basins_full)]
        )

        xx_zwally, yy_zwally, zwally_basins = get_zwally_basins(mockfile)

        expected_count = nx * ny - 1
        @test length(xx_zwally) == expected_count
        @test length(yy_zwally) == expected_count
        @test length(zwally_basins) == expected_count
        @test length(xx_zwally) == length(yy_zwally) == length(zwally_basins)
        @test all(zwally_basins .> 0)
        @test minimum(zwally_basins) > 0
        @test maximum(zwally_basins) <= (nx * ny - 1)
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

