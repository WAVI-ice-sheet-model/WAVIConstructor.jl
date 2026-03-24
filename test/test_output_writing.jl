using Test
using NCDatasets

# We test the OutputWriting module through the top-level WAVIConstructor module
# since OutputWriting is an internal sub-module.

@testset "OutputWriting" begin

    # ── Build small synthetic grids ──────────────────────────────────
    nx, ny = 4, 5
    xx = [Float64(i) for i in 1:nx, j in 1:ny] .* 1000.0
    yy = [Float64(j) for i in 1:nx, j in 1:ny] .* 1000.0
    mask = trues(nx, ny)

    Gh = (
        nx_clip = nx,
        ny_clip = ny,
        dx      = 1000.0,
        x0_clip = 0.0,
        y0_clip = 0.0,
        xx_clip = xx,
        yy_clip = yy,
        h_clip  = rand(nx, ny),
        s_clip  = rand(nx, ny),
        b_clip  = rand(nx, ny),
        mask_clip          = mask,
        basinID_clip       = ones(nx, ny),
        a_Arthern_clip     = rand(nx, ny),
        dhdt_clip          = rand(nx, ny),
        dhdtAccDataMask_clip = mask,
        levels = (
            sigma_full       = [0.0, 0.25, 0.5, 0.75, 1.0],
            temperature_clip = rand(nx, ny, 5),
        ),
    )

    Gu = (
        uData_clip         = rand(nx, ny),
        uDataMaskFull_clip = mask,
        uiszero_clip       = falses(nx, ny),
    )

    Gv = (
        vData_clip         = rand(nx, ny),
        vDataMaskFull_clip = mask,
        viszero_clip       = falses(nx, ny),
    )

    # ── Test binary output ───────────────────────────────────────────
    @testset "write_binary_files" begin
        dir = mktempdir()
        WAVIConstructor.OutputWriting.write_binary_files(Gh, Gu, Gv, dir)

        expected_files = [
            "thickness.bin", "surface.bin", "bed.bin", "h_mask.bin",
            "basinID.bin", "accumulation_data.bin", "dhdt_data.bin",
            "dhdt_acc_mask.bin", "udata.bin", "udata_mask.bin",
            "u_iszero.bin", "vdata.bin", "vdata_mask.bin", "v_iszero.bin",
            "temps.bin", "sigma_grid.bin",
        ]
        for f in expected_files
            @test isfile(joinpath(dir, f))
        end

        # Verify round-trip for a 2-D field
        data_read = Vector{Float64}(undef, nx * ny)
        read!(joinpath(dir, "thickness.bin"), data_read)
        @test reshape(data_read, nx, ny) ≈ Gh.h_clip
    end

    # ── Test NetCDF output ───────────────────────────────────────────
    @testset "write_netcdf_file" begin
        dir = mktempdir()
        ncpath = WAVIConstructor.OutputWriting.write_netcdf_file(Gh, Gu, Gv, dir;
                                                                  overwrite = true)
        @test isfile(ncpath)

        NCDataset(ncpath) do ds
            # Global attributes
            @test ds.attrib["Conventions"] == "CF-1.8"
            @test ds.attrib["source"]      == "WAVIConstructor.jl"
            @test ds.attrib["grid_spacing_m"] == 1000.0

            # Coordinate variables
            @test haskey(ds, "x")
            @test haskey(ds, "y")
            @test haskey(ds, "sigma")
            @test length(ds["x"]) == nx
            @test length(ds["y"]) == ny
            @test length(ds["sigma"]) == 5
            @test ds["sigma"][:] ≈ [0.0, 0.25, 0.5, 0.75, 1.0]

            # Check a few data variables and their attributes
            @test haskey(ds, "thickness")
            @test ds["thickness"].attrib["units"] == "m"
            @test ds["thickness"].attrib["long_name"] == "Ice thickness"
            @test Array(ds["thickness"]) ≈ Gh.h_clip

            @test haskey(ds, "temperature")
            @test ds["temperature"].attrib["units"] == "K"
            @test Array(ds["temperature"]) ≈ Gh.levels.temperature_clip

            @test haskey(ds, "udata")
            @test haskey(ds, "vdata")
        end
    end

    # ── Test write_output dispatcher ─────────────────────────────────
    @testset "write_output dispatcher" begin
        # :bin only
        dir_bin = mktempdir()
        WAVIConstructor.OutputWriting.write_output(Gh, Gu, Gv, dir_bin; format=:bin)
        @test isfile(joinpath(dir_bin, "thickness.bin"))
        @test !isfile(joinpath(dir_bin, "wavi_input.nc"))

        # :netcdf only
        dir_nc = mktempdir()
        WAVIConstructor.OutputWriting.write_output(Gh, Gu, Gv, dir_nc; format=:netcdf)
        @test !isfile(joinpath(dir_nc, "thickness.bin"))
        @test isfile(joinpath(dir_nc, "wavi_input.nc"))

        # :both
        dir_both = mktempdir()
        WAVIConstructor.OutputWriting.write_output(Gh, Gu, Gv, dir_both; format=:both)
        @test isfile(joinpath(dir_both, "thickness.bin"))
        @test isfile(joinpath(dir_both, "wavi_input.nc"))

        # invalid format
        @test_throws ErrorException WAVIConstructor.OutputWriting.write_output(
            Gh, Gu, Gv, mktempdir(); format=:csv)
    end

    # ── Test ConstructorParams output_format field ───────────────────
    @testset "ConstructorParams output_format" begin
        p_default = WAVIConstructor.default_constructor_params()
        @test p_default.output_format == :netcdf

        p_nc = WAVIConstructor.default_constructor_params(output_format = :netcdf)
        @test p_nc.output_format == :netcdf

        d = WAVIConstructor.to_dict(p_nc)
        @test d[:output_format] == :netcdf

        p_both = WAVIConstructor.default_constructor_params(output_format = :both)
        @test p_both.output_format == :both
    end
end
