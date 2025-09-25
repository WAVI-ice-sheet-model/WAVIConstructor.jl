using WAVIConstructor.DataLoading
using ArchGDAL
using Test

# Create a mock GeoTIFF file for testing
mktempdir() do tmpdir
    mockfile = joinpath(tmpdir, "mock_geotiff.tif")
    width, height = 4, 3
    dx, dy = 10.0, -10.0
    mapx, mapy = 100.0, 200.0

    ArchGDAL.create(mockfile, "GTiff", width, height, 1, ArchGDAL.GDT_Float32) do dataset
        ArchGDAL.setgeotransform!(dataset, [mapx, dx, 0.0, mapy, 0.0, dy])
        ArchGDAL.write!(ArchGDAL.getband(dataset, 1), ones(Float32, width, height))
    end

    println("Testing geotiff_read_axis_only on: ", mockfile)
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
