using WAVIConstructor.DataLoading
using ArchGDAL
using Test

# Create a mock GeoTIFF file for testing
mktempdir() do tmpdir
    mockfile = joinpath(tmpdir, "mock_smith_dhdt.tif")
    width, height = 4, 3
    dx, dy = 10.0, -10.0
    mapx, mapy = 100.0, 200.0

    # Create mock dhdt data
    mock_dhdt = Float32[1.5 2.5 3.5 4.5; -1.0 0.0 1.0 2.0; -2.5 -1.5 -0.5 0.5]

    ArchGDAL.create(mockfile, "GTiff", width, height, 1, ArchGDAL.GDT_Float32) do dataset
        ArchGDAL.setgeotransform!(dataset, [mapx, dx, 0.0, mapy, 0.0, dy])
        ArchGDAL.write!(ArchGDAL.getband(dataset, 1), mock_dhdt)
    end

    println("Testing get_smith_dhdt on: ", mockfile)
    xx, yy, dhdt = get_smith_dhdt(mockfile)

    @test size(xx) == (height, width)
    @test size(yy) == (height, width) 
    @test size(dhdt) == (height, width)

    # Test coordinate calculation (pixel center coordinates)
    expected_x = [mapx + dx * (i - 0.5) for i in 1:width]
    expected_y = [mapy - abs(dy) * (i - 0.5) for i in 1:height]

    @test all(xx[1, :] .≈ expected_x)
    @test all(yy[:, 1] .≈ reverse(expected_y))

    # Test that dhdt data matches the mock data
    @test all(dhdt .≈ mock_dhdt)

    println("dhdt size: ", size(dhdt))
    println("xx range: ", minimum(xx), " to ", maximum(xx))
    println("yy range: ", minimum(yy), " to ", maximum(yy))
    println("dhdt range: ", minimum(dhdt), " to ", maximum(dhdt))
end
