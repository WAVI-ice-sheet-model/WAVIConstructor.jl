using WAVIConstructor.DataLoading
using NCDatasets
using Test

# Create a mock ALBMAP NetCDF file for testing
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
		# Write data as (nx, ny) to match the variable dimensions
		v[:, :] .= ones(Float32, nx, ny)
	end
	close(ds)

	println("Testing get_albmap on: ", mockfile)
	Gh = get_albmap(mockfile)
	@test size(Gh[:h]) == (1, 12)
	@info "Updated size of Gh[:h]: $(size(Gh[:h]))"
	@test all(Gh[:h] .== -1.0f0)  # since usrf-lsrf-firn = 1-1-1 = -1
	println("h size: ", size(Gh[:h]))
	println("b size: ", size(Gh[:b]))
	println("Tma min/max: ", minimum(Gh[:Tma]), ", ", maximum(Gh[:Tma]))
end
