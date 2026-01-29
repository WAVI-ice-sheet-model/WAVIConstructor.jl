# WAVIConstructor.jl Examples

This directory contains example scripts demonstrating how to use WAVIConstructor.jl.

## Quick Start

**Start here:** `minimal_test.jl`
- ✅ Fast (30-60 seconds)
- ✅ Low memory usage
- ✅ Safe to run interactively
- Uses coarse grid (32 km) and minimal datasets

```julia
julia --project=. examples/minimal_test.jl
```

## Full Examples

**After testing:** `basic_usage.jl`
- ⚠️ Computationally expensive (5-30 minutes)
- ⚠️ High memory usage (several GB)
- ⚠️ Run on compute node, not interactively
- Shows complete workflows with all datasets

The examples in `basic_usage.jl` are **commented out by default**. Uncomment sections to run them.

## Reference

`quick_reference.jl` - Parameter reference and function documentation

## Computational Cost

| Grid Spacing | Grid Size | RAM Usage | Runtime | Recommendation |
|--------------|-----------|-----------|---------|----------------|
| 64 km        | ~200×200  | ~500 MB   | ~15 sec | Quick tests |
| 32 km        | ~400×400  | ~1-2 GB   | ~30 sec | Testing |
| 16 km        | ~800×800  | ~3-5 GB   | ~2 min  | Development |
| 10 km        | ~1300×1300| ~8-15 GB  | ~10 min | Production |
| 5 km         | ~2600×2600| ~30+ GB   | ~1 hour | High-res (HPC) |

**Tip:** Start with coarse grids (32-64 km) to test your workflow, then increase resolution.

## Dataset Requirements

The package can run with just BedMachine:
- **Minimal:** Only `BedMachineAntarctica-v3.nc` (required)
- **Full:** BedMachine + velocity + temperature + dhdt + accumulation + basins

Missing datasets will be replaced with zeros/defaults with a warning.

## Running on HPC

For production runs, submit as a batch job:

```bash
#!/bin/bash
#SBATCH --job-name=wavi_setup
#SBATCH --time=01:00:00
#SBATCH --mem=16GB
#SBATCH --cpus-per-task=1

module load julia
cd /path/to/WAVIConstructor.jl
julia --project=. examples/basic_usage.jl  # After uncommenting the examples
```
