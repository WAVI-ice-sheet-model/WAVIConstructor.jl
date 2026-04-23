# Contributing to WAVIConstructor.jl

Thank you for your interest in contributing! This guide covers setting up a development environment, running tests and linters, adding new datasets, and the workflow for branches, pull requests, and releases.

## Development Environment Setup

1. **Clone the repository:**
   ```sh
   git clone https://github.com/WAVI-ice-sheet-model/WAVIConstructor.jl.git
   cd WAVIConstructor.jl
   ```
2. **Activate the project and install dependencies:**
   ```julia
   julia --project
   pkg> instantiate
   ```
   Or, from the Julia REPL:
   ```julia
   using Pkg
   Pkg.activate(".")
   Pkg.instantiate()
   ```

## Running Tests

Run all tests from the project root:
```julia
julia --project
pkg> test
```
Or from the Julia REPL:
```julia
using Pkg
Pkg.test()
```


## Adding a New Dataset


To add a new dataset, you typically need to edit or add code in the following files:

- `src/DataSources.jl` — Define your new source type and its default path.
- `src/DataLoading.jl` — Implement the `load_data` method for your source.
- (If temperature source) `src/SetupData.jl` — Implement `interpolate_temperature` for your source.

**Example: Adding a new temperature dataset**

1. **Define a singleton type** in `src/DataSources.jl`:
  ```julia
  struct MyNewSource <: TemperatureSource end
  ```
2. **Set the default path** in `src/DataSources.jl`:
  ```julia
  WAVIConstructor.default_path(::MyNewSource) = "Data/my_file.nc"
  ```
3. **Implement the loader** in `src/DataLoading.jl`:
  ```julia
  function WAVIConstructor.load_data(::MyNewSource, file)
     # Load your data here (e.g., using NCDatasets)
     return (temps = rand(10,10,10), x = 1:10, y = 1:10, sigma = 1:10)
  end
  ```
4. **(Temperature sources only) Implement interpolation** in `src/SetupData.jl`:

If your temperature data is already on the WAVI grid, you can simply return it:
```julia
function WAVIConstructor.interpolate_temperature(::MyNewSource, data, Gh)
  return data.temps
end
```

If your data is on a different grid, you need to interpolate it to the WAVI grid. Here is a concrete example using the approach from the codebase (see BISICLESTemps and FrankTemps in `src/DataLoading.jl`):

```julia
function WAVIConstructor.interpolate_temperature(::MyNewSource, data, Gh)
  # Assume data.temps is (nz, ny, nx), data.x, data.y, data.sigma are the coordinates
  temperature = zeros(size(data.temps, 1), Gh.nx, Gh.ny)
  for i in 1:size(data.temps, 1)
    # Interpolate each sigma level to the WAVI grid
    temperature[i, :, :] = WAVIConstructor.interpolate_to_grid(
      data.x[:], data.y[:], data.temps[i, :, :][:], Gh.xx, Gh.yy
    )
  end
  return temperature, data.sigma
end
```

Replace `data.x`, `data.y`, `data.sigma`, and `data.temps` with the actual variable names from your loader, and `Gh.x`, `Gh.y`, `Gh.sigma` with the WAVI grid coordinates.

For more advanced examples, see how `interpolate_temperature` is implemented for other sources (BISICLESTemps, FranksTemps) in `src/DataLoading.jl`.

No changes to other files are needed. The dispatch system will pick up your new source automatically.

## Branch, PR Workflow, Versioning, and Releases

- **Branching:**
  - Fork the repo or create a feature branch from `main`.
  - Use descriptive branch names (e.g., `feature/new-dataset`, `fix/typo`).
- **Commits:**
  - Write clear, concise commit messages.
- **Pull Requests:**
  - Open a PR against `main`.
  - Ensure all tests pass and code is formatted.
  - Reference related issues in the PR description.
- **Reviews:**
  - At least one approval is required before merging.
- **Versioning:**
  - Follows [SemVer](https://semver.org/). Bump version in `Project.toml` for breaking changes, features, or fixes.
- **Releases:**
  - Tag a release on GitHub. CI will build docs and run tests.

## Questions?

Open an issue or discussion on GitHub if you need help!
