# Changelog

All notable changes to this project are documented here.
The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project follows [SemVer](https://semver.org/).

## [Unreleased]

## [0.1.0] - YYYY-MM-DD

First release.

### Added
- Singleton-type dispatch for data sources (`BedMachineV3`, `ALBMAPv1`, `FrankTemps`, `BISICLESTemps`, `MEaSUREs`, `ArthernAccumulation`, `SmithDhdt`, `ZwallyBasins`) with `load_data`, `default_path` and `interpolate_temperature`.
- `ConstructorParams`, `default_constructor_params` and `minimal_constructor_params` for configuring a run.
- `init_bedmachine` (grid initialisation), `select_domain_wavi` (basin-based domain selection) and `setup_wavi_data` (complete workflow).
- Output writing as NetCDF (default) or binary files via `write_output`, `write_netcdf_file` and `write_binary_files`.
- Documentation, examples and CI (tests on Julia 1.11 and 1.12, Runic formatting, docs build and deployment).
