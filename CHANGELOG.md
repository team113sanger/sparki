# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.2] - 2025-03-04

### Fixed
- Fixed the logging of the function `check_species_in_report()` as it was bailing out in some cases, specifically when the logging threshold was set to `debug` or `trace`.

### Changed
- `exceptions.R` is now `checks.R`
- The function `check_species_in_report()` was moved from `reports-subset.R` to `checks.R` and was refactored; a new function called `look_for_alternative_name()` was created in the refactoring process to make `check_species_in_report()` simpler.

### Added
- Implemented more tests around functions in `reports-subset.R` and added some end-to-end tests for the CLI (`cli.R`).

## [0.1.1] - 2025-02-25

### Warning
- Users are advised not to use `debug` and `trace` as the verbosity level in this version as they may come across a bug that was introduced in the logging of the function `check_species_in_report()`!

### Fixed
- Fixed part of the logging-related code in `cli()` to actually print the argument name when the CLI is launched with a missing required argument.
- `subsetReports()` is now an exported function.

### Changed
- Made `subsetReports()` more permissive so that it can handle both white space-delimited and underscore-delimited species names.

### Added
- Implemented tests around functions in `reports-subset.R`.

## [0.1.0] - 2024-12-10

### Changed
- Improved logging system; as part of this, the flag "--verbose" has been replaced with the option "--verbosity" through which the user can specify the verbosity level.
- Improved the documentation.

### Added
- Implemented tests around functions in `reports-add-info.R`, `exceptions.R`, `load-data.R`, and `main.R`.


## [0.0.4] - 2024-11-06

### Fixed
- Fixed paths in the CLI.

## ~~[0.0.3] - 2024-11-04~~

### Warning
- This tag has been removed because there was a small issue with the CLI.

### Changed
- Re-organised code in `R/`.
- Improved documentation.
- Updated `renv` to contain all packages required by SPARKI.

### Added
- Users can now specify their domains of interest (`--domain` via the CLI or by using the function `plotMinimisers_dotplot()` with the argument `domain`)

### Fixed

## [0.0.2] - 2024-08-27

### Changed
- Optimised and simplified the code to allow SPARKI to run faster.
- Several improvements made to the CLI.

### Added
- New function `mergeReports()` allows standard and MPA-style reports to be merged into a single dataframe, which now enables users to handle a single dataframe throughout the SPARKI workflow rather than two separate ones.
- New plotting functions (`plotClassificationProportion()`, `plotDistribution_histogram()`, `plotDistribution_violin()`, and `plotSignificanceSummary()`).
- Users can now specify the organism being analysed (`--organism` via the CLI or by using the function `subsetReports()` with the argument `species`).
- Users can now specify samples to be excluded from the SPARKI analysis (`--samples-to-remove` via the CLI or by using the function `mergeReports()` with the argument `samples_to_remove`).

## [0.0.1] - 2024-07-19

### Added
- Initial tag.
