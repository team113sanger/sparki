# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## ~~[0.0.3] - 2024-11-04~~

### Warning
- This tag has been removed because there was a small issue with the CLI.

### Changed 
- Re-organised code in `R/`
- Improved documentation
- Updated `renv` to contain all packages required by SPARKI.

### Added 
- Users can now specify their domains of interest (`--domain` via the CLI or by using the function `plotMinimisers_dotplot` with the argument `domain`)

### Fixed

## [0.0.2] - 2024-08-27

### Changed 
- Optimised and simplified the code to allow SPARKI to run faster
- Several improvements made to the CLI

### Added 
- New function `mergeReports` allows standard and MPA-style reports to be merged into a single dataframe, which now enables users to handle a single dataframe throughout the SPARKI workflow rather than two separate ones
- New plotting functions (`plotClassificationProportion`, `plotDistribution_histogram`, `plotDistribution_violin`, and `plotSignificanceSummary`)
- Users can now specify the organism being analysed (`--organism` via the CLI or by using the function `subsetReports` with the argument `species`)
- Users can now specify samples to be excluded from the SPARKI analysis (`--samples-to-remove` via the CLI or by using the function `mergeReports` with the argument `samples_to_remove`)

## [0.0.1] - 2024-07-19

### Added 
- Initial tag