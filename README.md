# :sparkles: SPARKI: a tool for the interpretation of pathogen identification results :sparkles:

This repository contains the code related to SPARKI (**S**tatistical **P**rocess **A**imed at **R**obust **K**raken2 **I**nterpretation), a framework developed in R to help collate, refine, visualise & interpret Kraken2 outputs.

## Installation
```
git clone git@gitlab.internal.sanger.ac.uk:team113_projects/jb62-projects/sparki.git
cd sparki
```

## Using SPARKI's command line interface (CLI)

### CLI arguments/options

#### *Required*
- `--std-reports`: path to a directory containing 'standard' reports.
- `--mpa-reports`: path to a directory containing MPA-style reports.
- `--organism`: species you are analysing (e.g. if working with human samples, this will be `Homo sapiens`).
- `--reference`: path to the `inspect.txt` file that is present in the reference database used to run Kraken2.
- `--domain`: domain(s) that the user is interested in (e.g. `Viruses`); note that if multiple domains are to be provided, they must be comma-separated (e.g. `Viruses,Bacteria`).
- `--outdir`: path to an output directory.

#### *Optional*
- `--metadata`: path to a metadata file containing sample-level information.
- `--sample-col`: if `--metadata` is provided, users should also specify the name the column that contains sample IDs.
- `--columns`: if `--metadata` is provided, users should also specify the names of the columns from the metadata table that should be used; the column names must be comma-separated.
- `--prefix`: prefix to be added to SPARKI's output files.
- `--samples-to-remove`: if the directories provided to `--std-reports` and `--mpa-reports` contain samples that should not be included in the final outputs, a list with those samples can be provided in a text file.
- `--verbose`: flag indicating that the user would like log messages to be printed out.
- `--include-eukaryotes`: flag indicating that the user would like eukaryotes to be included in all plots (note that using this flag may cause some plots to be too full; also note that some plots will have eukaryotes included regardless of this flag).
- `--include-sample-names`: flag indicating that the user would like sample names to be included in all plots (note that using this flag may cause some plots to be too full; also note that some plots will have sample names included regardless of this flag).
