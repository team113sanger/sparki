# ✨ SPARKI: a tool for the interpretation of pathogen identification results ✨

This repository contains the code related to SPARKI (**S**tatistical **P**rocess **A**imed at **R**obust **K**raken2 **I**nterpretation), a framework developed in R to help collate, refine, visualise & interpret Kraken2 outputs.

|                         Main                         |                         Develop                          |
| :----------------------------------------------------: | :------------------------------------------------------: |
| [![pipeline status][main-pipe-badge]][main-branch] | [![pipeline status][develop-pipe-badge]][develop-branch] |

[main-pipe-badge]: https://gitlab.internal.sanger.ac.uk/team113_projects/jb62-projects/sparki/badges/main/pipeline.svg
[main-branch]: https://gitlab.internal.sanger.ac.uk/team113_projects/jb62-projects/sparki/-/commits/main
[develop-pipe-badge]: https://gitlab.internal.sanger.ac.uk/team113_projects/jb62-projects/sparki/badges/develop/pipeline.svg
[develop-branch]: https://gitlab.internal.sanger.ac.uk/team113_projects/jb62-projects/sparki/-/commits/develop

## Table of contents
- [Quick start](#quick-start)
- [Installation](#installation)
  - [Install it into your R packages](#install-it-into-your-r-packages)
    - [Installing with `remotes`](#installing-with-remotes)
    - [Installing with `renv`](#installing-with-renv)
    - [Installing a specific tag, branch or commit of `SPARKI`](#installing-a-specific-tag-branch-or-commit-of-sparki)
  - [Using Docker](#using-docker)
- [Usage](#usage)
  - [Using `SPARKI` as a command line tool](#using-sparki-as-a-command-line-tool)
  - [Using `SPARKI` inside R](#using-sparki-inside-r)
- [For developers](#for-developers)
  - [Making a release](#making-a-release)
  - [Development-focused environment](#development-focused-environment)

## Quick start

SPARKI is an R package with a command line interface (CLI) that can be interacted with using `Rscript`:

```bash
# To confirm the installation & version!
Rscript -e "SPARKI::cli()" --version

# To see how SPARKI can be run!
Rscript -e "SPARKI::cli()" --help
```

## Installation

### Install it into your R packages

You can install SPARKI using the `remotes` or `renv` packages, as described below.

#### Installing with `remotes`

Start an R shell, install the `remotes` package if you haven't already, and then install `SPARKI`:

```R
remotes::install_github("team113sanger/sparki")
```

#### Installing with `renv`

If you are using `renv` to manage your R environment, you can add `SPARKI` to your project by running the following commands in an R shell:

```R
# It is assumed that you have already created an `renv` environment
renv::install("team113sanger/sparki")
```

#### Installing a specific tag, branch or commit of `SPARKI`

You can install a specific branch, tag or commit by adding a `@ref` to the end of the repo URL. For example, to install the `develop` branch:

With `remotes`:
```R
remotes::install_github("team113sanger/sparki@develop")
```

With `renv`:
```R
renv::install("team113sanger/sparki@develop")
```

### Using Docker

Alternatively, you can use a Docker container to interact with `SPARKI`'s CLI. In this repository we provide a `Dockerfile` that you can use to start a container and run a `SPARKI` analysis following the instructions below:

```bash
# 1 - Create an image from the Dockerfile.
docker build -t sparki:local -f Dockerfile .

# 2.1 - Check the SPARKI help.
docker run --rm sparki:local Rscript -e "SPARKI::cli()" --help

# 2.2 - Run SPARKI.
docker run --rm sparki:local -v path/to/dir:/opt/data/ Rscript -e "SPARKI::cli()" \
  --std-reports /opt/data/reports \
  --mpa-reports /opt/data/mpa \
  --organism "Homo sapiens" \
  --reference /opt/data/inspect.txt \
  --domain "Viruses" \
  --outdir /opt/data/
```

You can find more details on how to use `SPARKI`'s CLI below.

## Usage

### Using SPARKI as a command line tool (recommended!)

#### CLI arguments/options

##### *Required*
- `--std-reports`: path to a directory containing 'standard' reports.
- `--mpa-reports`: path to a directory containing MPA-style reports.
- `--organism`: species you are analysing (e.g. if working with human samples, this will be `Homo sapiens`).
- `--reference`: path to the `inspect.txt` file that is present in the reference database used to run Kraken2.
- `--domain`: domain(s) that the user is interested in (e.g. `Viruses`); note that if multiple domains are to be provided, they must be comma-separated (e.g. `Viruses,Bacteria`).
- `--outdir`: path to an output directory (note that it must be an empty directory!).

##### *Optional*
- `--metadata`: path to a metadata file containing sample-level information.
- `--sample-col`: if `--metadata` is provided, users should also specify the name of the column that contains sample IDs.
- `--columns`: if `--metadata` is provided, users should also specify the names of the columns from the metadata table that should be used; the column names must be comma-separated.
- `--prefix`: prefix to be added to SPARKI's output files (note: you do not need to include an underscore at the end of your prefix!).
- `--samples-to-remove`: if the directories provided to `--std-reports` and `--mpa-reports` contain samples that should not be included in the final outputs, a list with those samples can be provided in a text file; each line in the file should be a sample ID.
- `--verbosity`: verbosity level (one of `trace`/`t`, `debug`/`d`, `info`/`i`, `success`/`s`, `warn`/`w`, `error`/`e`, `fatal`/`f`, or `off`/`o`).
- `--include-eukaryotes`: flag indicating that the user would like eukaryotes to be included in all plots (note that using this flag may cause some plots to be too full; also note that some plots will have eukaryotes included regardless of this flag).
- `--include-sample-names`: flag indicating that the user would like sample names to be included in all plots (note that using this flag may cause some plots to be too full; also note that some plots will have sample names included regardless of this flag).

#### CLI usage example

```bash
PROJECTDIR="/lustre/scratch126/casm/team113da/users/jb62/projects/sparki"
module load sparki # If working on the farm!

Rscript -e "SPARKI::cli()" \
    --std-reports ${PROJECTDIR}/test/reports \
    --mpa-reports ${PROJECTDIR}/test/mpa \
    --organism "Homo sapiens" \
    --reference ${PROJECTDIR}/test/inspect.txt \
    --metadata ${PROJECTDIR}/test/metadata.csv \
    --sample-col Tumour_RNA \
    --columns Diagnosis_short,Site_group \
    --prefix SebT \
    --outdir ${PROJECTDIR}/test/outputs/ \
    --verbosity info \
    --domain Viruses,Bacteria \
    --samples-to-remove ${PROJECTDIR}/test/samples_remove.txt
```


### Using SPARKI as an R package

Please check out [this tutorial](https://gitlab.internal.sanger.ac.uk/team113_projects/jb62-projects/sparki/-/blob/develop/tutorials/SPARKI_basic_usage.pdf?ref_type=heads).

<details>

<summary> For Sanger users

As a convenience, SPARKI is installed as a module on the Sanger farm. To load the latest SPARKI module:

```bash
module load /software/team113/modules/modulefiles/sparki/default
```

</details>

<details>

<summary> For developers

### Making a release

To make a new release, please follow the steps below:

1. Start a release with HubFlow.
```bash
module load git
git hf release start 0.1.2 # Or whatever tag.
```

2. Make documentation commits.
```bash
git add DESCRIPTION
git commit -m "Bump package version to 0.1.2"
```
```bash
git add CHANGELOG.md
git commit -m "Update changelog"
```
```bash
git hf push
```

3. Finish the release.
```bash
git hf release finish 0.1.2 # Or whatever the name of the release.
```

4. Just follow the next steps that will be displayed in your terminal after running `git hf release finish` and that's it!

### Development-focused environment

If you are working on `SPARKI`, please do so using a development-focused Docker container that you can start by following the steps below:

```bash
# 1 - Create an image and start a container.
docker compose -f docker-compose.yml up -d --build

# 2 - Enter the container.
docker exec -it sparki bash
```

</details>
