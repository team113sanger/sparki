# ✨ SPARKI: a tool for the interpretation of pathogen identification results ✨

This repository contains the code related to SPARKI (**S**tatistical **P**rocess **A**imed at **R**obust **K**raken2 **I**nterpretation), a framework developed in R to help collate, refine, visualise & interpret [Kraken2](https://github.com/DerrickWood/kraken2) outputs.

|                         Main                         |                         Develop                          |
| :----------------------------------------------------: | :------------------------------------------------------: |
| [![pipeline status][main-pipe-badge]][main-branch] | [![pipeline status][develop-pipe-badge]][develop-branch] |

[main-pipe-badge]: https://gitlab.internal.sanger.ac.uk/team113_projects/jb62-projects/sparki/badges/main/pipeline.svg
[main-branch]: https://gitlab.internal.sanger.ac.uk/team113_projects/jb62-projects/sparki/-/commits/main
[develop-pipe-badge]: https://gitlab.internal.sanger.ac.uk/team113_projects/jb62-projects/sparki/badges/develop/pipeline.svg
[develop-branch]: https://gitlab.internal.sanger.ac.uk/team113_projects/jb62-projects/sparki/-/commits/develop

## Table of contents
- [Before you get started](#before-you-get-started)
  - [An overview of Kraken2](#an-overview-of-kraken2)
  - [How SPARKI works](#how-sparki-works)
- [Quick start](#quick-start)
  - [Basic usage example](#basic-usage-example)
- [Installation](#installation)
  - [Install `SPARKI` into your R packages](#install-sparki-into-your-r-packages)
    - [Installing with `remotes`](#installing-with-remotes)
    - [Installing with `renv`](#installing-with-renv)
    - [Installing a specific tag, branch or commit of `SPARKI`](#installing-a-specific-tag-branch-or-commit-of-sparki)
  - [Using Docker](#using-docker)
    - [Getting a Docker image](#getting-a-docker-image)
      - [Building an image from the Dockerfile](#building-an-image-from-the-dockerfile)
      - [Pulling an existing image](#pulling-an-existing-image)
- [Usage](#usage)
  - [Using `SPARKI` as a command line tool (recommended!)](#using-sparki-as-a-command-line-tool-recommended)
  - [Running `SPARKI` with a Docker image](#using-sparki-with-a-docker-image)
  - [Using `SPARKI` as an R package](#using-sparki-as-an-r-package)
- [Additional information](#additional-information)

## Before you get started

SPARKI is a tool to help collate and interpret the outputs produced by Kraken2. In this context, we have developed an end-to-end pathogen identification pipeline, `sparki-nf`, which integrates [Kraken2](https://github.com/DerrickWood/kraken2), [KrakenTools](https://github.com/jenniferlu717/KrakenTools), and SPARKI; additionally, we also provide an optional pipeline, `map-to-genome`, which users can leverage to further validate SPARKI hits.

| Tool | Repository | Goal |
| --- | --- | --- |
| SPARKI | [SPARKI repo](https://github.com/team113sanger/sparki) | Framework to help interpret Kraken2 outputs |
| sparki-nf | [sparki-nf repo](https://github.com/team113sanger/sparki-nf) | Pipeline integrating Kraken2, KrakenTools, and SPARKI |
| map-to-genome | [map-to-genome repo](https://github.com/team113sanger/map-to-genome) | Pipeline for validation of SPARKI hits |

Before using `SPARKI`, please ensure you have run [Kraken2](https://github.com/DerrickWood/kraken2) and [KrakenTools](https://github.com/jenniferlu717/KrakenTools), as their outputs are needed by `SPARKI`. Alternatively, you can use `sparki-nf` for a seamless end-to-end pathogen identification analysis.

### An overview of Kraken2 :bug:

Kraken2 [(Wood *et al*., 2019)](https://github.com/DerrickWood/kraken2) splits the sequencing data from a FASTQ file into *k*-mers, from which minimisers are obtained. By calculating a compact hash code for each minimiser, [Kraken2](https://github.com/DerrickWood/kraken2) is then able to efficiently access a database (DB) of taxa and assign each k-mer the appropriate lower common ancestor taxon. When run with the `--report` mode, Kraken2 generates a sample report (herein referred to as 'standard' report) containing all taxa, at different taxonomic ranks, which were identified in the sample; furthermore, if run with the flag `--report-minimizer-data`, the tool also outputs the number of unique minimisers associated with each taxon that were found in the sample. Alternatively, Kraken2 can be run with the `--report` mode and the flag `--use-mpa-style` to generate MetaPhlAn2 (MPA)-style reports. MPA-style reports can also be generated from 'standard' reports with KrakenTools [(Lu *et al*., 2022)](https://github.com/jenniferlu717/KrakenTools).

For more details on [Kraken2](https://github.com/DerrickWood/kraken2) and [KrakenTools](https://github.com/jenniferlu717/KrakenTools), please refer to their respective repositories and publications.

Pre-compiled [Kraken2](https://github.com/DerrickWood/kraken2) reference DBs are provided by the authors of the software and can be downloaded [here](https://benlangmead.github.io/aws-indexes/k2). Please note that DBs published from 2025 are not currently compatible with `SPARKI`, so we kindly ask that you download and use DBs that were published up to December 2024.

### How `SPARKI` works :sparkles:

`SPARKI` takes as input standard and MPA-style reports from one of more samples and collates all results into a single dataframe. After that, it will calculate, for each taxon, (i) the proportion of minimisers found in a sample (the higher the proportion of minimisers, the more of the organism’s genome is likely present in the sample) and (ii) p-values and adjusted p-values estimated based on a normal distribution and representing the statistical significance of the result.

To ensure reliability, we have implemented tests around the `SPARKI` codebase - please check the section [Additional information](#additional-information) for more details.

## Quick start

`SPARKI` has a command line interface (CLI) that can be interacted with using `Rscript`:

```bash
# To confirm the installation & version
Rscript -e "SPARKI::cli()" --version

# To see how SPARKI can be run
Rscript -e "SPARKI::cli()" --help
```

### Basic usage example

```bash
Rscript -e "SPARKI::cli()" \
  --std-reports reports/ \
  --mpa-reports mpa/ \
  --organism "Homo sapiens" \
  --reference inspect.txt \
  --domain "Viruses" \
  --outdir .
```

## Installation

### Install `SPARKI` into your R packages

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

#### Getting a Docker image

##### Building an image from the Dockerfile

```bash
docker build -t sparki:local -f Dockerfile .
```

##### Pulling an existing image

```bash
docker pull quay.io/team113sanger/sparki:latest
```

You can find more details on how to interact with `SPARKI`'s CLI using Docker in the section [Running SPARKI with a Docker image](#running-sparki-with-a-docker-image).

## Usage

### Using `SPARKI` as a command line tool (recommended!)

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
PROJECTDIR="path/to/project_dir"

Rscript -e "SPARKI::cli()" \
    --std-reports ${PROJECTDIR}/reports \
    --mpa-reports ${PROJECTDIR}/mpa \
    --organism "Homo sapiens" \
    --reference ${PROJECTDIR}/inspect.txt \
    --metadata ${PROJECTDIR}/metadata.csv \
    --sample-col sample_id \
    --columns column1,column2 \
    --prefix project \
    --outdir ${PROJECTDIR}/outputs/ \
    --verbosity info \
    --domain Viruses,Bacteria \
    --samples-to-remove ${PROJECTDIR}/samples_remove.txt
```

### Running `SPARKI` with a Docker image

After following the steps outlined in the section [Getting a Docker image](#getting-a-docker-image), you can interact with `SPARKI`'s CLI as follows:

```bash
docker run --rm sparki:local -v path/to/dir:/opt/data/ Rscript -e "SPARKI::cli()" \
  --std-reports /opt/data/reports \
  --mpa-reports /opt/data/mpa \
  --organism "Homo sapiens" \
  --reference /opt/data/inspect.txt \
  --domain "Viruses" \
  --outdir /opt/data/
```

Note that you will need to use a bind mount (as in `-v path/to/dir:/opt/data/`) to mount a directory from your machine into the container, so that your Kraken2 outputs and related files are visible to the container.

### Using `SPARKI` as an R package

Please check out the tutorial in this repository, which can be found at `vignettes/SPARKI_basic_usage.pdf`.

## Additional information

<details>

<summary>For Sanger users</summary>

As a convenience, SPARKI is installed as a module on the Sanger farm. To load the latest SPARKI module:

```bash
module load /software/team113/modules/modulefiles/sparki/default
```

</details>

<details>

<summary>For developers</summary>

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

### Testing

Unit and integration tests for `SPARKI` can be found in this repository at `tests/`. To run them, you can enter the [development-focused container created above](#development-focused-environment), start an R shell, and then run the following command:

```R
devtools::test()
```

</details>
