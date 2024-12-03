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
- [Installation - quick start](#installation---quick-start)
    - [Installing with `remotes`](#installing-with-remotes)
    - [Installing with `renv`](#installing-with-renv)
    - [Installing a specific tag, branch or commit of `SPARKI`](#installing-a-specific-tag-branch-or-commit-of-sparki)

## Installation - quick start

This package requires a **Personal Access Token** (PAT) to be able to install it from the **Sanger GitLab**. You can generate a PAT by following the instructions [here](https://gitlab.internal.sanger.ac.uk/-/user_settings/personal_access_tokens).

You can install SPARKI using the `remotes` or `renv` package.

If you are installing SPARKI on the Sanger farm you can use R version 4.2 or
newer. For convenience, you can load the R 4.4 module with the following command:

```bash
module load rocker/rver/4.4.0
```

### Installing with `remotes`

Start an R shell, install the `remotes` package if you haven't already, and then
install SPARKI:

```R
PAT <- "glpat-...<your-pat-goes-here>..."
install.packages("remotes")
remotes::install_gitlab(
  repo = "team113_projects/jb62-projects/sparki",
  host = "gitlab.internal.sanger.ac.uk",
  auth_token = PAT
)
```

### Installing with `renv`

If you are using `renv` to manage your R environment, you can add SPARKI to your
project by running the following commands in an R shell:

```R
# It is assumed that you have already created an `renv` environment
options(renv.config.gitlab.host = "gitlab.internal.sanger.ac.uk")
PAT <- "glpat-...<your-pat-goes-here>..."
Sys.setenv(GITLAB_PAT = PAT)
renv::install("gitlab::team113_projects/jb62-projects/sparki")
```

### Installing a specific tag, branch or commit of `SPARKI`

You can install a specific branch, tag or commit by adding a `@ref` to the end of the repo URL. For example, to install the `develop` branch:

With `remotes`:
```R
remotes::install_gitlab(
  repo = "team113_projects/jb62-projects/sparki@develop",
  host = "gitlab.internal.sanger.ac.uk",
  auth_token = "..."
)
```

With `renv`:
```R
options(renv.config.gitlab.host = "gitlab.internal.sanger.ac.uk")
renv::install("gitlab::team113_projects/jb62-projects/sparki@develop")
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
- `--sample-col`: if `--metadata` is provided, users should also specify the name of the column that contains sample IDs.
- `--columns`: if `--metadata` is provided, users should also specify the names of the columns from the metadata table that should be used; the column names must be comma-separated.
- `--prefix`: prefix to be added to SPARKI's output files (note: you do not need to include an underscore at the end of your prefix!).
- `--samples-to-remove`: if the directories provided to `--std-reports` and `--mpa-reports` contain samples that should not be included in the final outputs, a list with those samples can be provided in a text file.
- `--verbose`: flag indicating that the user would like log messages to be printed out.
- `--include-eukaryotes`: flag indicating that the user would like eukaryotes to be included in all plots (note that using this flag may cause some plots to be too full; also note that some plots will have eukaryotes included regardless of this flag).
- `--include-sample-names`: flag indicating that the user would like sample names to be included in all plots (note that using this flag may cause some plots to be too full; also note that some plots will have sample names included regardless of this flag).

### CLI usage example

```
PROJECTDIR="/lustre/scratch126/casm/team113da/users/jb62/projects/sparki"

/software/team113/dermatlas/R/R-4.2.2/bin/Rscript ${PROJECTDIR}/R/cli.R \
    --std-reports ${PROJECTDIR}/test/reports \
    --mpa-reports ${PROJECTDIR}/test/mpa \
    --organism "Homo sapiens" \
    --reference ${PROJECTDIR}/test/inspect.txt \
    --metadata ${PROJECTDIR}/test/metadata.csv \
    --sample-col Tumour_RNA \
    --columns Diagnosis_short,Site_group \
    --prefix SebT \
    --outdir ${PROJECTDIR}/test/outputs/ \
    --verbose \
    --domain Viruses,Bacteria \
    --samples-to-remove ${PROJECTDIR}/test/samples_remove.txt
```


## Using SPARKI inside R

Please check out [this tutorial](https://gitlab.internal.sanger.ac.uk/team113_projects/jb62-projects/sparki/-/blob/develop/tutorials/SPARKI_basic_usage.pdf?ref_type=heads).
