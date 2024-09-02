# SPARKI :sparkles:

This repository contains the code related to SPARKI (**S**tatistical **P**rocess **A**imed at **R**obust **K**raken2 **I**nterpretation), a framework developed in R to help refine, visualise & interpret Kraken2 outputs.

## Installation
```
git clone git@gitlab.internal.sanger.ac.uk:team113_projects/jb62-projects/sparki.git
cd sparki
```

## Tutorial

```
rmarkdown::render(
    "/lustre/scratch126/casm/team113da/users/jb62/projects/sparki/tutorials/SPARKI_basic_usage.Rmd",
    output_dir = "/lustre/scratch126/casm/team113da/users/jb62/projects/sparki/tutorials",
    output_file = "SPARKI_basic_usage.pdf"
)
```

# CLI usage
```
PROJECTDIR="/lustre/scratch126/casm/team113da/users/jb62/projects/sparki"
/software/team113/dermatlas/R/R-4.2.2/bin/Rscript ${PROJECTDIR}/src/cli.R \
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
    --domain Viruses,Bacteria
```
