# SPARKI (Statistical Process Aimed at Robust Kraken2 Interpretation) :sparkles:

This repository contains the code related to SPARKI (**S**tatistical **P**rocess **A**imed at **R**obust **K**raken2 **I**nterpreation), a framework developed in R to help refine, visualise & interpret Kraken2 outputs.

## Usage

### Import functionalities
```
source("library/utilities.R")
source("library/helper.R")
source("library/plotting.R")
```

### Read input data
```
merged_mpa <- load_mpa("test/mpa")
merged_reports <- load_mpa("test/reports")
```