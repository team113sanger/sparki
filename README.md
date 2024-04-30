# SPARKI (Statistical Process Aimed at Robust Kraken2 Interpretation) :sparkles:

This repository contains the code related to SPARKI (**S**tatistical **P**rocess **A**imed at **R**obust **K**raken2 **I**nterpreation), a framework developed in R to help refine, visualise & interpret Kraken2 outputs. This package is tailored for the identification of viruses

## Tutorial

### Import functionalities
```
source("library/utilities.R")
source("library/helper.R")
source("library/plotting.R")
```

### Read Kraken2 results
```
merged_mpa <- load_mpa("test/mpa", verbose = FALSE)
merged_reports <- load_reports("test/reports", verbose = FALSE)
```

### Read metadata (optional) and Kraken2's reference database
```
mdata <- load_metadata("test/metadata.csv") # Optional
ref_db <- load_reference("test/inspect.txt")
```

### Add metadata to dataframes containing Kraken2 results *(optional)*
```
mdata_columns <- c("Diagnosis_short", "Site_group")

merged_mpa <- add_metadata(merged_mpa, mdata, mdata_columns)
merged_reports <- add_metadata(merged_reports, mdata, mdata_columns)
```

### Prepare standard reports and MPA-style reports for downstream analysis
```
merged_mpa <- mpa_determine_last_rank(merged_mpa)
merged_mpa <- mpa_add_rank_columns(merged_mpa)

merged_reports <- subset_reports(merged_reports)
merged_reports <- report_add_domains(merged_reports, merged_mpa)
```

