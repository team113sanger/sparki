# SPARKI (Statistical Process Aimed at Robust Kraken2 Interpretation) :sparkles:

This repository contains the code related to SPARKI (**S**tatistical **P**rocess **A**imed at **R**obust **K**raken2 **I**nterpreation), a framework developed in R to help refine, visualise & interpret Kraken2 outputs. This package is **tailored/fine-tuned for the interpretation of results related to viruses**.

## Installation
```
git clone git@gitlab.internal.sanger.ac.uk:team113_projects/jb62-projects/sparki.git
cd sparki
```

## Tutorial

### Import functionalities
```
source("R/utilities.R")
source("R/helper.R")
source("R/plotting.R")
source("R/constants.R")
```

### Read Kraken2 results
```
mpa_reports <- load_mpa_reports("test/mpa", verbose = FALSE)
std_reports <- load_std_reports("test/reports", verbose = FALSE)
```

### Read metadata (optional) and Kraken2's reference database
```
mdata <- load_metadata("test/metadata.csv") # Optional
ref_db <- load_reference("test/inspect.txt")
```

### Add metadata to dataframes containing Kraken2 results *(optional)*
```
mdata_columns <- c("Diagnosis_short", "Site_group")

mpa_reports <- addMetadata(mpa_reports, mdata, mdata_columns)
std_reports <- addMetadata(std_reports, mdata, mdata_columns)
```

### Process MPA-style reports
```
mpa_reports <- addRank(mpa_reports)
mpa_reports <- addConciseTaxon(mpa_reports)
mpa_reports <- transfer_ncbiID(mpa_reports, std_reports)
mpa_reports <- transfer_nReads(mpa_reports, std_reports) # To be run after: std_reports <- add_nReads(std_reports)
mpa_reports <- add_DBinfo(mpa_reports, ref_db)
```

### Process standard reports
```
std_reports <- add_nReads(std_reports)
std_reports <- add_DBinfo(std_reports, ref_db)
std_reports <- transferDomains(std_reports, mpa_reports)
```

### Create auxiliary dataframes
```
nDomainReads <- get_nDomainReads(mpa_reports) # Or nDomainReads <- get_nDomainReads(std_reports)
```

n_ranks_in_samples <- mpa_get_n_ranks_in_samples(merged_mpa)

merged_reports <- report_assess_n_taxa_per_sample(merged_reports) # Must be run before subset_reports()
merged_reports <- report_assess_n_reads_per_sample(merged_reports) # Must be run before subset_reports()
class_unclass_df <- get_class_unclass_numbers(merged_reports) # Must be run before subset_reports()

merged_reports <- subset_reports(merged_reports, include_human = FALSE) # Keep only family-, genus-, and species-level results
merged_reports <- report_add_domains(merged_reports, merged_mpa)
merged_reports <- report_add_db_info(merged_reports, ref_db)
merged_reports <- report_assess_n_taxa_specific_rank_per_sample(merged_reports, rank = c("F", "G", "S"))
merged_reports <- report_assess_statistical_significance(merged_reports, ref_db)
```

### Create summary plots
```
# Number of classified and unclassified reads per sample.
plotClassificationSummary_dotplot()




plot_classified_vs_unclassified_summary(class_unclass_df, "outputs/")
plot_classified_vs_unclassified_per_sample(class_unclass_df)
```

