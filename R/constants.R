######################
##  HELPER CONSTANTS ##
#####################################################################

#### CLI help text ####
CLI_DESCRIPTION <- "SPARKI (Statistical Process Aimed at Robust Kraken2 Interpretation) is a framework to help collate, refine, visualise & interpret Kraken2 outputs."
CLI_PROGRAM_NAME <- "SPARKI"

####  Reference database ####
COLNAME_REF_DB_PCT_FRAG_CLADE <- "pct_fragments_clade"
COLNAME_REF_DB_MINIMISERS_CLADE <- "n_minimisers_clade"
COLNAME_REF_DB_MINIMISERS_TAXON <- "n_minimisers_taxon"
COLNAME_REF_DB_RANK <- "rank"
COLNAME_REF_DB_NCBI_ID <- "ncbi_id"
COLNAME_REF_DB_TAXON <- "taxon"

####  Standard report ####
COLNAME_STD_SAMPLE <- "sample"
COLNAME_STD_PCT_FRAG_CLADE <- "pct_fragments_clade"
COLNAME_STD_N_FRAG_CLADE <- "n_fragments_clade"
COLNAME_STD_LOG_N_FRAG_CLADE <- "log10_n_fragments_clade"
COLNAME_STD_N_FRAG_TAXON <- "n_fragments_taxon"
COLNAME_STD_MINIMISERS <- "n_minimisers"
COLNAME_STD_UNIQ_MINIMISERS <- "n_distinct_minimisers"
COLNAME_STD_RANK <- "rank"
COLNAME_STD_NCBI_ID <- "ncbi_id"
COLNAME_STD_TAXON <- "taxon"
COLNAME_STD_HIERARCHY <- "info_from_mpa_report"
COLNAME_STD_SAMPLE_SIZE <- "sample_size"
COLNAME_STD_DOMAIN <- "domain"
COLNAME_STD_DB_MINIMISERS_CLADE <- "db_n_minimisers_clade"
COLNAME_STD_DB_MINIMISERS_TAXON <- "db_n_minimisers_taxon"
COLNAME_STD_RATIO_CLADE <- "ratio_clade"
COLNAME_STD_RATIO_TAXON <- "ratio_taxon"
COLNAME_STD_P_CLADE_IN_DB <- "p_clade_in_db"
COLNAME_STD_PVALUE <- "pvalue"
COLNAME_STD_PADJ <- "padj"
COLNAME_STD_SIGNIF <- "significance"
COLNAME_STD_N_TAXA_RANK <- "n_taxa_in_rank"

####  MPA-style report ####
COLNAME_MPA_SAMPLE <- "sample"
COLNAME_MPA_N_FRAG_CLADE <- "n_fragments_clade"
COLNAME_MPA_TAXON_HIERARCHY <- "taxon_hierarchy"
COLNAME_MPA_TAXON_LEAF <- "taxon_leaf"
COLNAME_MPA_RANK <- "rank"
COLNAME_MPA_DOMAIN <- "domain"
COLNAME_MPA_KINGDOM <- "kingdom"
COLNAME_MPA_PHYLUM <- "phylum"
COLNAME_MPA_CLASS <- "class"
COLNAME_MPA_ORDER <- "order"
COLNAME_MPA_FAMILY <- "family"
COLNAME_MPA_GENUS <- "genus"
COLNAME_MPA_SPECIES <- "species"


####  Auxiliary dataframes ####

# Number of classified/unclassified reads per sample.
COLNAME_CLASSIF_SUMMARY_SAMPLE <- "sample"
COLNAME_CLASSIF_SUMMARY_READ_TYPE <- "type"
COLNAME_CLASSIF_SUMMARY_N_FRAG <- "n_fragments_clade"
COLNAME_CLASSIF_SUMMARY_LOG_N_FRAG <- "log10_n_fragments_clade"

# Proportion of classified reads per sample.
COLNAME_PROP_SAMPLE <- "sample"
COLNAME_PROP_CLASSIFIED <- "proportion_classified"

# Number of domain reads per sample.
COLNAME_DOMAIN_READS_SAMPLE <- "sample"
COLNAME_DOMAIN_READS_TAXON <- "taxon"
COLNAME_DOMAIN_READS_N_FRAG <- "n_fragments_clade"
COLNAME_DOMAIN_READS_LOG_N_FRAG <- "log10_n_fragments_clade"

#  Significance summary.
COLNAME_SIGNIF_SUMMARY_SAMPLE <- "sample"
COLNAME_SIGNIF_SUMMARY_RANK <- "rank"
COLNAME_SIGNIF_SUMMARY_N_TAXA <- "n_taxa"
COLNAME_SIGNIF_SUMMARY_SIGNIF <- "significance"

#  Significance per sample.
COLNAME_N_SAMPLES_SIGNIF_RANK <- "rank"
COLNAME_N_SAMPLES_SIGNIF_SIGNIF <- "significance"
COLNAME_N_SAMPLES_SIGNIF_N_SAMPLES <- "n_samples"

#### Info/Warning/Error messages ####

LOAD_MPA_ERROR_NO_REPORTS_FOUND <- "LOG ERROR: No MPA-style reports were found in the directory provided. Please review your input:"
LOAD_MPA_WARNING_NO_SAMPLES_REMOVED <- "LOG WARNING: No samples were filtered out from the collated MPA-style reports table."
LOAD_MPA_INFO_SUCCESS <- "LOG INFO: MPA-style reports loaded successfully!"

LOAD_STD_ERROR_NO_REPORTS_FOUND <- "LOG ERROR: No standard reports were found in the directory provided. Please review your input:"
LOAD_STD_WARNING_NO_SAMPLES_REMOVED <- "LOG WARNING: No samples were filtered out from the collated standard reports table."
LOAD_STD_INFO_SUCCESS <- "LOG INFO: Standard reports loaded successfully!"

LOAD_SAMPLES_TO_REMOVE_WARNING <- "LOG WARNING: The following samples will not be included in the SPARKI analysis:"
LOAD_SAMPLES_TO_REMOVE_SUCCESS <- "LOG INFO: Samples-to-remove loaded successfully!"

CLI_ERROR_ARGUMENT_NOT_PROVIDED <- "LOG ERROR: Please review your input! The following argument is required but has not been provided:"
CLI_INFO_ARGUMENTS <- "LOG INFO: Running SPARKI with the following arguments:"
CLI_INFO_ARGUMENT_STD_REPORTS <- "LOG INFO:\t- Standard reports directory:"
CLI_INFO_ARGUMENT_MPA_REPORTS <- "LOG INFO:\t- MPA-style reports directory:"
CLI_INFO_ARGUMENT_ORGANISM <- "LOG INFO:\t- Organism set to:"
CLI_INFO_ARGUMENT_REFERENCE <- "LOG INFO:\t- Reference DB:"
CLI_INFO_ARGUMENT_METADATA <- "LOG INFO:\t- Metadata:"
CLI_INFO_ARGUMENT_SAMPLE_COL <- "LOG INFO:\t- Metadata sample column is:"
CLI_INFO_ARGUMENT_COLUMNS <- "LOG INFO:\t- Metadata columns are:"
CLI_INFO_ARGUMENT_OUTDIR <- "LOG INFO:\t- Output directory set to:"
CLI_INFO_ARGUMENT_PREFIX <- "LOG INFO:\t- Prefix set to: "
CLI_INFO_ARGUMENT_VERBOSE <- "LOG INFO:\t- Verbose?"
CLI_INFO_ARGUMENT_INC_EUKARYOTES <- "LOG INFO:\t- Include eukaryotes?"
CLI_INFO_ARGUMENT_INC_SAMPLE_NAMES <- "LOG INFO:\t- Include sample names?"
CLI_INFO_ARGUMENT_DOMAIN <- "LOG INFO:\t- Domain(s) of interest is(are):"
CLI_INFO_ARGUMENT_SAMPLES_TO_REMOVE <- "LOG INFO:\t- Samples to remove:"

LOAD_METADATA_INFO_SUCCESS <- "LOG INFO: Metadata loaded successfully."

CHECK_EMPTY_FILE_WARNING <- "LOG WARNING: An empty report was found! Carrying on...\nThe following sample will not be included in the SPARKI analysis:"

MAIN_WARNING_PREFIX_NOT_PROVIDED <- "LOG WARNING: No prefix has been provided. Carrying on..."
MAIN_WARNING_NO_SAMPLES_TO_REMOVE_PROVIDED <- "LOG WARNING: No list of samples to be removed has been provided. Carrying on..."
MAIN_WARNING_METADATA_NOT_PROVIDED <- "LOG WARNING: No metadata file has been provided. Carrying on..."
MAIN_WARNING_NO_METADATA_ADDED <- "LOG WARNING: No metadata was added to the reports. Carrying on..."
MAIN_ERROR_NO_COLUMNS_SPECIFIED <- "LOG ERROR: A metadata table has been provided, but no columns have been specified. Please review your input!"
MAIN_ERROR_NO_SAMPLE_COL_SPECIFIED <- "LOG ERROR: A metadata table has been provided, but no sample column has been specified. Please review your input!"
