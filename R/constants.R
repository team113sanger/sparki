##--------------##
## PACKAGE INFO ##
##--------------##

# CLI help text.
CLI_DESCRIPTION <- "SPARKI (Statistical Process Aimed at Robust Kraken2 Interpretation) is a framework to help collate, refine, visualise & interpret Kraken2 outputs."
CLI_PROGRAM_NAME <- "SPARKI"

##--------------##
## COLUMN NAMES ##
##--------------##

#### Main dataframes ####

# Dataframe with data from the reference database 'inspect.txt' file.
COLNAME_REF_DB_PCT_FRAG_CLADE <- "pct_fragments_clade"
COLNAME_REF_DB_MINIMISERS_CLADE <- "n_minimisers_clade"
COLNAME_REF_DB_MINIMISERS_TAXON <- "n_minimisers_taxon"
COLNAME_REF_DB_RANK <- "rank"
COLNAME_REF_DB_NCBI_ID <- "ncbi_id"
COLNAME_REF_DB_TAXON <- "taxon"

# Dataframe with collated standard reports.
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

# Dataframe with collated MPA-style reports.
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

#### Auxiliary dataframes ####

# Dataframe with the number of classified/unclassified reads per sample.
COLNAME_CLASSIF_SUMMARY_SAMPLE <- "sample"
COLNAME_CLASSIF_SUMMARY_READ_TYPE <- "type"
COLNAME_CLASSIF_SUMMARY_N_FRAG <- "n_fragments_clade"
COLNAME_CLASSIF_SUMMARY_LOG_N_FRAG <- "log10_n_fragments_clade"

# Dataframe with the proportion of classified reads per sample.
COLNAME_PROP_SAMPLE <- "sample"
COLNAME_PROP_CLASSIFIED <- "proportion_classified"

# Dataframe with the number of domain-level reads per sample.
COLNAME_DOMAIN_READS_SAMPLE <- "sample"
COLNAME_DOMAIN_READS_TAXON <- "taxon"
COLNAME_DOMAIN_READS_N_FRAG <- "n_fragments_clade"
COLNAME_DOMAIN_READS_LOG_N_FRAG <- "log10_n_fragments_clade"

# Dataframe with the significance summary.
COLNAME_SIGNIF_SUMMARY_SAMPLE <- "sample"
COLNAME_SIGNIF_SUMMARY_RANK <- "rank"
COLNAME_SIGNIF_SUMMARY_N_TAXA <- "n_taxa"
COLNAME_SIGNIF_SUMMARY_SIGNIF <- "significance"

# Dataframe with the significance per sample.
COLNAME_N_SAMPLES_SIGNIF_RANK <- "rank"
COLNAME_N_SAMPLES_SIGNIF_SIGNIF <- "significance"
COLNAME_N_SAMPLES_SIGNIF_N_SAMPLES <- "n_samples"

##-----------------##
## LOGGER MESSAGES ##
##-----------------##

#### load-data.R ####

# load_MPAreport()
LOAD_MPA_NO_REPORTS_FOUND <- "No MPA-style reports were found in the directory provided. Please review your input:"
LOAD_MPA_NO_SAMPLES_REMOVED <- "No samples were filtered out from the collated MPA-style reports table."
LOAD_MPA_SUCCESS <- "MPA-style reports loaded successfully!"

# load_STDreport()
LOAD_STD_NO_REPORTS_FOUND <- "No standard reports were found in the directory provided. Please review your input:"
LOAD_STD_NO_SAMPLES_REMOVED <- "No samples were filtered out from the collated standard reports table."
LOAD_STD_SUCCESS <- "Standard reports loaded successfully!"

# check_for_empty_files()
CHECK_EMPTY_FILE_WARNING <- "An empty report was found! Carrying on... The following sample will not be included in the SPARKI analysis:"

# loadMetadata()
LOAD_METADATA_SUCCESS <- "Metadata loaded successfully."

LOAD_SAMPLES_TO_REMOVE_WARNING <- "The following samples will not be included in the SPARKI analysis:"
LOAD_SAMPLES_TO_REMOVE_SUCCESS <- "Samples-to-remove loaded successfully!"

#### cli.R ####

# cli()
CLI_VERBOSITY_NOT_VALID <- "The verbosity argument provided is not valid, therefore the verbosity level will be set automatically."
CLI_ARGUMENT_NOT_PROVIDED <- "Please review your input! The following argument is required but has not been provided: "
CLI_ARGUMENTS <- "Running SPARKI with the following arguments: "
CLI_ARGUMENT_STD_REPORTS <- "Standard reports directory: "
CLI_ARGUMENT_MPA_REPORTS <- "MPA-style reports directory: "
CLI_ARGUMENT_ORGANISM <- "Organism set to: "
CLI_ARGUMENT_REFERENCE <- "Reference DB: "
CLI_ARGUMENT_METADATA <- "Metadata: "
CLI_ARGUMENT_SAMPLE_COL <- "Metadata sample column is: "
CLI_ARGUMENT_COLUMNS <- "Metadata columns are: "
CLI_ARGUMENT_OUTDIR <- "Output directory set to: "
CLI_ARGUMENT_PREFIX <- "Prefix set to: "
CLI_ARGUMENT_VERBOSE <- "Verbosity level set to: "
CLI_ARGUMENT_INC_EUKARYOTES <- "Include eukaryotes? "
CLI_ARGUMENT_INC_SAMPLE_NAMES <- "Include sample names? "
CLI_ARGUMENT_DOMAIN <- "Domain(s) of interest is(are): "
CLI_ARGUMENT_SAMPLES_TO_REMOVE <- "Samples to remove: "

#### main.R ####

# check_inputs()
MAIN_PREFIX_NOT_PROVIDED <- "No prefix has been provided. Carrying on..."
MAIN_NO_SAMPLES_TO_REMOVE_PROVIDED <- "No list of samples to be removed has been provided. Carrying on..."
MAIN_METADATA_NOT_PROVIDED <- "No metadata file has been provided. Carrying on..."
MAIN_NO_STD_PROVIDED <- "No standard reports directory has been specified. Please review your input!"
MAIN_NO_MPA_PROVIDED <- "No MPA-style reports directory has been specified. Please review your input!"
MAIN_NO_REF_PROVIDED <- "No reference database file has been specified. Please review your input!"
MAIN_NO_OUTDIR_PROVIDED <- "No output directory has been specified. Please review your input!"
MAIN_NO_DOMAIN_PROVIDED <- "No domain has been specified. Please review your input!"
MAIN_NO_COLUMNS_SPECIFIED <- "A metadata table has been provided, but no columns have been specified. Please review your input!"
MAIN_NO_SAMPLE_COL_SPECIFIED <- "A metadata table has been provided, but no sample column has been specified. Please review your input!"

# load_data()
MAIN_NO_METADATA_ADDED <- "No metadata was added to the reports. Carrying on..."

#### reports-merge.R ####

# mergeReports()
MERGE_REPORTS_SUCCESS <- "Standard and MPA-style reports merged successfully!"

#### exceptions.R ####

# check_directory()
EXCEPTIONS_CHECKDIR_DOES_NOT_EXIST <- "The directory provided does not exist. Please review your input! The directory path specified was: "
EXCEPTIONS_CHECKDIR_NOT_EMPTY <- "The directory provided is empty but it shouldn't be. Please review your input! The directory path specified was: "
EXCEPTIONS_CHECKDIR_EMPTY <- "The directory provided is not empty but it should be. Please review your input! The directory path specified was: "

# check_report_directory()
EXCEPTIONS_CHECK_REPORT_DIR_NO_MPA <- "The directory provided does not contain any MPA-style reports. Please review your input! The directory path specified was: "
EXCEPTIONS_CHECK_REPORT_DIR_NO_STD <- "The directory provided does not contain any standard reports. Please review your input! The directory path specified was: "
EXCEPTIONS_CHECK_REPORT_DIR_STD_IN_MPA <- "The MPA-style reports directory provided also contains standard reports. Please review your input! The directory path specified was: "
EXCEPTIONS_CHECK_REPORT_DIR_MPA_IN_STD <- "The standard reports directory provided also contains MPA-style reports. Please review your input! The directory path specified was: "
EXCEPTIONS_CHECK_REPORT_DIR_INVALID <- "The report format specified is not valid. Valid formats are 'std' or 'mpa', but got: "
EXCEPTIONS_CHECK_REPORT_DIR_BOTH_FORMATS <- "Both standard and MPA-style reports were found in the same directory, and this is not supported. The directory specified was: "

# check_file()
EXCEPTIONS_CHECKFILE_DOES_NOT_EXIST <- "The file provided does not exist. Please review your input! The file path specified was: "
EXCEPTIONS_CHECKFILE_EMPTY <- "The file provided is empty. Please review your input! The file path specified was: "

# check_columns()
EXCEPTIONS_CHECK_COLUMNS <- "The column provided does not exist in the metadata table provided. Please review your input! The column specified was: "

# check_prefix()
EXCEPTIONS_CHECK_PREFIX_NOT_CHARACTER <- "The prefix must be of class 'character', but got: "

# check_domain()
EXCEPTIONS_CHECK_DOMAIN_NOT_VALID <- "The domain specified is not valid. Please choose from 'Viruses', 'Bacteria', 'Archaea' or 'Eukaryota'. The value provided was: "

# check_organism()
EXCEPTIONS_CHECK_ORGANISM_STD <- "The organism specified was not found in the standard reports. Please review your input! The organism provided was: "
EXCEPTIONS_CHECK_ORGANISM_MPA <- "The organism specified was not found in the MPA-style reports. Please review your input! The organism provided was: "

#### reports-add-info.R ####

# addSampleSize()
ADD_INFO_SAMPLE_SIZE_SUCCESS <- "The sample sizes were calculated and added to the report."

# addMinimiserData()
ADD_INFO_MINIMISERS_SUCCESS <- "The minimiser info was retrieved and added to the report."
ADD_INFO_METADATA_SUCCESS <- "The metadata columns specified were added to the report."

#### plot-classification-summary.R ####

# plotClassificationSummary_violin()
PLOT_CLASS_SUMM_VIOLIN_SUCCESS <- "Created violin plot with classication summary."

# plotClassificationSummary_barplot()
PLOT_CLASS_SUMM_BARPLOT_SUCCESS <- "Created bar plot with classication summary."

# plotClassificationProportion()
PLOT_CLASS_PROPORTION_SUCCESS <- "Created violin plot with classication proportion."

##-------------------##
## OUTPUT FILE NAMES ##
##-------------------##

PLOT_CLASSIF_SUMMARY_VIOLIN <- "nReads_classified_vs_unclassified_absNumbers_violinPlot"
PLOT_CLASSIF_SUMMARY_BARPLOT <- "nReads_classified_vs_unclassified_proportion_perSample_barPlot"
PLOT_CLASSIF_PROPORTION <- "nReads_classifiedProportion_violinPlot"
