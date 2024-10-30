###############################################
## FUNCTION FOR LOADING A REFERENCE DATABASE ##
#######################################################################################################

#' LOAD INFORMATION FROM KRAKEN2'S REFERENCE DATABASE
#' 
#' This function takes the path to an 'inspect.txt' file inside a Kraken2 reference database
#' and reads the table. The first 7 lines are skipped, as they correspond to header lines.
#' Column names are added before the reference database dataframe is returned.
#' 
#' @param reference_path Path to an 'inspect.txt' file inside a Kraken2 reference database.
#' @return A dataframe containing the information from the 'inspect.txt' file. 
#' @export
#' 
loadReference <- function(reference_path, n_header_lines = 7) {

    # Read Kraken2 reference file (inspect.txt).
    ref <- read.table(
        reference_path,
        row.names = NULL, 
        header = FALSE,
        sep = "\t",
        skip = n_header_lines, # Skip header lines.
        quote = "",
        stringsAsFactors = FALSE
    )
    # Add column names.
    colnames(ref) <- c(
        COLNAME_REF_DB_PCT_FRAG_CLADE, 
        COLNAME_REF_DB_MINIMISERS_CLADE,
        COLNAME_REF_DB_MINIMISERS_TAXON,
        COLNAME_REF_DB_RANK, 
        COLNAME_REF_DB_NCBI_ID, 
        COLNAME_REF_DB_TAXON
    )

    ref[, COLNAME_REF_DB_TAXON] <- gsub("^[[:space:]]\\s*(.*?)", "", ref[, COLNAME_REF_DB_TAXON], perl = TRUE)

    return(ref)
}

#############################################
## UTILITY FUNCTIONS FOR HANDLING METADATA ##
#######################################################################################################

#' LOAD METADATA TABLE
#' 
#' This function takes a path to a metadata file and reads the table.
#' and reads the table. The first 7 lines are skipped, as they correspond to header lines.
#' Column names are added before the reference database dataframe is returned.
#' 
#' @param mdata_path Path to a metadata table.
#' @return 
#' @export
#' 
loadMetadata <- function(mdata_path) {

    # Read metadata file.
    mdata <- readr::read_delim(mdata_path)

    return(mdata)
}

#' ADD COLUMNS FROM ONE DATAFRAME TO ANOTHER
#' 
#' This function takes as input a dataframe (df1) and adds metadata to it from another dataframe (df2). Both 
#' dataframes must have the same sample IDs present, as these IDs will be used as an "anchor" between the tables. 
#' Sample IDs can be specified as either row names or a standard column in the dataframes. The columns from df2
#' that should be added to df1 are specified by the user.
#' 
#' @param df1 Dataframe 1.
#' @param df1_sample_col Name of column in df1 that contains sample IDs; alternatively, this can be "rownames".
#' @param df2 Dataframe 2.
#' @param df2_sample_col Name of column in df2 that contains sample IDs; alternatively, this can be "rownames".
#' @param categories Categories of interest from df2 that should be added to df1.
#' @return An updated version of df1.
#' @export
#' 
addMetadata <- function(report, metadata, metadata_sample_col, metadata_columns) {

    # Check report format.
    report_colname_sample <- ifelse(is_mpa(report), COLNAME_MPA_SAMPLE, COLNAME_STD_SAMPLE)

    # Process metadata dataframe (tibble).
    metadata <- metadata |> 
        # Select metadata columns that will be added to the report.
        dplyr::select(metadata_sample_col, metadata_columns) |>
        # Rename column with sample IDs to keep its name consistent
        # with the report.
        dplyr::rename(!!report_colname_sample := metadata_sample_col)
    
    # Replace spaces (if any) with underscores.
    colnames(metadata) <- gsub(" ", "_", colnames(metadata))
    metadata_columns <- gsub(" ", "_", metadata_columns)

    # Add metadata columns to the report.
    report <- dplyr::left_join(
        report, metadata,
        by = dplyr::join_by(!!report_colname_sample)
    ) |>
        # Get names of columns that contain results (and not sample names / metadata) and
        # reorder the columns so that the metadata stays between the sample IDs and the
        # Kraken2 results.
        dplyr::relocate(metadata_columns, .after = !!report_colname_sample)

    return(report)
}

###################################################################
## UTILITY FUNCTIONS FOR READING, MERGING AND SUBSETTING REPORTS ##
#######################################################################################################

check_mpa_lines <- function(mpa_reports) {

    expected_ranks <- c("d__", "k__", "p__", "c__", "o__", "f__", "g__", "s__")

    updated_hierarchies <- mpa_reports |>
        # Create temporary column to store the original hierarchy column.
        dplyr::mutate(temp = !!as.name(COLNAME_MPA_TAXON_HIERARCHY)) |>
        # Split hierarchies across 
        tidyr::separate_rows(!!COLNAME_MPA_TAXON_HIERARCHY, sep = "\\|") |>
        dplyr::mutate(
            rank_ = substr(!!as.name(COLNAME_MPA_TAXON_HIERARCHY), 1, 3),
            taxon_ = substr(!!as.name(COLNAME_MPA_TAXON_HIERARCHY), 4, nchar(!!as.name(COLNAME_MPA_TAXON_HIERARCHY)))
        ) |>
        dplyr::select(!!as.name(COLNAME_MPA_SAMPLE), temp, rank_, taxon_) |>
        dplyr::group_by(!!as.name(COLNAME_MPA_SAMPLE), temp) |>
        tidyr::complete(rank_ = expected_ranks, fill = list(taxon_ = NA)) |>
        dplyr::arrange(factor(rank_, levels = expected_ranks), .by_group = TRUE) |> 
        dplyr::mutate(rank_taxon = paste0(rank_, taxon_)) |>
        dplyr::summarise(updated_hierarchy = paste(rank_taxon, collapse = "|"), .groups = "keep") |>
        dplyr::rename(!!COLNAME_MPA_TAXON_HIERARCHY := temp)

    mpa_reports <- dplyr::left_join( 
        mpa_reports, updated_hierarchies,
        by = dplyr::join_by(!!COLNAME_MPA_SAMPLE, !!COLNAME_MPA_TAXON_HIERARCHY)
    ) |>
        dplyr::select(!!COLNAME_MPA_SAMPLE, updated_hierarchy, !!COLNAME_MPA_N_FRAG_CLADE, !!COLNAME_MPA_RANK) |>
        dplyr::rename(!!COLNAME_MPA_TAXON_HIERARCHY := updated_hierarchy)

    return(mpa_reports)
}

#' LOAD MPA-STYLE REPORTS 
#' 
#' This function takes a path to a directory containing MPA-style reports, reading and processing
#' all reports into a single dataframe.
#' 
#' @param mpa_reports_dir Path to a directory containing MPA-style reports.
#' @param verbose Whether to output log messages.
#' @return A dataframe (tibble) with the content of all MPA-style reports from the specified directory.
#' @export
#'
load_MPAreports <- function(mpa_reports_dir, verbose = TRUE) {

    # Get paths to MPA-style reports in a specified directory.
    mpa_files <- fs::dir_ls(mpa_reports_dir, glob = "*.mpa$")

    # Check if directory really has any MPA-style reports...
    if (length(mpa_files) == 0) {
        stop(paste0("No MPA-style reports were found at ", mpa_reports_dir, ". Please review your input."))
    } 

    # Create vector with taxonomic rank names.
    taxonomy <- c("domain", "kingdom", "phylum", "class", "order", "family", "genus", "species")

    # Create a dataframe (tibble).
    mpa_reports <- readr::read_tsv(
        mpa_files, 
        col_names = c(COLNAME_MPA_TAXON_HIERARCHY, COLNAME_MPA_N_FRAG_CLADE), 
        id = get("COLNAME_MPA_SAMPLE")
    ) 

    # Add rank column.
    mpa_reports[, COLNAME_MPA_RANK] <- lapply(mpa_reports[, COLNAME_MPA_TAXON_HIERARCHY], addRank)

    # Double check each line to ensure all lines have complete cases (i.e. all ranks between domain and species).
    mpa_reports <- check_mpa_lines(mpa_reports)

    # Further process the dataframe.
    mpa_reports <- mpa_reports |>
        # Split rows by "|" into primary taxa.
        tidyr::separate(col = COLNAME_MPA_TAXON_HIERARCHY, into = taxonomy, sep = "\\|") |>
        # Cleanup the names in the taxonomy columns.
        dplyr::mutate(dplyr::across(taxonomy, ~ stringr::str_remove(.x, pattern = "[a-z]__"))) |>
        # Simplify sample IDs.
        dplyr::mutate(!!COLNAME_MPA_SAMPLE := stringr::str_remove(basename(sample), ".kraken.mpa")) |>
        # Convert character "NA"s to proper NA values.
        dplyr::mutate(dplyr::across(where(is.character), ~dplyr::na_if(., "NA"))) |>
        # Collect the rightmost non-NA item in each row.
        dplyr::mutate(
            !!COLNAME_MPA_TAXON_LEAF := dplyr::coalesce(species, genus, family, order, class, phylum, kingdom, domain),
            .before = "domain"
        ) |>
        # Replace underscores with spaces in taxon names.
        dplyr::mutate(
            !!COLNAME_MPA_TAXON_LEAF := stringr::str_replace_all(taxon_leaf, pattern = "_", replacement = " ")
        ) |>
        # Put rank column right after the taxon leaf column.
        dplyr::relocate(!!COLNAME_MPA_RANK, !!COLNAME_MPA_N_FRAG_CLADE, .after = !!COLNAME_MPA_TAXON_LEAF)

    if (verbose) cat("MPA-style reports loaded successfully.\n")
    
    return(mpa_reports) 
}

#' LOAD STANDARD REPORTS 
#' 
#' This function takes a path to a directory containing standard reports, reading and processing
#' all reports into a single dataframe.
#' 
#' @param mpa_reports_dir Path to a directory containing standard reports.
#' @param verbose Whether to output log messages.
#' @return A dataframe (tibble) with the content of all standard reports from the specified directory.
#' @export
#'
load_STDreports <- function(std_reports_dir, verbose = TRUE) {

    # Get paths to standard reports in a specified directory.
    std_files <- fs::dir_ls(std_reports_dir, glob = "*.kraken$")

    # Check if directory really has any standard reports...
    if (length(std_files) == 0) {
        stop(paste0("No standard reports were found at ", std_reports_dir, ". Please review your input."))
    }

    # Create a dataframe (tibble) and process.
    std_reports <- readr::read_tsv(
        std_files,
        col_names = c(COLNAME_STD_PCT_FRAG_CLADE, COLNAME_STD_N_FRAG_CLADE, COLNAME_STD_N_FRAG_TAXON, 
            COLNAME_STD_MINIMISERS, COLNAME_STD_UNIQ_MINIMISERS, COLNAME_STD_RANK, COLNAME_STD_NCBI_ID, 
            COLNAME_STD_TAXON),
        id = get("COLNAME_STD_SAMPLE")
    ) |>
        # Remove the subranks.
        dplyr::filter(!grepl(!!as.name(COLNAME_STD_RANK), pattern = "[0-9]")) |>
        # Simplify sample IDs.
        dplyr::mutate(!!COLNAME_STD_SAMPLE := stringr::str_remove(basename(sample), ".kraken"))

    if (verbose) cat("Standard reports loaded successfully.\n")

    return(std_reports) 
}

#' LOAD LIST OF SAMPLES TO BE EXCLUDED FROM SPARKI ANALYSIS
#' 
#' This function takes the path to a tab-delimited file, in which each line should
#' be a sample ID, and returns a list with these IDs. Note that the sample IDs need
#' to match those that are present in the merged reports dataframe.
#' 
#' @param filepath Path to a tab-delimited file in which each line is a sample ID.
#' @return A list containing the samples IDs from the input file.
#' @export
#'
loadSamplesToRemove <- function(filepath) {

    samples_to_remove <- readr::read_delim(filepath, col_names = "sample")
    return(samples_to_remove[["sample"]])
}
