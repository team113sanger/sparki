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

    message("LOG INFO: Reference database info loaded successfully!")

    return(ref)
}

#' LOAD METADATA TABLE
#'
#' This function takes the path to a metadata file and reads the metadata table using readr::read_delim().
#'
#' @param metadata Path to a metadata table.
#' @return A dataframe containing the metadata.
#' @export
#'
#' @examples
#' loadMetadata("metadata.csv")
#' loadMetadata(metadata = "metadata.csv")
#'
loadMetadata <- function(metadata, verbose) {

    # Read metadata file.
    mdata <- readr::read_delim(metadata)

    if (verbose == TRUE) {
        message("LOG INFO: Metadata loaded successfully.")
    }

    return(mdata)
}

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

check_for_empty_files <- function(file_list, verbose) {

    # Check if files are empty...
    for (file in file_list) {

        is_empty <- (file.size(file) == 0L)

        # If the file is empty...
        if (is_empty) {

            if (verbose) {
                message(
                    "LOG WARNING: The file ", file, " is empty, so this sample ",
                    "will not be included in the SPARKI analysis."
                )
            }

            # Remove given file path from the list of paths.
            file_list <- file_list[(file_list != file)]
        }
    }

    return(file_list)
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
load_MPAreports <- function(mpa_reports_dir, samples_to_remove, verbose) {

    # Get paths to MPA-style reports in a specified directory.
    mpa_files <- fs::dir_ls(mpa_reports_dir, glob = "*.mpa$")

    # Check if directory really has any MPA-style reports...
    if (length(mpa_files) == 0) {
        stop("No MPA-style reports were found at ", mpa_reports_dir, ". Please review your input.")
    }

    mpa_files <- check_for_empty_files(mpa_files, verbose = verbose)

    # Create vector with taxonomic rank names.
    taxonomy <- c("domain", "kingdom", "phylum", "class", "order", "family", "genus", "species")

    # Create a dataframe (tibble).
    mpa_reports <- readr::read_tsv(
        mpa_files,
        col_names = c(COLNAME_MPA_TAXON_HIERARCHY, COLNAME_MPA_N_FRAG_CLADE),
        id = get("COLNAME_MPA_SAMPLE"),
        show_col_types = FALSE # Supressing messages about column types when the dataframe is created.
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
        dplyr::mutate(dplyr::across(dplyr::all_of(taxonomy), ~ stringr::str_remove(.x, pattern = "[a-z]__"))) |>
        # Simplify sample IDs.
        dplyr::mutate(!!COLNAME_MPA_SAMPLE := stringr::str_remove(basename(sample), ".kraken.mpa")) |>
        # Convert character "NA"s to proper NA values.
        dplyr::mutate(dplyr::across(where(is.character), ~dplyr::na_if(., "NA"))) |>
        # Collect the rightmost non-NA item in each row.
        dplyr::mutate(
            !!COLNAME_MPA_TAXON_LEAF := dplyr::coalesce(
                !!as.name(COLNAME_MPA_SPECIES),
                !!as.name(COLNAME_MPA_GENUS),
                !!as.name(COLNAME_MPA_FAMILY),
                !!as.name(COLNAME_MPA_ORDER),
                !!as.name(COLNAME_MPA_CLASS),
                !!as.name(COLNAME_MPA_PHYLUM),
                !!as.name(COLNAME_MPA_KINGDOM),
                !!as.name(COLNAME_MPA_DOMAIN)
            ),
            .before = "domain"
        ) |>
        # Replace underscores with spaces in taxon names.
        dplyr::mutate(
            !!COLNAME_MPA_TAXON_LEAF := stringr::str_replace_all(taxon_leaf, pattern = "_", replacement = " ")
        ) |>
        # Put rank column right after the taxon leaf column.
        dplyr::relocate(!!COLNAME_MPA_RANK, !!COLNAME_MPA_N_FRAG_CLADE, .after = !!COLNAME_MPA_TAXON_LEAF)

    if (!missing(samples_to_remove)) {
        mpa_reports <- mpa_reports |> dplyr::filter(!(!!as.name(COLNAME_MPA_SAMPLE) %in% samples_to_remove))
    } else {
        if (verbose == TRUE) {
            message("LOG WARNING: No samples were filtered out from the collated MPA-style reports table.")
        }
    }

    if (verbose == TRUE) message("LOG INFO: MPA-style reports loaded successfully!")

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
load_STDreports <- function(std_reports_dir, samples_to_remove, verbose = FALSE) {

    # Get paths to standard reports in a specified directory.
    std_files <- fs::dir_ls(std_reports_dir, glob = "*.kraken$")

    # Check if directory really has any standard reports...
    if (length(std_files) == 0) {
        stop("No standard reports were found at ", std_reports_dir, ". Please review your input.")
    }

    std_files <- check_for_empty_files(std_files, verbose = verbose)

    # Create a dataframe (tibble) and process.
    std_reports <- readr::read_tsv(
        std_files,
        col_names = c(COLNAME_STD_PCT_FRAG_CLADE, COLNAME_STD_N_FRAG_CLADE, COLNAME_STD_N_FRAG_TAXON,
            COLNAME_STD_MINIMISERS, COLNAME_STD_UNIQ_MINIMISERS, COLNAME_STD_RANK, COLNAME_STD_NCBI_ID,
            COLNAME_STD_TAXON),
        id = get("COLNAME_STD_SAMPLE"),
        show_col_types = FALSE # Supressing messages about column types when the dataframe is created.
    ) |>
        # Remove the subranks.
        dplyr::filter(!grepl(!!as.name(COLNAME_STD_RANK), pattern = "[0-9]")) |>
        # Simplify sample IDs.
        dplyr::mutate(!!COLNAME_STD_SAMPLE := stringr::str_remove(basename(sample), ".kraken"))

    if (!missing(samples_to_remove)) {
        std_reports <- std_reports |> dplyr::filter(!(!!as.name(COLNAME_STD_SAMPLE) %in% samples_to_remove))
    } else {
        if (verbose) message("LOG WARNING: No samples were filtered out from the collated standard reports table.")
    }

    if (verbose) message("LOG INFO: Standard reports loaded successfully!")

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
loadSamplesToRemove <- function(filepath, verbose) {

    samples_to_remove <- readr::read_table(filepath, col_names = "sample")

    message(
        "LOG WARNING: Samples ", paste(samples_to_remove[["sample"]], collapse = ","),
        " will be removed from the SPARKI analysis."
    )

    if (verbose) message("LOG INFO: Samples-to-remove loaded successfully.")

    return(samples_to_remove[["sample"]])
}
