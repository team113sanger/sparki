
#' Load MPA-style reports
#'
#' This function takes a path to a directory containing MPA-style reports, reading and collating
#' all sample-level reports into a single dataframe.
#'
#' @param mpa_reports_dir Path to a directory containing MPA-style reports.
#' @param samples_to_remove A list of samples to be removed, loaded with loadSamplesToRemove().
#'
#' @return A dataframe (tibble) with the content of all MPA-style reports from the specified directory.
#'
#' @export
#'
load_MPAreports <- function(mpa_reports_dir, samples_to_remove) {
  # Get paths to MPA-style reports in a specified directory.
  mpa_files <- fs::dir_ls(mpa_reports_dir, glob = "*.mpa$")

  #  Check if directory really has any MPA-style reports...
  if (length(mpa_files) == 0) stop(LOAD_MPA_NO_REPORTS_FOUND, mpa_reports_dir)

  mpa_files <- check_for_empty_files(mpa_files)

  #  Create vector with taxonomic rank names.
  taxonomy <- c("domain", "kingdom", "phylum", "class", "order", "family", "genus", "species")

  # Create a dataframe (tibble).
  mpa_reports <- readr::read_tsv(
    mpa_files,
    col_names = c(COLNAME_MPA_TAXON_HIERARCHY, COLNAME_MPA_N_FRAG_CLADE),
    id = get("COLNAME_MPA_SAMPLE"),
    col_types = "cic",
    show_col_types = FALSE # Supressing messages about column types when the dataframe is created.
  )

  #  Add rank column.
  mpa_reports[, COLNAME_MPA_RANK] <- lapply(mpa_reports[, COLNAME_MPA_TAXON_HIERARCHY], addRank)

  #  Double check each line to ensure all lines have complete cases (i.e. all ranks between domain and species).
  mpa_reports <- check_mpa_lines(mpa_reports)

  # Further process the dataframe.
  mpa_reports <- mpa_reports |>
    # Split rows by "|" into primary taxa.
    tidyr::separate(col = COLNAME_MPA_TAXON_HIERARCHY, into = taxonomy, sep = "\\|") |>
    # Cleanup the names in the taxonomy columns.
    dplyr::mutate(dplyr::across(dplyr::all_of(taxonomy), ~ stringr::str_remove(.x, pattern = "[a-z]__"))) |>
    # Simplify sample IDs.
    dplyr::mutate(!!COLNAME_MPA_SAMPLE := stringr::str_remove(basename(sample), ".kraken.mpa")) |>
    #  Convert character "NA"s to proper NA values.
    dplyr::mutate(dplyr::across(where(is.character), ~ dplyr::na_if(., "NA"))) |>
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
    #  Put rank column right after the taxon leaf column.
    dplyr::relocate(!!COLNAME_MPA_RANK, !!COLNAME_MPA_N_FRAG_CLADE, .after = !!COLNAME_MPA_TAXON_LEAF)

  if (!missing(samples_to_remove)) {
    mpa_reports <- mpa_reports |> dplyr::filter(!(!!as.name(COLNAME_MPA_SAMPLE) %in% samples_to_remove))
  } else {
    logger::log_info(LOAD_MPA_NO_SAMPLES_REMOVED)
  }

  logger::log_info(LOAD_MPA_SUCCESS)

  return(mpa_reports)
}

#' Load standard reports
#'
#' This function takes a path to a directory containing standard reports, reading and collating
#' all sample-level reports into a single dataframe.
#'
#' @param std_reports_dir Path to a directory containing standard reports.
#' @param samples_to_remove A list of samples to be removed, loaded with loadSamplesToRemove().
#'
#' @return A dataframe (tibble) with the content of all standard reports from the specified directory.
#'
#' @export
#'
load_STDreports <- function(std_reports_dir, samples_to_remove) {
  #  Get paths to standard reports in a specified directory.
  std_files <- fs::dir_ls(std_reports_dir, glob = "*.kraken$")

  #  Check if directory really has any standard reports...
  if (length(std_files) == 0) stop(LOAD_STD_NO_REPORTS_FOUND, std_reports_dir)

  std_files <- check_for_empty_files(std_files)

  # Create a dataframe (tibble) and process.
  std_reports <- readr::read_tsv(
    std_files,
    col_names = c(
      COLNAME_STD_PCT_FRAG_CLADE, COLNAME_STD_N_FRAG_CLADE, COLNAME_STD_N_FRAG_TAXON,
      COLNAME_STD_MINIMISERS, COLNAME_STD_UNIQ_MINIMISERS, COLNAME_STD_RANK, COLNAME_STD_NCBI_ID,
      COLNAME_STD_TAXON
    ),
    id = get("COLNAME_STD_SAMPLE"),
    col_types = "diiiicccc",
    show_col_types = FALSE # Supressing messages about column types when the dataframe is created.
  ) |>
    # Remove the subranks.
    dplyr::filter(!grepl(!!as.name(COLNAME_STD_RANK), pattern = "[0-9]")) |>
    # Simplify sample IDs.
    dplyr::mutate(!!COLNAME_STD_SAMPLE := stringr::str_remove(basename(sample), ".kraken"))

  if (!missing(samples_to_remove)) {
    std_reports <- std_reports |> dplyr::filter(!(!!as.name(COLNAME_STD_SAMPLE) %in% samples_to_remove))
  } else {
    logger::log_info(LOAD_STD_NO_SAMPLES_REMOVED)
  }

  logger::log_info(LOAD_STD_SUCCESS)

  return(std_reports)
}

#' Load list of samples to be removed from the SPARKI analysis
#'
#' This function takes the path to a tab-delimited file, in which each line should
#' be a sample ID, and returns a list with these IDs. Note that the sample IDs need
#' to match those that are present in the reports to be loaded with load_STDreports()
#' and load_MPAreports().
#'
#' @param filepath Path to a tab-delimited file in which each line is a sample ID.
#'
#' @return A list containing the samples IDs from the input file.
#'
#' @export
#'
loadSamplesToRemove <- function(filepath) {
  samples_to_remove <- readr::read_table(
    filepath,
    col_names = "sample",
    col_types = "c",
    show_col_types = FALSE # Supressing messages about column types when the dataframe is created.
  )

  logger::log_info(paste0(LOAD_SAMPLES_TO_REMOVE_WARNING, " ", paste(samples_to_remove[["sample"]], collapse = ",")))
  logger::log_info(LOAD_SAMPLES_TO_REMOVE_SUCCESS)

  return(samples_to_remove[["sample"]])
}

#' Load the 'inspect.txt' file from a Kraken2 reference database
#'
#' This function takes the path to an 'inspect.txt' file inside a Kraken2 reference database
#' and reads the table. The first 7 lines are skipped, as they correspond to header lines.
#' Column names are added before the reference database dataframe is returned.
#'
#' @param reference_path Path to an 'inspect.txt' file inside a Kraken2 reference database.
#'
#' @return A dataframe containing the information from the 'inspect.txt' file.
#'
#' @export
#'
loadReference <- function(reference_path) {
  # Read Kraken2 reference file (inspect.txt).
  ref <- utils::read.table(
    reference_path,
    row.names = NULL,
    header = FALSE,
    sep = "\t",
    skip = 7, #  Skip header lines.
    quote = "",
    stringsAsFactors = FALSE
  )

  #  Add column names.
  colnames(ref) <- c(
    COLNAME_REF_DB_PCT_FRAG_CLADE,
    COLNAME_REF_DB_MINIMISERS_CLADE,
    COLNAME_REF_DB_MINIMISERS_TAXON,
    COLNAME_REF_DB_RANK,
    COLNAME_REF_DB_NCBI_ID,
    COLNAME_REF_DB_TAXON
  )

  ref[, COLNAME_REF_DB_TAXON] <- gsub("^[[:space:]]\\s*(.*?)", "", ref[, COLNAME_REF_DB_TAXON], perl = TRUE)
  ref[, COLNAME_REF_DB_NCBI_ID] <- as.character(ref[, COLNAME_REF_DB_NCBI_ID])

  logger::log_info("Reference database info loaded successfully!")

  return(ref)
}

#' Load a metadata table
#'
#' This function takes the path to a metadata file and reads the metadata table using
#' readr::read_delim().
#'
#' @param metadata Path to a metadata table.
#'
#' @return A dataframe containing the metadata.
#'
#' @export
#'
loadMetadata <- function(metadata) {
  #  Read metadata file.
  mdata <- readr::read_delim(
    metadata,
    show_col_types = FALSE # Supressing messages about column types when the dataframe is created.
  )

  logger::log_info(LOAD_METADATA_SUCCESS)

  return(mdata)
}

#' Ensure taxon hierarchies in an MPA-style report are "complete"
#'
#' @param mpa_reports A report loaded by load_MPAreports().
#'
#' @return An updated version of the report provided, with all taxon
#'  hierarchies being "complete".
#'
check_mpa_lines <- function(mpa_reports) {
  expected_ranks <- c("d__", "k__", "p__", "c__", "o__", "f__", "g__", "s__")

  updated_hierarchies <- mpa_reports |>
    #  Create temporary column to store the original hierarchy column.
    dplyr::mutate(temp = !!as.name(COLNAME_MPA_TAXON_HIERARCHY)) |>
    #  Split hierarchies across
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

#' Exclude empty reports from the analysis
#'
#' @param file_list A list containing paths to standard or MPA-style reports.
#'
#' @return An updated version of file list, without empty reports.
#'
check_for_empty_files <- function(file_list) {
  #  Check if files are empty...
  for (file in file_list) {
    is_empty <- (file.size(file) == 0L)

    #  If the file is empty...
    if (is_empty) {
      # Remove given file path from the list of paths.
      logger::log_info(paste0(CHECK_EMPTY_FILE_WARNING, glue::double_quote(file)))
      file_list <- file_list[(file_list != file)]
    }
  }

  return(file_list)
}
