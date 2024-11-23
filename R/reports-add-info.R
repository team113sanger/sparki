#####################################################
##  UTILITY FUNCTIONS FOR ADDING COLUMNS TO REPORTS ##
#######################################################################################################

addRank <- function(taxon) {
  # Define rank prefixes and corresponding letters.
  rank_prefixes <- c("d__", "k__", "p__", "c__", "o__", "f__", "g__", "s__")
  rank_letters <- c("D", "K", "P", "C", "O", "F", "G", "S")

  # Create a single regex pattern.
  pattern <- paste0("(", paste(rank_prefixes, collapse = "|"), ")")

  # Extract all rank prefixes from each taxon string.
  all_matches <- stringr::str_extract_all(taxon, pattern)

  # Find the last (most specific) rank for each taxon.
  last_ranks <- sapply(all_matches, function(x) tail(x, 1))

  # Map rank prefixes to letters.
  rank <- rank_letters[match(last_ranks, rank_prefixes)]

  return(rank)
}

#' CALCULATE SAMPLE SIZES
#'
#' @param report A report.
#' @return An updated version of the input dataframe, with a new column containing sample sizes.
#' @export
#'
addSampleSize <- function(report, verbose) {
  #  If report is in MPA format...
  if (is_mpa(report)) {
    stop(paste0(
      "This function is not applicable to MPA-style reports. The purpose of this function ",
      "is to assess the total number of reads (classified + unclassified) assessed in each ",
      "sample (i.e., sample size), but the MPA-style reports do not contain information on ",
      "unclassified reads. If you would like to have the total number of reads in your MPA-style ",
      "report, please do:\n\nstd_reports <- add_nReads(std_reports)\n", "mpa_reports <- ",
      "transfer_nReads(mpa_reports, std_reports)\n\n"
    ))
  }

  subset <- report |>
    #  Select columns of interest.
    dplyr::select(!!COLNAME_STD_SAMPLE, !!COLNAME_STD_TAXON, !!COLNAME_STD_N_FRAG_CLADE) |>
    # Select rows of interest.
    dplyr::filter(!!as.name(COLNAME_STD_TAXON) %in% c("unclassified", "root")) |>
    #  Reformat dataframe into a more convenient shape.
    tidyr::pivot_wider(
      names_from = !!COLNAME_STD_TAXON,
      values_from = !!COLNAME_STD_N_FRAG_CLADE
    ) |>
    #  Rename columns.
    dplyr::rename(
      Unclassified = unclassified,
      Classified = root
    ) |>
    #  Create column with sample sizes (i.e. sum of unclasified and classified reads).
    dplyr::mutate(!!COLNAME_STD_SAMPLE_SIZE := rowSums(dplyr::across((where(is.numeric))))) |>
    #  Select final columns.
    dplyr::select(!!COLNAME_STD_SAMPLE, !!COLNAME_STD_SAMPLE_SIZE)

  #  Add sample sizes to report.
  report <- dplyr::left_join(report, subset, by = dplyr::join_by(!!COLNAME_STD_SAMPLE)) |>
    #  Reorder columns.
    dplyr::relocate(!!COLNAME_STD_SAMPLE_SIZE, .after = !!COLNAME_STD_SAMPLE)

  if (verbose) message("LOG INFO: Sample size successfully added to the report!")

  return(report)
}

addMinimiserData <- function(report, ref_db, verbose) {
  if (is_mpa(report)) stop("This function does not support MPA-style reports.")

  #  Select relevant columns.
  ref_db <- ref_db |>
    dplyr::select(
      !!COLNAME_REF_DB_NCBI_ID,
      !!COLNAME_REF_DB_MINIMISERS_TAXON,
      !!COLNAME_REF_DB_MINIMISERS_CLADE
    )

  report <- report |>
    #  Add DB info to the report.
    dplyr::left_join(ref_db, by = dplyr::join_by(!!COLNAME_STD_NCBI_ID)) |>
    # Rename columns.
    dplyr::rename(
      !!COLNAME_STD_DB_MINIMISERS_CLADE := !!COLNAME_REF_DB_MINIMISERS_CLADE,
      !!COLNAME_STD_DB_MINIMISERS_TAXON := !!COLNAME_REF_DB_MINIMISERS_TAXON
    )

  if (verbose) message("LOG INFO: Minimiser info successfully added to the report!")

  return(report)
}

add_nTaxaInRank <- function(report) {
  n_taxa_in_rank <- report |>
    #  Create column with sample+rank combinations.
    dplyr::mutate(
      temp = paste(!!as.name(COLNAME_STD_SAMPLE), !!as.name(COLNAME_STD_RANK), sep = "_")
    ) |>
    #  Count sample+rank occurrences (this will tell how many taxa are present in each rank
    #  for each sample).
    dplyr::count(temp)

  report <- report |>
    #  Create temporary column with sample+rank combinations.
    dplyr::mutate(
      temp = paste(!!as.name(COLNAME_STD_SAMPLE), !!as.name(COLNAME_STD_RANK), sep = "_")
    ) |>
    #  Merge dataframes.
    dplyr::left_join(n_taxa_in_rank, by = dplyr::join_by(temp)) |>
    #  Remove the temporary column.
    dplyr::select(!temp) |>
    #  Rename column.
    dplyr::rename(!!COLNAME_STD_N_TAXA_RANK := n)

  report[[COLNAME_STD_N_TAXA_RANK]] <- as.numeric(report[[COLNAME_STD_N_TAXA_RANK]])

  return(report)
}


#' ADD COLUMNS FROM ONE DATAFRAME TO ANOTHER
#'
#' This function takes as input a dataframe (df1) and adds metadata to it from another dataframe (df2). Both
#' dataframes must have the same sample IDs present, as these IDs will be used as an "anchor" between the tables.
#' Sample IDs can be specified as either row names or a standard column in the dataframes. The columns from df2
#' that should be added to df1 are specified by the user.
#'
#'  @param df1 Dataframe 1.
#'  @param df1_sample_col Name of column in df1 that contains sample IDs; alternatively, this can be "rownames".
#'  @param df2 Dataframe 2.
#'  @param df2_sample_col Name of column in df2 that contains sample IDs; alternatively, this can be "rownames".
#' @param categories Categories of interest from df2 that should be added to df1.
#' @return An updated version of df1.
#' @export
#'
addMetadata <- function(report, metadata, metadata_sample_col, metadata_columns, verbose) {
  #  Check report format.
  report_colname_sample <- ifelse(is_mpa(report), COLNAME_MPA_SAMPLE, COLNAME_STD_SAMPLE)

  #  Process metadata dataframe (tibble).
  metadata <- metadata |>
    # Select metadata columns that will be added to the report.
    dplyr::select(dplyr::all_of(metadata_sample_col), dplyr::all_of(metadata_columns)) |>
    # Rename column with sample IDs to keep its name consistent
    #  with the report.
    dplyr::rename(!!report_colname_sample := dplyr::all_of(metadata_sample_col))

  #  Replace spaces (if any) with underscores.
  colnames(metadata) <- gsub(" ", "_", colnames(metadata))
  metadata_columns <- gsub(" ", "_", metadata_columns)

  #  Add metadata columns to the report.
  report <- dplyr::left_join(
    report, metadata,
    by = dplyr::join_by(!!report_colname_sample)
  ) |>
    # Get names of columns that contain results (and not sample names / metadata) and
    #  reorder the columns so that the metadata stays between the sample IDs and the
    # Kraken2 results.
    dplyr::relocate(dplyr::all_of(metadata_columns), .after = !!report_colname_sample)

  if (verbose) message("LOG INFO: Metadata successfully added to the report!")

  return(report)
}
