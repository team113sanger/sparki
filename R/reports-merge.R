#' Merge standard and MPA-style reports
#'
#' This function takes a standard and an MPA-style dataframe, loaded with load_STDreports()
#' and load_MPAreports() respectively, and merges the latter into the former.
#'
#' @param std_reports A dataframe with standard reports, loaded with load_STDreports().
#' @param mpa_reports A dataframe with MPA-style reports, loaded with load_MPAreports().
#'
#' @return An updated version of the dataframe with standard reports, now containing the
#' information from the dataframe with MPA-style reports.
#'
#' @export
#'
mergeReports <- function(std_reports, mpa_reports) {
  mpa_reports <- mpa_reports |>
    # Rename column for consistency with the standard report.
    dplyr::rename(!!COLNAME_STD_TAXON := !!COLNAME_MPA_TAXON_LEAF) |>
    # Drop columns that are already present in the standard report.
    dplyr::select(!c(!!COLNAME_MPA_N_FRAG_CLADE, !!COLNAME_MPA_RANK))

  # The primary keys that will be used for joining are "taxon leaf"/"name" and "sample_id"
  merged_reports <- dplyr::left_join(
    std_reports, mpa_reports,
    by = dplyr::join_by(!!COLNAME_STD_TAXON, !!COLNAME_STD_SAMPLE)
  )

  logger::log_success(MERGE_REPORTS_SUCCESS)

  return(merged_reports)
}
