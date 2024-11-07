#' Helper function to get path to test data directory
#'
#' @returns Path to test data directory
get_test_mpa_report_dir <- function() {
  testthat::test_path("testdata", "mpa_reports")
}
