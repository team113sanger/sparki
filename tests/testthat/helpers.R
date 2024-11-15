#' Helper function to get path to test data directory (MPA-style reports)
#'
#' @returns Path to test data directory
get_test_mpa_report_dir <- function(case_type) {
  if (!(case_type %in% c("usual", "empty"))) {
    stop(paste0(
      case_type, " is not a valid case type for the function get_test_mpa_report_dir! ",
      "Please choose from 'usual' or 'empty'."
    ))
  }
  testthat::test_path("testdata", "mpa_reports", paste0("case_with_", case_type, "_reports"))
}

#' Helper function to get path to test data directory (standard reports)
#'
#' @returns Path to test data directory
get_test_std_report_dir <- function(case_type) {
  if (!(case_type %in% c("usual", "empty"))) {
    stop(paste0(
      case_type, " is not a valid case type for the function get_test_std_report_dir! ",
      "Please choose from 'usual' or 'empty'."
    ))
  }
  testthat::test_path("testdata", "std_reports", paste0("case_with_", case_type, "_reports"))
}

#' Helper function to get path to inspect.txt file in Kraken2's reference database
#'
#' @returns Path to Kraken2 reference database
# get_test_std_report_dir <- function(case_type) {
#   if (!(case_type %in% c("usual", "empty"))) {
#     stop(paste0(
#       case_type, " is not a valid case type for the function get_test_std_report_dir! ",
#       "Please choose from 'usual' or 'empty'."
#     ))
#   }
#   testthat::test_path("testdata", "std_reports", paste0("case_with_", case_type, "_reports"))
# }
