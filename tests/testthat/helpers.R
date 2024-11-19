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
  return(testthat::test_path("testdata", "mpa_reports", paste0("case_with_", case_type, "_reports")))
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
  return(testthat::test_path("testdata", "std_reports", paste0("case_with_", case_type, "_reports")))
}

#' Helper function to get path to inspect.txt file in Kraken2's reference database
#'
#' @returns Path to Kraken2 reference database
get_test_reference <- function() {
  return(testthat::test_path("testdata", "inspect.txt"))
}

#' Helper function to get path to a metadata file
#'
#' @returns Path to metadata file
get_test_metadata <- function() {
  return(testthat::test_path("testdata", "metadata.txt"))
}

#' Helper function to get path to a samples-to-remove file
#'
#' @returns Path to a file with samples to be removed from the SPARKI analysis
get_test_samples_to_remove <- function() {
  return(testthat::test_path("testdata", "samples_to_remove.txt"))
}

#' Create a temporary directory that auto-cleans after test
#'
#' @description
#' Creates a temporary directory for test files and automatically cleans it up
#' when the enclosing test or function exits, similar to pytest's tmp_path fixture.
#'
#' @param env The environment to bind cleanup behavior to. Defaults to parent.frame()
#'
#' @return Path to the created temporary directory as a character string
#'
#' @examples
#' test_that("files are processed correctly", {
#'   tmp_dir <- local_tmp_dir()
#'   writeLines("test", file.path(tmp_dir, "test.txt"))
#'   # Directory and contents are automatically cleaned up after test
#' })
#'
get_local_tmp_dir <- function(env = parent.frame()) {
  tmp_dir <- file.path(tempdir(), paste0("test-", format(Sys.time(), "%Y%m%d-%H%M%S-"), sample(1000, 1)))
  dir.create(tmp_dir, recursive = TRUE)

  withr::defer({
    unlink(tmp_dir, recursive = TRUE)
  }, env)

  return(tmp_dir)
}

# get_test_std_report_dir <- function(case_type) {
#   if (!(case_type %in% c("usual", "empty"))) {
#     stop(paste0(
#       case_type, " is not a valid case type for the function get_test_std_report_dir! ",
#       "Please choose from 'usual' or 'empty'."
#     ))
#   }
#   testthat::test_path("testdata", "std_reports", paste0("case_with_", case_type, "_reports"))
# }
