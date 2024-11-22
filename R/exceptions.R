#' CARRY OUT CHECKS ON A GIVEN DIRECTORY
#'
#' This function takes the path to a given directory and performs a few checks:
#' 1) First, it checks if the directory exists.
#' 2) Second, it checks if the directory is empty.
#' 3) Last, it ensures that the last character in the path is a slash symbol ("/").
#'
#' @param dirpath Path to a given directory.
#'
#' @seealso [check_report_directory()] for another function that carries out
#'  checks based on a given directory path.
#'
#' @examples
#' check_directory(dirpath = "path/to/directory")
#'
#' @return Checked path to the given directory.
#'
check_directory <- function(dirpath) {
  # Check that the directory exists.
  if (!(dir.exists(dirpath))) {
    stop(EXCEPTIONS_CHECKDIR_DOES_NOT_EXIST, dirpath)

  # Check that the directory is not empty.
  } else if (length(list.files(dirpath)) == 0) {
    stop(EXCEPTIONS_CHECKDIR_EMPTY, dirpath)

  # Ensure that the directory path has a slash at the end.
  } else if (substr(dirpath, nchar(dirpath), nchar(dirpath)) != "/") {
    dirpath <- paste0(dirpath, "/")
  }

  return(dirpath)
}

#' CARRY OUT CHECKS ON DIRECTORIES CONTAINING SPECIFICALLY STANDARD OR MPA-STYLE REPORTS
#'
#' This function takes the path to a directory containing either standard or
#' MPA-style reports and performs a few checks:
#' - If the directory is expected to contain MPA-style reports, the function will check
#' whether MPA-style reports are indeed present.
#' - If the directory is expected to contain standard reports, the function will check
#' whether standard reports are indeed present.
#' - The function also ensures that the directory provided has only one report format,
#' throwing an error otherwise.
#'
#' @param dirpath Path to a given directory.
#' @param report_format One of 'std' (standard) or 'mpa' (MPA-style), indicating
#'  the report format that is expected in the directory specified.
#'
#' @seealso [check_directory()] for another function that carries out
#'  checks based on a given directory path.
#'
#' @examples
#' check_report_directory(dirpath = "path/to/directory", report_format = "std")
#'
#' @return Checked path to the given directory.
#'
check_report_directory <- function(dirpath, report_format) {
  # Check that the directory exists and is not empty.
  dirpath <- check_directory(dirpath)

  # Get the number of standard or MPA-style reports in the directory.
  n_mpa_reports <- length(list.files(dirpath, pattern = "mpa$"))
  n_std_reports <- length(list.files(dirpath, pattern = "kraken$"))

  # Check that the directory contains MPA reports if it is supposed to.
  if (report_format == "mpa" && n_mpa_reports == 0) {
    stop(EXCEPTIONS_CHECK_REPORT_DIR_NO_MPA, dirpath)

  # Check that the directory contains standard reports if it is supposed to.
  } else if (report_format == "std" && n_std_reports == 0) {
    stop(EXCEPTIONS_CHECK_REPORT_DIR_NO_STD, dirpath)

  # Ensure that the report format specified is valid.
  } else if (!(report_format %in% c("std", "mpa"))) {
    stop(EXCEPTIONS_CHECK_REPORT_DIR_INVALID, report_format)

  # Check that the directory does not contain both report formats.
  } else if (n_mpa_reports > 0 && n_std_reports > 0) {
    stop(EXCEPTIONS_CHECK_REPORT_DIR_BOTH_FORMATS, dirpath)
  }

  return(dirpath)
}

#' CARRY OUT CHECKS ON A GIVEN FILE
#'
#' This function takes the path to a given file and performs a few checks:
#' 1) First, it checks if the file exists.
#' 2) Second, it checks if the file is empty.
#'
#' @param file_path Path to a given file.
#'
#' @examples
#' check_file(file_path = "path/to/file")
#'
check_file <- function(file_path) {
  if (!(file.exists(file_path))) stop(EXCEPTIONS_CHECKFILE_DOES_NOT_EXIST, file_path)
  else if (file.size(file_path) == 0L) stop(EXCEPTIONS_CHECKFILE_EMPTY, file_path)
}

#' CHECK IF A SET OF COLUMN NAMES IS PRESENT IN A GIVEN DATAFRAME
#'
#' This function takes a dataframe and a vector containing column names that
#' are expected to be found in the given dataframe; then it checks if each
#' column is indeed present, throwing an error otherwise.
#'
#' @param df Dataframe where given column names will be looked for.
#' @param columns A vector with column names that are expected to be present
#'  in the given dataframe.
#'
#' @examples
#' check_columns(df = df, columns = c("col1", "col2", "col3"))
#'
check_columns <- function(df, columns) {
  for (column in columns) {
    if (!(column %in% colnames(df))) stop(EXCEPTIONS_CHECK_COLUMNS, column)
  }
}

#' CARRY OUT CHECKS ON A PREFIX
#'
#' This function takes a prefix and performs a few checks:
#' 1) First, it checks if the prefix corresponds to the appropriate class ("character").
#' 2) Second, it ensures that there is an underscore symbol ("_") at the end of the prefix.
#'
#' @param prefix A given prefix that will later be added to output file names.
#' @return Checked prefix.
#'
#' @examples
#' check_prefix(prefix = "prefix")
#'
check_prefix <- function(prefix) {
  if (!(is.character(prefix))) stop(EXCEPTIONS_CHECK_PREFIX_NOT_CHARACTER, class(prefix))
  else if (substr(prefix, nchar(prefix), nchar(prefix)) != "_") prefix <- paste0(prefix, "_")

  return(prefix)
}

#' CHECK IF THE DOMAIN SPECIFIED IS VALID
#'
#' This function takes a domain and checks if it is valid, i.e. one of "Viruses",
#' "Bacteria", "Archaea", or "Eukaryota".
#'
#' @param domain A domain ("Viruses", "Bacteria", "Archaea", or "Eukaryota").
#'
#' @examples
#' check_domain(domain = "Viruses")
#'
check_domain <- function(domain) {
  if (!(domain %in% c("Viruses", "Bacteria", "Archaea", "Eukaryota"))) {
    stop(EXCEPTIONS_CHECK_DOMAIN_NOT_VALID, domain)
  }
}
