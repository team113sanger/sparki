#' Carry out checks on a given directory
#'
#' This function takes the path to a given directory and performs a few checks:
#' - It checks if the directory exists.
#' - It checks if the directory is empty or not depending on the expectation.
#' - It ensures that the last character in the path is a slash symbol ("/").
#'
#' @param dirpath Path to a given directory.
#' @param expectation One of "empty" or "not_empty".
#'
#' @seealso [check_report_directory()] for another function that carries out
#'  checks based on a given directory path.
#'
#' @return Checked path to the given directory.
#'
check_directory <- function(dirpath, expectation) {
  # Check that the directory exists.
  if (!(dir.exists(dirpath))) {
    error_message <- paste0(EXCEPTIONS_CHECKDIR_DOES_NOT_EXIST, glue::double_quote(dirpath))
    logger::log_fatal(error_message)
    stop(error_message)

  # If the directory is not supposed to be empty, check that the directory is not empty.
  } else if (expectation == "not_empty" && length(list.files(dirpath)) == 0) {
    error_message <- paste0(EXCEPTIONS_CHECKDIR_NOT_EMPTY, glue::double_quote(dirpath))
    logger::log_fatal(error_message)
    stop(error_message)

  # If the directory is supposed to be empty, check that the directory is empty.
  } else if (expectation == "empty" && length(list.files(dirpath)) > 0) {
    error_message <- paste0(EXCEPTIONS_CHECKDIR_EMPTY, glue::double_quote(dirpath))
    logger::log_fatal(error_message)
    stop(error_message)

  # Ensure that the directory path has a slash at the end.
  } else if (substr(dirpath, nchar(dirpath), nchar(dirpath)) != "/") {
    dirpath <- paste0(dirpath, "/")
  }

  return(dirpath)
}

#' Carry out checks on directories that contain standard or MPA-style reports
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
#' @return Checked path to the given directory.
#'
check_report_directory <- function(dirpath, report_format) {
  # Check that the directory exists and is not empty.
  dirpath <- check_directory(dirpath, expectation = "not_empty")

  # Get the number of standard or MPA-style reports in the directory.
  n_mpa_reports <- length(list.files(dirpath, pattern = "mpa$"))
  n_std_reports <- length(list.files(dirpath, pattern = "kraken$"))

  # Check that the directory contains MPA reports if it is supposed to.
  if (report_format == "mpa" && n_mpa_reports == 0) {
    error_message <- paste0(EXCEPTIONS_CHECK_REPORT_DIR_NO_MPA, glue::double_quote(dirpath))
    logger::log_fatal(error_message)
    stop(error_message)

  # Check that the MPA-style reports directory does not contain standard reports.
  } else if (report_format == "mpa" && n_std_reports > 0) {
    error_message <- paste0(EXCEPTIONS_CHECK_REPORT_DIR_STD_IN_MPA, glue::double_quote(dirpath))
    logger::log_fatal(error_message)
    stop(error_message)

  # Check that the directory contains standard reports if it is supposed to.
  } else if (report_format == "std" && n_std_reports == 0) {
    error_message <- paste0(EXCEPTIONS_CHECK_REPORT_DIR_NO_STD, glue::double_quote(dirpath))
    logger::log_fatal(error_message)
    stop(error_message)

  # Check that the standard reports directory does not contain MPA-style reports.
  } else if (report_format == "std" && n_mpa_reports > 0) {
    error_message <- paste0(EXCEPTIONS_CHECK_REPORT_DIR_MPA_IN_STD, glue::double_quote(dirpath))
    logger::log_fatal(error_message)
    stop(error_message)

  # Ensure that the report format specified is valid.
  } else if (!(report_format %in% c("std", "mpa"))) {
    error_message <- paste0(EXCEPTIONS_CHECK_REPORT_DIR_INVALID, glue::double_quote(report_format))
    logger::log_fatal(error_message)
    stop(error_message)

  # Check that the directory does not contain both report formats.
  } else if (n_mpa_reports > 0 && n_std_reports > 0) {
    error_message <- paste0(EXCEPTIONS_CHECK_REPORT_DIR_BOTH_FORMATS, glue::double_quote(dirpath))
    logger::log_fatal(error_message)
    stop(error_message)
  }

  return(dirpath)
}



#' Carry out checks on a given file
#'
#' This function takes the path to a given file and performs a few checks:
#' 1) First, it checks if the file exists.
#' 2) Second, it checks if the file is empty.
#'
#' @param file_path Path to a given file.
#'
check_file <- function(file_path) {
  if (!(file.exists(file_path))) {
    error_message <- paste0(EXCEPTIONS_CHECKFILE_DOES_NOT_EXIST, glue::double_quote(file_path))
    logger::log_fatal(error_message)
    stop(error_message)
  }
  else if (file.size(file_path) == 0L) {
    error_message <- paste0(EXCEPTIONS_CHECKFILE_EMPTY, glue::double_quote(file_path))
    logger::log_fatal(error_message)
    stop(error_message)
  }
}

#' Check if a set of column names is present in a given dataframe
#'
#' This function takes a dataframe and a vector containing column names that
#' are expected to be found in the given dataframe; then it checks if each
#' column is indeed present, throwing an error otherwise.
#'
#' @param df Dataframe where given column names will be looked for.
#' @param columns A vector with column names that are expected to be present
#'  in the given dataframe.
#'
check_columns <- function(df, columns) {
  for (column in columns) {
    if (!(column %in% colnames(df))) {
      error_message <- paste0(EXCEPTIONS_CHECK_COLUMNS, glue::double_quote(column))
      logger::log_fatal(error_message)
      stop(error_message)
    }
  }
}

#' Carry out checks on a prefix
#'
#' This function takes a prefix and performs a few checks:
#' 1) First, it checks if the prefix corresponds to the appropriate type ("character").
#' 2) Second, it ensures that there is an underscore symbol ("_") at the end of the prefix.
#'
#' @param prefix A given prefix that will later be added to output file names.
#' @return Checked prefix.
#'
check_prefix <- function(prefix) {
  if (!(is.character(prefix))) {
    error_message <- paste0(EXCEPTIONS_CHECK_PREFIX_NOT_CHARACTER, glue::double_quote(typeof(prefix)))
    logger::log_fatal(error_message)
    stop(error_message)
  }
  else if (substr(prefix, nchar(prefix), nchar(prefix)) != "_") prefix <- paste0(prefix, "_")

  logger::log_trace("")

  return(prefix)
}

#' Check if the domain specified is valid
#'
#' This function takes a domain and checks if it is valid, i.e. one of "Viruses",
#' "Bacteria", "Archaea", or "Eukaryota".
#'
#' @param domain A domain ("Viruses", "Bacteria", "Archaea", or "Eukaryota").
#'
check_domain <- function(domain) {
  if (!(domain %in% c("Viruses", "Bacteria", "Archaea", "Eukaryota"))) {
    error_message <- paste0(EXCEPTIONS_CHECK_DOMAIN_NOT_VALID, glue::single_quote(domain))
    logger::log_fatal(error_message)
    stop(error_message)
  }
}

#' Check if a given species is present in the report
#'
#' This function takes a species name (e.g. "Homo sapiens") and checks if it is present
#' in the report.
#'
#' @param report A report loaded with load_STDreports() or generated by mergeReports().
#' @param species The name of the species being analysed. For example, if the input data
#'  for Kraken2 comes from human samples, then the species name should be "Homo sapiens".
#'
check_species_in_report <- function(report, species) {
  # Here the species name will be specified by the user, so it might be that it will not
  # be present in the report exactly as provided by the user.

  # If the species name specified by the user is not present in the report as is...
  if (species %in% report[[COLNAME_STD_TAXON]]) {
    logger::log_debug(paste0(
      "The species ", species, " was found in the report at rows ",
      paste(which(report[[COLNAME_STD_TAXON]] == species), collapse = ", "), "."
    ))

    return(species)

  } else if (!(species %in% report[[COLNAME_STD_TAXON]])) {

    logger::log_debug(
      paste0("Species ", glue::single_quote(species), " was not found in the report."
    ))

    # Check if the species name specified by the user has an underscore...
    if (grepl("_", species)) {
      species_alt <- look_for_alternative_name(report, species, "_")

    # Alternatively, check if the species name specified by the user has a white space...
    } else if (grepl("[[:space:]]", species)) {
      species_alt <- look_for_alternative_name(report, species, " ")

    # If none of the above works, raise an error as the species name specified was not
    # found at all in the report (bearing in mind that the delimiters supported here
    # are "_" [underscore] and " " [white space]).
    } else {
      error_message <- paste0(
        "The species name specified, ", glue::single_quote(species), ", does not contain any ",
        "of the supported delimiters (white space [' '] or underscore ['_']). Please review your input!"
      )
      logger::log_fatal(error_message)
      stop(error_message)
    }

    return(species_alt)
  }
}

look_for_alternative_name <- function(report, species, delimiter) {
  logger::log_debug(
    paste0("Species ", glue::single_quote(species), " is delimited with ", delimiter, "."
  ))

  # Replace the delimiter.
  if (delimiter == "_") {
    species_alt <- gsub("_", " ", species)
  } else if (delimiter == " ") {
    species_alt <- gsub(" ", "_", species)
  } else {
    stop("Delimiter not supported.")
  }

  logger::log_debug(paste0(
    "The altenative version of ", glue::single_quote(species), " is ",
    glue::single_quote(species_alt), "; checking if it is present in the report..."
  ))
  # Look for alternative species name in report...
  if (species_alt %in% report[[COLNAME_STD_TAXON]]) {
    logger::log_warn(paste0(
      glue::single_quote(species), " was not found in the report, but ",
      glue::single_quote(species_alt), " was. Carrying on with ", glue::single_quote(species_alt),
      "..."
    ))

    # Return the alternative species name if it is present.
    return(species_alt)

  # Raise an error if the alternative name is not found.
  } else {
    error_message <- paste0(
      "Neither ", glue::single_quote(species), " or ", glue::single_quote(species_alt),
      " were found in the report. Please review your input!"
    )
    logger::log_fatal(error_message)
    stop(error_message)
  }
}
