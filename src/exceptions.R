check_report_directory <- function(dirpath, report_format) {

    # Check that the directory exists and is not empty.
    dirpath <- check_directory(dirpath)

    # Check that the directory contains MPA reports if it 
    # is supposed to.
    if (
        report_format == "mpa" && 
        length(list.files(dirpath, pattern = "mpa$")) == 0
    ) {
        stop(paste0(
            "The directory ", dirpath, " does not contain ",
            "any MPA-style reports. Please review your input!"
        ))

    # Check that the directory contains standard reports if it 
    # is supposed to.
    } else if (
        report_format == "std" && 
        length(list.files(dirpath, pattern = "kraken$")) == 0
    ) {
        stop(paste0(
            "The directory ", dirpath, " does not contain ",
            "any standard format reports. Please review your input!"
        ))

    # Double-check that the report format specified is valid.
    } else if (!(report_format %in% c("std", "mpa"))) {
        stop(paste0("The format ", report_format, " is not valid."))
    
    # Check that the directory is not empty.
    } else if (length(list.files(dirpath)) == 0) {
        stop(paste0(
            "The directory ", dirpath, " is empty. ",
            "Please review your input!"
        ))

    return(dirpath)
}

check_file <- function(file_path) {

    if (!(file.exists(file_path))) {
        stop(paste0(
            "The file ", file_path, " does not exist. ",
            "Please review your input!"
        ))
    } else if (file.size(file_path) == 0L) {
        stop(paste0(
            "The file ", file_path, " is empty. ",
            "Please review your input!"
        ))
    }
}

check_directory <- function(dirpath) {

    # Check that the directory exists.
    if (!(dir.exists(dirpath)))  {
        stop(paste0(
            "The directory ", dirpath, " does not exist. ",
            "Please review your input!"
        ))
    
    # Ensure that the directory path has a slash at the end.
    } else if (substr(dirpath, nchar(dirpath), nchar(dirpath)) != "/") {
        dirpath <- paste0(dirpath, "/")
    }

    return(dirpath)
}

check_columns <- function(df, columns) {
    
    for (column in columns) {
        if (!(column %in% colnames(df))) {
            stop(paste0(
                "The column ", column, " does not exist ",
                "in the metadata table provided. Please ",
                "review your input!"
            ))
        }
    }
}

check_prefix <- function(prefix) {

    if (substr(prefix, nchar(prefix), nchar(prefix)) != "_") {
        prefix <- paste0(prefix, "_")
    }

    return(prefix)
}