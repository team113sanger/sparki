check_report_directory <- function(directory_path, report_format) {

    if (!(dir.exists(directory_path)))  {
        stop(paste0(
            "The directory ", directory_path, " does not exist. ",
            "Please review your input!"
        ))
    } else if (length(list.files(directory_path)) == 0) {
        stop(paste0(
            "The directory ", directory_path, " is empty. ",
            "Please review your input!"
        ))
    } else if (
        report_format == "MPA" && 
        length(list.files(directory_path, pattern = "mpa$")) == 0
    ) {
        stop(paste0(
            "The directory ", directory_path, " does not contain ",
            "any MPA-style reports. Please review your input!"
        ))
    } else if (
        report_format == "STD" && 
        length(list.files(directory_path, pattern = "kraken$")) == 0
    ) {
        stop(paste0(
            "The directory ", directory_path, " does not contain ",
            "any standard format reports. Please review your input!"
        ))
    }
}

check_reference <- function(reference_path) {

    if (!(file.exists(reference_path)))  {
        stop(paste0(
            "The file ", reference_path, " does not exist. ",
            "Please review your input!"
        ))
    } 
}