#' PREPARE DATA
#'
#' @param std_reports_path Path to directory containing standard Kraken2 reports.
#' @param mpa_reports_path Path to directory containing MPA-style reports.
#' @param reference_path Path to file containing information on a Kraken2 reference database.
#' @param metadata_path Path to metadata file.
#' @param verbose Verbosity level.
#' @return Processed data ready for downstream steps.
#' @export
prepare_data <- function(std_reports_path, mpa_reports_path, reference_path, metadata_path, verbose) {

    # Check the integrity of the directories provided & read reports.
    reports_paths <- list(
        "std" = std_reports_path, 
        "mpa" = mpa_reports_path
    )
    for (report_format in names(reports_paths)) {
        reports_paths[[path]] <- check_report_directory(
            dirpath = reports_paths[[path]], 
            report_format = report_format
        )
    }

    # Check the integrity of the reference and metadata (if provided) files.
    if (is.na(metadata_path)) warning("No metadata file has been provided.")
    for (file_path in c(reference_path, metadata_path)) { 
        check_file(file_path = file_path) 
    }

    # Read standard reports.
    std_reports <- load_STDreports(reports_paths[[1]], verbose = FALSE)
    
    # Read MPA-style reports.
    mpa_reports <- load_MPAreports(reports_paths[[2]], verbose = FALSE)
    
    # Read reference database data.
    ref_db <- loadReference(reference_path)

    # Read metadata (if any).
    mdata <- NA
    if (!(is.na(metadata_path))) mdata <- loadMetadata(metadata_path)
    
    return(list(std_reports, mpa_reports, ref_db, mdata))

}

#' PROCESS DATA
#'
#' @param std_reports_path Path to directory containing standard Kraken2 reports.
#' @param mpa_reports_path Path to directory containing MPA-style reports.
#' @param reference_path Path to file containing information on a Kraken2 reference database.
#' @param metadata_path Path to metadata file.
#' @param verbose Verbosity level.
#' @return Processed reports ready for downstream analysis.
#' @export
process_kraken2 <- function(
    std_reports_path, 
    mpa_reports_path, 
    reference_path, 
    metadata_path, 
    metadata_columns, 
    verbose
) {

    prepared_data <- prepare_data(
        std_reports_path = std_reports_path,
        mpa_reports_path = mpa_reports_path,
        reference_path = reference_path,
        metadata_path = metadata_path,
        verbose = verbose
    )

    std_reports <- prepared_data[[1]]
    mpa_reports <- prepared_data[[2]]
    ref_db <- prepared_data[[3]]
    mdata <- prepared_data[[4]]
    
    mdata_columns <- c("Diagnosis_short", "Site_group")

    mpa_reports <- addMetadata(mpa_reports, mdata, mdata_columns)
    std_reports <- addMetadata(std_reports, mdata, mdata_columns)


}