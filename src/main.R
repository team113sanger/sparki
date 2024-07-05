process_kraken2 <- function(std_reports_path, mpa_reports_path, metadata_path, ref_db_path) {

    ############################
    # Check report directories #
    ############################

    report_paths <- c(std_reports_path, mpa_reports_path)
    names(report_paths) <- c("STD", "MPA")

    for (path in report_paths) {
        check_report_directory(directory_path = path, report_format = names(report_paths)[report_paths == path])
    }

    #######################################
    # Check metadata & reference DB files #
    #######################################

    check_reference(ref_db_path)


    ########################
    # Load Kraken2 reports #
    ########################

    std_reports <- load_STDreports(std_reports_dir = std_reports_path, verbose = FALSE)
    mpa_reports <- load_MPAreports(mpa_reports_dir = mpa_reports_path, verbose = FALSE)

    ###########################
    # Process Kraken2 reports #
    ###########################




}