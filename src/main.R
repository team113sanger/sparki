#' PREPARE DATA
#'
#' @param std_reports_path Path to directory containing standard Kraken2 reports.
#' @param mpa_reports_path Path to directory containing MPA-style reports.
#' @param reference_path Path to file containing information on a Kraken2 reference database.
#' @param metadata_path Path to metadata file.
#' @param verbose Verbosity level.
#' @return Processed data ready for downstream steps.
#' @export
prepare_data <- function(
    std_reports_path, 
    mpa_reports_path, 
    reference_path, 
    metadata_path, 
    metadata_sample_col,
    metadata_columns,
    outdir_path,
    prefix,
    verbose,
    domain,
    remove
) {

    ######################
    # CHECKS FOR REPORTS #
    ######################

    # Check the integrity of the directories provided & read reports.
    reports_paths <- list(
        "std" = std_reports_path, 
        "mpa" = mpa_reports_path
    )

    for (report_format in names(reports_paths)) {

        reports_paths[[report_format]] <- check_report_directory(
            dirpath = reports_paths[[report_format]], 
            report_format = report_format
        )
    }

    #######################
    # CHECKS FOR METADATA #
    #######################

    # Check the integrity of the metadata file (if provided).
    if (is.na(metadata_path)) {
        warning("No metadata file has been provided.")
    } else {
        check_file(file_path = metadata_path)
        if (is.na(metadata_columns)) {
            stop(paste0(
                "A metadata table has been provided, but ",
                "no columns have been specified. Please ",
                "review your input!"
            ))
        } else if (is.na(metadata_sample_col)) {
            stop(paste0(
                "A metadata table has been provided, but ",
                "no sample column has been specified. Please ",
                "review your input!"
            ))
        } else {
            metadata_columns <- parse_delimited_list(
                del_list = metadata_columns, 
                delimiter = ","
            )
        }
    }

    ########################
    # CHECKS FOR REFERENCE #
    ########################

    # Check the integrity of the reference file.
    check_file(reference_path)

    #####################
    # CHECKS FOR PREFIX #
    #####################

    # Check prefix (if provided).
    if(is.na(prefix)) {
        warning("No prefix has been provided.")
    } else {
        prefix <- check_prefix(prefix)
    }

    ###############################
    # CHECKS FOR OUTPUT DIRECTORY #
    ###############################

    # Check the integrity of the output directory.
    outdir_path <- check_directory(outdir_path)

    #####################
    # CHECKS FOR DOMAIN #
    #####################

    # Check the integrity of the domain specified.
    check_domain(domain)

    ####################################
    # CHECKS FOR SAMPLES TO BE REMOVED #
    ####################################

    if (!is.na(remove)) {

        # Check the integrity of the samples-to-remove file specified.
        check_file(remove)

    } else { 
        warning("No list of samples to be removed has been provided.")
    }

    #############
    # LOAD DATA #
    #############

    # Read standard reports.
    std_reports <- load_STDreports(reports_paths[[1]], verbose = FALSE)
    
    # Read MPA-style reports.
    mpa_reports <- load_MPAreports(reports_paths[[2]], verbose = FALSE)

    samples_to_remove <- NA
    if (!is.na(remove)) {

        # Load samples-to-remove file.
        samples_to_remove <- loadSamplesToRemove(remove)   
    }

    # Merge standard and MPA-style reports.
    merged_reports <- mergeReports(std_reports, mpa_reports, samples_to_remove)
    
    # Read reference database data.
    ref_db <- loadReference(reference_path)

    # Read metadata (if any).
    metadata <- NA
    if (!(is.na(metadata_path))) {

        # Load metadata table.
        metadata <- loadMetadata(metadata_path)

        # Check columns are present in metadata table.
        check_columns(
            df = metadata,
            columns = c(metadata_sample_col, metadata_columns)
        )

        # Add metadata to merged reports.
        merged_reports <- addMetadata(
            merged_reports, metadata, 
            metadata_sample_col, metadata_columns
        )

    } else {
        warning("No metadata was added to the reports.")
    }

    return(list(
        merged_reports, 
        ref_db, 
        metadata, 
        metadata_sample_col,
        metadata_columns, 
        outdir_path, 
        prefix, 
        domain
    ))

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
    metadata_sample_col,
    metadata_columns, 
    outdir_path,
    prefix,
    verbose,
    include_eukaryotes,
    include_sample_names,
    domain,
    remove
) {

    prepared_data <- prepare_data(
        std_reports_path = std_reports_path,
        mpa_reports_path = mpa_reports_path,
        reference_path = reference_path,
        metadata_path = metadata_path,
        metadata_sample_col = metadata_sample_col,
        metadata_columns = metadata_columns,
        outdir_path = outdir_path,
        prefix = prefix,
        verbose = verbose,
        domain = domain,
        remove = remove
    )

    merged_reports <- prepared_data[[1]]
    ref_db <- prepared_data[[2]]
    mdata <- prepared_data[[3]]
    sample_col <- prepared_data[[4]]
    columns <- prepared_data[[5]]
    outdir <- prepared_data[[6]]
    prefix <- prepared_data[[7]]
    domain <- prepared_data[[8]]

    merged_reports <- addSampleSize(merged_reports)
    merged_reports <- addMinimiserData(merged_reports, ref_db)

    plotClassificationSummary_violin(
        merged_reports, 
        return_plot = FALSE,
        outdir = outdir,
        prefix = prefix
    )

    plotClassificationSummary_barplot(
        merged_reports, 
        include_sample_names = include_sample_names, 
        orientation = "horizontal",
        return_plot = FALSE,
        outdir = outdir,
        prefix = prefix
    )

    plotDomainReads_violin(
        merged_reports, 
        include_eukaryotes = include_eukaryotes, 
        return_plot = FALSE,
        outdir = outdir,
        prefix = prefix
    )

    plotDomainReads_barplot(
        merged_reports, 
        include_eukaryotes = include_eukaryotes, 
        include_sample_names = include_sample_names, 
        orientation = "horizontal", 
        return_plot = FALSE,
        outdir = outdir,
        prefix = prefix
    )

    write.csv(
        merged_reports,
        paste0(outdir, prefix, "pre_filtering_and_statistics.csv")
    )

    merged_reports <- subsetReports(merged_reports, include_human = FALSE)
    merged_reports <- assessMinimiserRatio(merged_reports)
    merged_reports <- assessStatistics(merged_reports, ref_db, verbose = TRUE)

    plotMinimisers_dotplot(
        merged_reports, 
        domain = domain, 
        return_plot = FALSE, 
        fig_width = 25, 
        fig_height = 15, 
        outdir = outdir,
        prefix = prefix
    )   

    write.csv(
        merged_reports,
        paste0(outdir, prefix, "final_table_with_pvalues.csv")
    )
}