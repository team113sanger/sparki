#' Check input arguments before running the SPARKI analysis
#'
#' @param std_reports_path Path to a directory containing standard reports.
#' @param mpa_reports_path Path to a directory containing MPA-style reports.
#' @param organism Species-level name of the organism being analysed.
#' @param reference_path Path to 'inspect.txt' file in a Kraken2 reference database.
#' @param metadata_path Path to metadata table file.
#' @param metadata_sample_col Name of column containing sample IDs in the metadata table.
#' @param metadata_columns Name of columns of interest in the metadata table.
#' @param outdir_path Path to output directory.
#' @param prefix Prefix to be added to output file names.
#' @param domain Domain(s) of interest.
#' @param samples_to_remove_path Path to file containing a list of samples that should
#'  be removed from the SPARKI analysis.
#'
#' @return A list containing the input arguments after checks.
#'
check_inputs <- function(
    std_reports_path,
    mpa_reports_path,
    organism,
    reference_path,
    metadata_path,
    metadata_sample_col,
    metadata_columns,
    outdir_path,
    prefix,
    domain,
    samples_to_remove_path) {

  ######################
  # CHECKS FOR REPORTS #
  ######################

  # Carry out checks on the directories with standard and MPA-style reports.
  if (is.na(std_reports_path)) stop(MAIN_NO_STD_PROVIDED)
  std_reports_path <- check_report_directory(dirpath = std_reports_path, report_format = "std")
  if (is.na(mpa_reports_path)) stop(MAIN_NO_MPA_PROVIDED)
  mpa_reports_path <- check_report_directory(dirpath = mpa_reports_path, report_format = "mpa")

  #######################
  # CHECKS FOR ORGANISM #
  #######################
  check_organism(std_reports_path, mpa_reports_path, organism)

  #######################
  # CHECKS FOR METADATA #
  #######################

  # Check if a metadata file has been provided.
  if (is.na(metadata_path)) {
    logger::log_info(MAIN_METADATA_NOT_PROVIDED)
  }
  else {
    # Check if column names have been specified.
    if (is.na(metadata_columns)) {
      logger::log_fatal(MAIN_NO_COLUMNS_SPECIFIED)
      stop(MAIN_NO_COLUMNS_SPECIFIED)

    # Check if a sample column has been specified.
    } else if (is.na(metadata_sample_col)) {
      logger::log_fatal(MAIN_NO_SAMPLE_COL_SPECIFIED)
      stop(MAIN_NO_SAMPLE_COL_SPECIFIED)

    # Carry out further checks and processing on the metadata file path.
    } else {
      # Check if the metadata file exists and is not empty.
      check_file(file_path = metadata_path)

      # Parse the comma-delimited column list.
      metadata_columns <- parse_delimited_list(del_list = metadata_columns, delimiter = ",")

      # Load the metadata table and check if the user-defined columns are present.
      metadata <- loadMetadata(metadata_path)
      check_columns(df = metadata, columns = c(metadata_sample_col, metadata_columns))
    }
  }

  ########################
  # CHECKS FOR REFERENCE #
  ########################

  # Check the integrity of the reference file.
  if (is.na(reference_path)) {
    logger::log_fatal(MAIN_NO_REF_PROVIDED)
    stop(MAIN_NO_REF_PROVIDED)
  }
  check_file(reference_path)

  #####################
  # CHECKS FOR PREFIX #
  #####################

  # Check prefix.
  if (prefix == "") {
    logger::log_info(MAIN_PREFIX_NOT_PROVIDED)
  }
  else prefix <- check_prefix(prefix)

  ###############################
  # CHECKS FOR OUTPUT DIRECTORY #
  ###############################

  # Check the output directory.
  if (is.na(outdir_path)) {
    logger::log_fatal(MAIN_NO_OUTDIR_PROVIDED)
    stop(MAIN_NO_OUTDIR_PROVIDED)
  }
  outdir_path <- check_directory(outdir_path, expectation = "empty")

  #####################
  # CHECKS FOR DOMAIN #
  #####################

  # Check the user-specified domain(s).
  if (is.na(domain)) {
    logger::log_fatal(MAIN_NO_DOMAIN_PROVIDED)
    stop(MAIN_NO_DOMAIN_PROVIDED)
  }
  domain <- parse_delimited_list(domain, ",")
  for (d in domain) check_domain(d)

  ####################################
  # CHECKS FOR SAMPLES TO BE REMOVED #
  ####################################

  # Check the integrity of the samples-to-remove file specified.
  if (!is.na(samples_to_remove_path)) {
    check_file(samples_to_remove_path)
  }
  else {
    logger::log_info(MAIN_NO_SAMPLES_TO_REMOVE_PROVIDED)
  }

  return(
    list(
      std_reports_path,
      mpa_reports_path,
      reference_path,
      metadata_path,
      metadata_sample_col,
      metadata_columns,
      outdir_path,
      prefix,
      domain,
      samples_to_remove_path)
  )
}

#' Load and prepare data for SPARKI analysis
#'
#' @param std_reports_path Path to a directory containing standard reports.
#' @param mpa_reports_path Path to a directory containing MPA-style reports.
#' @param organism Species-level name of the organism being analysed.
#' @param reference_path Path to 'inspect.txt' file in a Kraken2 reference database.
#' @param metadata_path Path to metadata table file.
#' @param metadata_sample_col Name of column containing sample IDs in the metadata table.
#' @param metadata_columns Name of columns of interest in the metadata table.
#' @param outdir_path Path to output directory.
#' @param prefix Prefix to be added to output file names.
#' @param domain Domain(s) of interest.
#' @param samples_to_remove_path Path to file containing a list of samples that should
#'  be removed from the SPARKI analysis.
#'
#' @return A list containing the input arguments ready for the SPARKI analysis.
#'
load_data <- function(
    std_reports_path,
    mpa_reports_path,
    organism,
    reference_path,
    metadata_path,
    metadata_sample_col,
    metadata_columns,
    outdir_path,
    prefix,
    domain,
    samples_to_remove_path) {

  # Check the user-defined inputs to ensure everything is ok
  # before starting loading and processing the data.
  checked_inputs <- check_inputs(
    std_reports_path = std_reports_path,
    mpa_reports_path = mpa_reports_path,
    organism = organism,
    reference_path = reference_path,
    metadata_path = metadata_path,
    metadata_sample_col = metadata_sample_col,
    metadata_columns = metadata_columns,
    outdir_path = outdir_path,
    prefix = prefix,
    domain = domain,
    samples_to_remove_path = samples_to_remove_path
  )
  std_reports_path <- checked_inputs[[1]]
  mpa_reports_path <- checked_inputs[[2]]
  reference_path <- checked_inputs[[3]]
  metadata_path <- checked_inputs[[4]]
  metadata_sample_col <- checked_inputs[[5]]
  metadata_columns <- checked_inputs[[6]]
  outdir_path <- checked_inputs[[7]]
  prefix <- checked_inputs[[8]]
  domain <- checked_inputs[[9]]
  samples_to_remove_path <- checked_inputs[[10]]

  #  Load the samples-to-remove file.
  samples_to_remove <- NA
  if (!is.na(samples_to_remove_path)) {
    samples_to_remove <- loadSamplesToRemove(filepath = samples_to_remove_path)
  }

  # Load standard and MPA-style reports, and merge them afterwards.
  std_reports <- load_STDreports(std_reports_path, samples_to_remove = samples_to_remove)
  mpa_reports <- load_MPAreports(mpa_reports_path, samples_to_remove = samples_to_remove)
  merged_reports <- mergeReports(std_reports, mpa_reports)

  # Load reference database data.
  ref_db <- loadReference(reference_path)

  # Load the metadata table.
  metadata <- NA
  if (!(is.na(metadata_path))) {
    metadata <- loadMetadata(metadata_path)

    # Ensure that the user-specified columns are present in the metadata table provided.
    check_columns(df = metadata, columns = c(metadata_sample_col, metadata_columns))

    # Add metadata to merged reports.
    merged_reports <- addMetadata(merged_reports, metadata, metadata_sample_col, metadata_columns)

  # Throw a warning message if no metadata table has been provided.
  } else logger::log_info(MAIN_NO_METADATA_ADDED)

  return(
    list(merged_reports, ref_db, metadata, metadata_sample_col, metadata_columns, outdir_path, prefix, domain)
  )
}

#' Run an end-to-end SPARKI analysis
#'
#' @param std_reports_path Path to a directory containing standard reports.
#' @param mpa_reports_path Path to a directory containing MPA-style reports.
#' @param organism Species-level name of the organism being analysed.
#' @param reference_path Path to 'inspect.txt' file in a Kraken2 reference database.
#' @param metadata_path Path to metadata table file.
#' @param metadata_sample_col Name of column containing sample IDs in the metadata table.
#' @param metadata_columns Name of columns of interest in the metadata table.
#' @param outdir_path Path to output directory.
#' @param prefix Prefix to be added to output file names.
#' @param include_eukaryotes Whether eukaryotes should be included in all plots (they are
#'  included in some of the plots regardless of this flag).
#' @param include_sample_names Whether sample names should be included in all plots (they
#'  are included in some of the plots regardless of this flag).
#' @param domain Domain(s) of interest.
#' @param samples_to_remove_path Path to file containing a list of samples that should
#'  be removed from the SPARKI analysis.
#'
#' @return Processed reports ready for downstream analysis.
#' @export
#'
run_analysis <- function(
    std_reports_path,
    mpa_reports_path,
    organism,
    reference_path,
    metadata_path,
    metadata_sample_col,
    metadata_columns,
    outdir_path,
    prefix,
    include_eukaryotes,
    include_sample_names,
    domain,
    samples_to_remove_path) {

  # Load data.
  loaded_data <- load_data(
    std_reports_path = std_reports_path,
    mpa_reports_path = mpa_reports_path,
    organism = organism,
    reference_path = reference_path,
    metadata_path = metadata_path,
    metadata_sample_col = metadata_sample_col,
    metadata_columns = metadata_columns,
    outdir_path = outdir_path,
    prefix = prefix,
    domain = domain,
    samples_to_remove_path = samples_to_remove_path
  )
  merged_reports <- loaded_data[[1]]
  ref_db <- loaded_data[[2]]
  mdata <- loaded_data[[3]]
  sample_col <- loaded_data[[4]]
  columns <- loaded_data[[5]]
  outdir <- loaded_data[[6]]
  prefix <- loaded_data[[7]]
  domain <- loaded_data[[8]]

  # Add sample sizes and minimiser data to the merged reports table.
  merged_reports <- addSampleSize(merged_reports)
  merged_reports <- addMinimiserData(merged_reports, ref_db)
  merged_reports <- add_nTaxaInRank(merged_reports)

  # Make a violin plot with the classification summary.
  plotClassificationSummary_violin(
    merged_reports,
    return_plot = FALSE,
    outdir = outdir,
    prefix = prefix
  )

  # Make a barplot with the classification summary.
  plotClassificationSummary_barplot(
    merged_reports,
    include_sample_names = include_sample_names,
    orientation = "horizontal",
    return_plot = FALSE,
    outdir = outdir,
    prefix = prefix
  )

  plotClassificationProportion(
    merged_reports,
    return_plot = FALSE,
    outdir = outdir,
    prefix = prefix
  )

  plotDistribution_histogram(
    merged_reports,
    return_plot = FALSE,
    outdir = outdir,
    prefix = prefix
  )

  plotDistribution_violin(
    merged_reports,
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

  readr::write_csv(
    merged_reports,
    append = FALSE,
    col_names = TRUE,
    quote = "needed",
    paste0(outdir, prefix, "pre_filtering_and_statistics.csv")
  )

  merged_reports <- subsetReports(merged_reports, species = organism)
  merged_reports <- assessMinimiserRatio(merged_reports)
  merged_reports <- assessStatistics(merged_reports, ref_db)

  plotSignificanceSummary(
    merged_reports,
    return_plot = FALSE,
    outdir = outdir,
    prefix = prefix
  )

  for (d in domain) {
    plotMinimisers_dotplot(
      merged_reports,
      domain = d,
      return_plot = FALSE,
      fig_width = 25,
      fig_height = 15,
      outdir = outdir,
      prefix = prefix
    )
  }

  readr::write_csv(
    merged_reports,
    append = FALSE,
    col_names = TRUE,
    quote = "needed",
    paste0(outdir, prefix, "final_table_with_pvalues.csv")
  )
}
