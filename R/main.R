check_inputs <- function(
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
    samples_to_remove_path) {

  ######################
  # CHECKS FOR REPORTS #
  ######################

  # Carry out checks on the directories with standard and MPA-style reports.
  std_reports_path <- check_report_directory(dirpath = std_reports_path, report_format = "std")
  mpa_reports_path <- check_report_directory(dirpath = mpa_reports_path, report_format = "mpa")

  #######################
  # CHECKS FOR METADATA #
  #######################

  # Check if a metadata file has been provided.
  if (is.na(metadata_path)) message(MAIN_WARNING_METADATA_NOT_PROVIDED)
  else {
    # Check if column names have been specified.
    if (is.na(metadata_columns)) {
      stop(MAIN_ERROR_NO_COLUMNS_SPECIFIED)

    # Check if a sample column has been specified.
    } else if (is.na(metadata_sample_col)) {
      stop(MAIN_ERROR_NO_SAMPLE_COL_SPECIFIED)

    # Carry out further checks and processing on the metadata file path.
    } else {
      # Check if the metadata file exists and is not empty.
      check_file(file_path = metadata_path)

      # Parse the comma-delimited column list.
      metadata_columns <- parse_delimited_list(del_list = metadata_columns, delimiter = ",")

      # Load the metadata table and check if the user-defined columns are present.
      metadata <- loadMetadata(metadata_path, verbose = FALSE)
      check_columns(df = metadata, columns = c(metadata_sample_col, metadata_columns))
    }
  }

  ########################
  # CHECKS FOR REFERENCE #
  ########################

  # Check the integrity of the reference file.
  check_file(reference_path)

  #####################
  # CHECKS FOR PREFIX #
  #####################

  # Check prefix.
  if (prefix == "") message(MAIN_WARNING_PREFIX_NOT_PROVIDED)
  else prefix <- check_prefix(prefix)

  ###############################
  # CHECKS FOR OUTPUT DIRECTORY #
  ###############################

  # Check the output directory.
  outdir_path <- check_directory(outdir_path, expectation = "empty")

  #####################
  # CHECKS FOR DOMAIN #
  #####################

  # Check the user-specified domain(s).
  domain <- parse_delimited_list(domain, ",")
  for (d in domain) check_domain(d)

  ####################################
  # CHECKS FOR SAMPLES TO BE REMOVED #
  ####################################

  # Check the integrity of the samples-to-remove file specified.
  if (!is.na(samples_to_remove_path)) check_file(samples_to_remove_path)
  else message(MAIN_WARNING_NO_SAMPLES_TO_REMOVE_PROVIDED)

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

#' PREPARE DATA
#'
#' @param std_reports_path Path to directory containing standard Kraken2 reports.
#' @param mpa_reports_path Path to directory containing MPA-style reports.
#' @param reference_path Path to file containing information on a Kraken2 reference database.
#' @param metadata_path Path to metadata file.
#' @param verbose Verbosity level.
#' @return Processed data ready for downstream steps.
#'
load_data <- function(
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
    samples_to_remove_path) {

  # Check the user-defined inputs to ensure everything is ok
  # before starting loading and processing the data.
  checked_inputs <- check_inputs(
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
    samples_to_remove <- loadSamplesToRemove(filepath = samples_to_remove_path, verbose = verbose)
  }

  # Load standard and MPA-style reports, and merge them afterwards.
  std_reports <- load_STDreports(std_reports_path, samples_to_remove = samples_to_remove, verbose = verbose)
  mpa_reports <- load_MPAreports(mpa_reports_path, samples_to_remove = samples_to_remove, verbose = verbose)
  merged_reports <- mergeReports(std_reports, mpa_reports, verbose = verbose)

  # Load reference database data.
  ref_db <- loadReference(reference_path)

  # Load the metadata table.
  metadata <- NA
  if (!(is.na(metadata_path))) {
    metadata <- loadMetadata(metadata_path, verbose = TRUE)

    # Ensure that the user-specified columns are present in the metadata table provided.
    check_columns(df = metadata, columns = c(metadata_sample_col, metadata_columns))

    # Add metadata to merged reports.
    merged_reports <- addMetadata(merged_reports, metadata, metadata_sample_col, metadata_columns, verbose = verbose)

  # Throw a warning message if no metadata table has been provided.
  } else message(MAIN_WARNING_NO_METADATA_ADDED)

  return(
    list(merged_reports, ref_db, metadata, metadata_sample_col, metadata_columns, outdir_path, prefix, domain)
  )
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
run_sparki <- function(
    std_reports_path,
    mpa_reports_path,
    organism,
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
    samples_to_remove) {

  # Load data.
  loaded_data <- load_data(
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
    samples_to_remove = samples_to_remove
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
  merged_reports <- addSampleSize(merged_reports, verbose = verbose)
  merged_reports <- addMinimiserData(merged_reports, ref_db, verbose = verbose)

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

  write.csv(merged_reports, paste0(outdir, prefix, "pre_filtering_and_statistics.csv"))

  merged_reports <- subsetReports(merged_reports, species = organism, verbose = verbose)
  merged_reports <- assessMinimiserRatio(merged_reports, verbose = verbose)
  merged_reports <- assessStatistics(merged_reports, ref_db, verbose = verbose)

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

  write.csv(merged_reports, paste0(outdir, prefix, "final_table_with_pvalues.csv"))
}
