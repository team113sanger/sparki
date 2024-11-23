describe("Functions in main.R", {
  test_that("check_inputs() returns valid outputs", {
    # Run SPARKI's check_inputs()
    out <- SPARKI:::check_inputs(
      std_reports_path = get_test_std_report_dir("usual"),
      mpa_reports_path = get_test_mpa_report_dir("usual"),
      reference_path = get_test_reference(),
      metadata_path = get_test_metadata(),
      metadata_sample_col = "sample",
      metadata_columns = "type,date,status",
      outdir_path = get_local_tmp_dir(),
      prefix = "test",
      verbose = FALSE,
      domain = "Viruses",
      samples_to_remove_path = get_test_samples_to_remove()
    )

    ############################################################
    ## Carry out tests on the standard reports directory path ##
    ############################################################
    std_reports_path <- out[[1]]

    # Carry out tests.
    expect_true(dir.exists(std_reports_path)) # Ensure directory exists.
    expect_gt(length(list.files(std_reports_path)), 0) # Ensure directory is not empty.
    expect_gt(length(list.files(std_reports_path, pattern = "kraken$")), 0) # Ensure there are standard reports.
    expect_equal(length(list.files(std_reports_path, pattern = "mpa$")), 0) # Ensure there are no MPA-style reports.

    #############################################################
    ## Carry out tests on the MPA-style reports directory path ##
    #############################################################
    mpa_reports_path <- out[[2]]

    # Carry out tests.
    expect_true(dir.exists(mpa_reports_path)) # Ensure directory exists.
    expect_gt(length(list.files(mpa_reports_path)), 0) # Ensure directory is not empty.
    expect_gt(length(list.files(mpa_reports_path, pattern = "mpa$")), 0) # Ensure there are MPA-stype reports.
    expect_equal(length(list.files(mpa_reports_path, pattern = "kraken$")), 0) # Ensure there are no standard reports.

    #########################################################
    ## Carry out tests on the reference database file path ##
    #########################################################
    reference_path <- out[[3]]

    # Carry out tests
    expect_true(file.exists(reference_path)) # Ensure file exists.
    expect_false((file.size(reference_path) == 0L)) # Ensure file is not empty.

    # Carry out tests on the metadata-related outputs.
    metadata_path <- out[[4]]
    metadata_sample_col <- out[[5]]
    metadata_columns <- out[[6]]

    expect_true(file.exists(metadata_path)) # Ensure file exists.
    expect_false((file.size(metadata_path) == 0L)) # Ensure file is not empty.

    metadata <- loadMetadata(metadata_path, verbose = FALSE) # Load metadata to double-check the specified columns are present.
    for (column in c(metadata_sample_col, metadata_columns)) {
      expect_true((column %in% colnames(metadata)))
    }

    # Carry out tests on the output directory path.
    outdir_path <- out[[7]]

    expect_true(dir.exists(outdir_path)) # Ensure directory exists.
    expect_equal(length(list.files(outdir_path)), 0) # Ensure directory is empty.

    # Carry out tests on the specified prefix.
    prefix <- out[[8]]

    expect_type(prefix, "character")
    expect_true((substr(prefix, nchar(prefix), nchar(prefix)) == "_"))

    # Carry out tests on the specified domain.
    domain <- out[[9]]

    for (d in domain) {
      expect_type(d, "character")
      expect_true((d %in% c("Viruses", "Bacteria", "Archaea", "Eukaryota")))
    }

    # Carry out tests on the samples-to-remove file path.
    samples_to_remove_path <- out[[10]]

    expect_true(file.exists(samples_to_remove_path)) # Ensure file exists.
    expect_false((file.size(samples_to_remove_path) == 0L)) # Ensure file is not empty.

  })

  test_that("load_data() returns valid outputs", {
    # Run SPARKI's load_data()
    out <- SPARKI:::load_data(
      std_reports_path = get_test_std_report_dir("usual"),
      mpa_reports_path = get_test_mpa_report_dir("usual"),
      reference_path = get_test_reference(),
      metadata_path = get_test_metadata(),
      metadata_sample_col = "sample",
      metadata_columns = "type,date,status",
      outdir_path = get_local_tmp_dir(),
      prefix = "test",
      verbose = FALSE,
      domain = "Viruses",
      samples_to_remove = get_test_samples_to_remove()
    )

    ###########################################
    ## Carry out tests on the merged reports ##
    ###########################################
    merged_reports <- out[[1]]
    columns <- out[[5]]

    # Get the expected number of columns in the merged reports dataframe.
    # This will be:
    # 9 (# columns from the standard report) +
    # 8 (# columns from the MPA-style report, after excluding the ones that overlap with the standard report) +
    # The number of metadata columns.
    expected_number_of_columns <- (9 + 8 + length(columns))

    # Define the expected column names.
    expected_colnames <- c(
      SPARKI:::COLNAME_STD_SAMPLE,
      columns,
      SPARKI:::COLNAME_STD_PCT_FRAG_CLADE,
      SPARKI:::COLNAME_STD_N_FRAG_CLADE,
      SPARKI:::COLNAME_STD_N_FRAG_TAXON,
      SPARKI:::COLNAME_STD_MINIMISERS,
      SPARKI:::COLNAME_STD_UNIQ_MINIMISERS,
      SPARKI:::COLNAME_STD_RANK,
      SPARKI:::COLNAME_STD_NCBI_ID,
      SPARKI:::COLNAME_STD_TAXON,
      SPARKI:::COLNAME_MPA_DOMAIN,
      SPARKI:::COLNAME_MPA_KINGDOM,
      SPARKI:::COLNAME_MPA_PHYLUM,
      SPARKI:::COLNAME_MPA_CLASS,
      SPARKI:::COLNAME_MPA_ORDER,
      SPARKI:::COLNAME_MPA_FAMILY,
      SPARKI:::COLNAME_MPA_GENUS,
      SPARKI:::COLNAME_MPA_SPECIES
    )

    # Get the actual number of rows and columns in the merged reports dataframe.
    actual_number_of_rows <- nrow(merged_reports)
    actual_number_of_columns <- ncol(merged_reports)

    # Get the actual column names.
    actual_colnames <- colnames(merged_reports)

    # Carry out tests.
    expect_s3_class(merged_reports, "data.frame") # Ensure the class corresponds to "data.frame".
    expect_false(actual_number_of_rows == 0) # Ensure the number of rows is not zero.
    expect_false(actual_number_of_columns == 0) # Ensure the number of columns is not zero.
    expect_false(is.null(dim(merged_reports))) # Ensure the dataframe is valid.
    expect_equal(expected_number_of_columns, actual_number_of_columns) # Ensure the number of columns is right.
    for (column in columns) expect_true((column %in% colnames(merged_reports))) # Ensure each metadata column is present.
    expect_equal(expected_colnames, actual_colnames) # Ensure all expected column names are present.

  })

})
