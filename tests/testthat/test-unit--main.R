describe("Functions in main.R", {
  test_that("check_inputs() returns valid outputs", {
    # Run SPARKI's check_inputs()
    out <- SPARKI:::check_inputs(
      std_reports_path = get_test_std_report_dir("usual"),
      mpa_reports_path = get_test_mpa_report_dir("usual"),
      organism = "Homo sapiens",
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

    ##--------------------------------------------------------##
    ## Carry out tests on the standard reports directory path ##
    ##--------------------------------------------------------##
    std_reports_path <- out[[1]]

    # Carry out tests.
    expect_true(dir.exists(std_reports_path)) # Ensure directory exists.
    expect_gt(length(list.files(std_reports_path)), 0) # Ensure directory is not empty.
    expect_gt(length(list.files(std_reports_path, pattern = "kraken$")), 0) # Ensure there are standard reports.
    expect_equal(length(list.files(std_reports_path, pattern = "mpa$")), 0) # Ensure there are no MPA-style reports.
    expect_false(is.na(std_reports_path)) # Ensure path is not NA.
    expect_false(is.null(std_reports_path)) # Ensure path is not NULL.

    ##---------------------------------------------------------##
    ## Carry out tests on the MPA-style reports directory path ##
    ##---------------------------------------------------------##
    mpa_reports_path <- out[[2]]

    # Carry out tests.
    expect_true(dir.exists(mpa_reports_path)) # Ensure directory exists.
    expect_gt(length(list.files(mpa_reports_path)), 0) # Ensure directory is not empty.
    expect_gt(length(list.files(mpa_reports_path, pattern = "mpa$")), 0) # Ensure there are MPA-stype reports.
    expect_equal(length(list.files(mpa_reports_path, pattern = "kraken$")), 0) # Ensure there are no standard reports.
    expect_false(is.na(mpa_reports_path)) # Ensure path is not NA.
    expect_false(is.null(mpa_reports_path)) # Ensure path is not NULL.

    ##-----------------------------------------------------##
    ## Carry out tests on the reference database file path ##
    ##-----------------------------------------------------##
    reference_path <- out[[3]]

    # Carry out tests.
    expect_true(file.exists(reference_path)) # Ensure file exists.
    expect_false((file.size(reference_path) == 0L)) # Ensure file is not empty.
    expect_false(is.na(reference_path)) # Ensure path is not NA.
    expect_false(is.null(reference_path)) # Ensure path is not NULL.

    ##-------------------------------------------------##
    ## Carry out tests on the metadata-related outputs ##
    ##-------------------------------------------------##
    metadata_path <- out[[4]]
    metadata_sample_col <- out[[5]]
    metadata_columns <- out[[6]]

    # Carry out tests.
    expect_true(file.exists(metadata_path)) # Ensure file exists.
    expect_false((file.size(metadata_path) == 0L)) # Ensure file is not empty.
    expect_false(is.na(metadata_path)) # Ensure path is not NA.
    expect_false(is.null(metadata_path)) # Ensure path is not NULL.

    # Load metadata to double-check the specified columns are present.
    metadata <- loadMetadata(metadata_path, verbose = FALSE)
    for (column in c(metadata_sample_col, metadata_columns)) {
      expect_true((column %in% colnames(metadata)))
    }

    ##----------------------------------------------##
    ## Carry out tests on the output directory path ##
    ##----------------------------------------------##
    outdir_path <- out[[7]]

    # Carry out tests.
    expect_true(dir.exists(outdir_path)) # Ensure directory exists.
    expect_equal(length(list.files(outdir_path)), 0) # Ensure directory is empty.
    expect_false(is.na(outdir_path)) # Ensure path is not NA.
    expect_false(is.null(outdir_path)) # Ensure path is not NULL.

    ##-------------------------------##
    ## Carry out tests on the prefix ##
    ##-------------------------------##
    prefix <- out[[8]]

    # Carry out tests.
    expect_type(prefix, "character") # Ensure the prefix corresponds to the right type.
    expect_true((substr(prefix, nchar(prefix), nchar(prefix)) == "_")) # Ensure the last character is a slash symbol.
    expect_false(is.na(prefix)) # Ensure prefix is not NA.
    expect_false(is.null(prefix)) # Ensure prefix is not NULL.

    ##----------------------------------##
    ## Carry out tests on the domain(s) ##
    ##----------------------------------##
    domain <- out[[9]]

    # Carry out tests.
    for (d in domain) {
      expect_type(d, "character") # Ensure the domain corresponds to the right type.
      expect_true((d %in% c("Viruses", "Bacteria", "Archaea", "Eukaryota"))) # Ensure the domain specified is valid.
      expect_false(is.na(domain)) # Ensure domain is not NA.
      expect_false(is.null(domain)) # Ensure domain is not NULL.
    }

    ##----------------------------------------------------##
    ## Carry out tests on the samples-to-remove file path ##
    ##----------------------------------------------------##
    samples_to_remove_path <- out[[10]]

    # Carry out tests.
    expect_true(file.exists(samples_to_remove_path)) # Ensure file exists.
    expect_false((file.size(samples_to_remove_path) == 0L)) # Ensure file is not empty.
    expect_false(is.na(samples_to_remove_path)) # Ensure path is not NA.
    expect_false(is.null(samples_to_remove_path)) # Ensure path is not NULL.

  })

  test_that("load_data() returns valid outputs", {
    # Run SPARKI's load_data()
    out <- SPARKI:::load_data(
      std_reports_path = get_test_std_report_dir("usual"),
      mpa_reports_path = get_test_mpa_report_dir("usual"),
      organism = "Homo sapiens",
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

    ##---------------------------------------##
    ## Carry out tests on the merged reports ##
    ##---------------------------------------##
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
    expect_false(is.null(merged_reports)) # Ensure it is not NULL.
    expect_s3_class(merged_reports, "data.frame") # Ensure the class corresponds to "data.frame".
    expect_false(actual_number_of_rows == 0) # Ensure the number of rows is not zero.
    expect_false(actual_number_of_columns == 0) # Ensure the number of columns is not zero.
    expect_false(is.null(dim(merged_reports))) # Ensure the dataframe is valid.
    expect_equal(actual_number_of_columns, expected_number_of_columns) # Ensure the number of columns is right.
    for (column in columns) expect_true((column %in% colnames(merged_reports))) # Ensure each metadata column is present.
    expect_identical(actual_colnames, expected_colnames) # Ensure all expected column names are present.
    # Ensure each column corresponds to the appropriate type.
    expect_type(merged_reports[[SPARKI:::COLNAME_STD_SAMPLE]], "character")
    expect_type(merged_reports[[SPARKI:::COLNAME_STD_PCT_FRAG_CLADE]], "double")
    expect_type(merged_reports[[SPARKI:::COLNAME_STD_N_FRAG_CLADE]], "integer")
    expect_type(merged_reports[[SPARKI:::COLNAME_STD_N_FRAG_TAXON]], "integer")
    expect_type(merged_reports[[SPARKI:::COLNAME_STD_MINIMISERS]], "integer")
    expect_type(merged_reports[[SPARKI:::COLNAME_STD_UNIQ_MINIMISERS]], "integer")
    expect_type(merged_reports[[SPARKI:::COLNAME_STD_RANK]], "character")
    expect_type(merged_reports[[SPARKI:::COLNAME_STD_NCBI_ID]], "character")
    expect_type(merged_reports[[SPARKI:::COLNAME_STD_TAXON]], "character")
    expect_type(merged_reports[[SPARKI:::COLNAME_MPA_DOMAIN]], "character")
    expect_type(merged_reports[[SPARKI:::COLNAME_MPA_KINGDOM]], "character")
    expect_type(merged_reports[[SPARKI:::COLNAME_MPA_PHYLUM]], "character")
    expect_type(merged_reports[[SPARKI:::COLNAME_MPA_CLASS]], "character")
    expect_type(merged_reports[[SPARKI:::COLNAME_MPA_ORDER]], "character")
    expect_type(merged_reports[[SPARKI:::COLNAME_MPA_FAMILY]], "character")
    expect_type(merged_reports[[SPARKI:::COLNAME_MPA_GENUS]], "character")
    expect_type(merged_reports[[SPARKI:::COLNAME_MPA_SPECIES]], "character")

    ##-----------------------------------------------------##
    ## Carry out tests on the reference database dataframe ##
    ##-----------------------------------------------------##
    ref_db <- out[[2]]

    # Get the expected number of columns in the reference database dataframe.
    expected_number_of_columns <- 6

    # Define the expected column names.
    expected_colnames <- c(
      SPARKI:::COLNAME_REF_DB_PCT_FRAG_CLADE,
      SPARKI:::COLNAME_REF_DB_MINIMISERS_CLADE,
      SPARKI:::COLNAME_REF_DB_MINIMISERS_TAXON,
      SPARKI:::COLNAME_REF_DB_RANK,
      SPARKI:::COLNAME_REF_DB_NCBI_ID,
      SPARKI:::COLNAME_REF_DB_TAXON
    )

    # Get the actual number of rows and columns in the merged reports dataframe.
    actual_number_of_rows <- nrow(ref_db)
    actual_number_of_columns <- ncol(ref_db)

    # Get the actual column names.
    actual_colnames <- colnames(ref_db)

    # Carry out tests.
    expect_false(is.null(ref_db)) # Ensure it is not NULL.
    expect_s3_class(ref_db, "data.frame") # Ensure the class corresponds to "data.frame".
    expect_false(actual_number_of_rows == 0) # Ensure the number of rows is not zero.
    expect_false(actual_number_of_columns == 0) # Ensure the number of columns is not zero.
    expect_false(is.null(dim(ref_db))) # Ensure the dataframe is valid.
    expect_equal(actual_number_of_columns, expected_number_of_columns) # Ensure the number of columns is right.
    expect_identical(actual_colnames, expected_colnames) # Ensure all expected column names are present.
    # Ensure each column corresponds to the appropriate type.
    expect_type(ref_db[[SPARKI:::COLNAME_REF_DB_PCT_FRAG_CLADE]], "double")
    expect_type(ref_db[[SPARKI:::COLNAME_REF_DB_MINIMISERS_CLADE]], "double")
    expect_type(ref_db[[SPARKI:::COLNAME_REF_DB_MINIMISERS_TAXON]], "integer")
    expect_type(ref_db[[SPARKI:::COLNAME_REF_DB_RANK]], "character")
    expect_type(ref_db[[SPARKI:::COLNAME_REF_DB_NCBI_ID]], "character")
    expect_type(ref_db[[SPARKI:::COLNAME_REF_DB_TAXON]], "character")
  })

  test_that("run_sparki() generates all the output files it should", {

    # Define the output directory and ensure it is empty.
    outdir <- get_local_tmp_dir()
    expect_true((length(list.files(outdir)) == 0))

    # Run SPARKI's run_sparki()
    out <- SPARKI:::run_sparki(
      std_reports_path = get_test_std_report_dir("usual"),
      mpa_reports_path = get_test_mpa_report_dir("usual"),
      organism = "Homo sapiens",
      reference_path = get_test_reference(),
      metadata_path = get_test_metadata(),
      metadata_sample_col = "sample",
      metadata_columns = "type,date,status",
      outdir_path = outdir,
      prefix = "test",
      verbose = FALSE,
      include_eukaryotes = FALSE,
      include_sample_names = FALSE,
      domain = "Viruses",
      samples_to_remove = get_test_samples_to_remove()
    )

    ##------------------------------------------------##
    ## ENSURE THE EXPECTED OUTPUT FILES ARE GENERATED ##
    ##------------------------------------------------##

    # Check that the output directory is not empty after the
    # function run_sparki() is run.
    expect_true((length(list.files(outdir)) > 0))

    # Define the expected number of output PDF files in the
    # output directory. In this case, it will be 9 PDF files
    # because only one domain has been specified (for 2 domains
    # this number would have been 10, for example).
    expected_n_pdf_files <- 9

    # Get the actual number of PDF files in the output directory.
    actual_n_pdf_files <- length(fs::dir_ls(outdir, glob = "*.pdf$"))

    # Check that the expected amount of PDF files has been generated.
    expect_equal(actual_n_pdf_files, expected_n_pdf_files)

    # Define the expected number of output csv files in the
    # output directory.
    expected_n_csv_files <- 2

    # Get the actual number of csv files in the output directory.
    actual_n_csv_files <- length(fs::dir_ls(outdir, glob = "*.csv$"))

    # Check that the expected amount of csv files has been generated.
    expect_equal(actual_n_csv_files, expected_n_csv_files)

    # Ensure the file sizes are not zero.
    all_outfiles <- c(fs::dir_ls(outdir, glob = "*.pdf$"), fs::dir_ls(outdir, glob = "*.csv$"))
    for (file in all_outfiles) expect_gt(file.size(file), 0)

    ##--------------------------------------------------------------##
    ## ENSURE THE UNFILTERED RESULTS TABLE IS GENERATED AS EXPECTED ##
    ##--------------------------------------------------------------##

    # Load unfiltered results table and carry out tests.
    unfiltered_res_table <- readr::read_csv(
      fs::dir_ls(outdir, glob = "*pre_filtering_and_statistics.csv$"),
      show_col_types = FALSE # Supressing messages about column types when the dataframe is created.
    )

    # Define the expected/actual number of columns in the unfiltered results table
    # and carry out test.
    expected_ncol <- 23 # 20 fixed columns + 3 metadata columns as specified above.
    actual_ncol <- ncol(unfiltered_res_table)
    expect_equal(actual_ncol, expected_ncol)

    # Define the expected/actual column names in the unfiltered results table
    # and carry out test.
    expected_colnames <- c(
      SPARKI:::COLNAME_STD_SAMPLE,
      SPARKI:::COLNAME_STD_SAMPLE_SIZE,
      "type", "date", "status", # Metadata columns specified via "metadata_columns".
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
      SPARKI:::COLNAME_MPA_SPECIES,
      SPARKI:::COLNAME_STD_DB_MINIMISERS_TAXON,
      SPARKI:::COLNAME_STD_DB_MINIMISERS_CLADE
    )
    actual_colnames <- colnames(unfiltered_res_table)
    expect_identical(actual_colnames, expected_colnames)

    ##--------------------------------------------------------------##
    ## ENSURE THE FILTERED RESULTS TABLE IS GENERATED AS EXPECTED ##
    ##--------------------------------------------------------------##

    # Load filtered results table and carry out tests.
    filtered_res_table <- readr::read_csv(
      fs::dir_ls(outdir, glob = "*final_table_with_pvalues.csv$"),
      show_col_types = FALSE # Supressing messages about column types when the dataframe is created.
    )

    # Define the expected/actual number of columns in the unfiltered results table
    # and carry out test.
    expected_ncol <- 29 # 26 fixed columns + 3 metadata columns as specified above.
    actual_ncol <- ncol(filtered_res_table)
    expect_equal(actual_ncol, expected_ncol)

    # Define the expected/actual column names in the unfiltered results table
    # and carry out test.
    expected_colnames <- c(
      SPARKI:::COLNAME_STD_SAMPLE,
      SPARKI:::COLNAME_STD_SAMPLE_SIZE,
      "type", "date", "status", # Metadata columns specified via "metadata_columns".
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
      SPARKI:::COLNAME_MPA_SPECIES,
      SPARKI:::COLNAME_STD_DB_MINIMISERS_TAXON,
      SPARKI:::COLNAME_STD_DB_MINIMISERS_CLADE,
      SPARKI:::COLNAME_STD_RATIO_TAXON,
      SPARKI:::COLNAME_STD_RATIO_CLADE,
      SPARKI:::COLNAME_STD_PVALUE,
      SPARKI:::COLNAME_STD_N_TAXA_RANK,
      SPARKI:::COLNAME_STD_PADJ,
      SPARKI:::COLNAME_STD_SIGNIF
    )
    actual_colnames <- colnames(filtered_res_table)
    expect_identical(actual_colnames, expected_colnames)

  })

})
