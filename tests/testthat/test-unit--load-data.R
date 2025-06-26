describe("load_MPAreports()", {
  test_that("load_MPAreports() returns a dataframe with the appropriate dimensions", {
    logger::log_threshold(logger::OFF)

    #  Define the expected number of columns.
    expected_number_of_columns <- 12

    #  Get path to test directory.
    path_to_mpa <- get_test_mpa_report_dir("valid")

    # Load the MPA-style reports and obtain the actual numbers of rows and columns.
    out <- SPARKI::load_MPAreports(path_to_mpa)
    actual_number_of_rows <- nrow(out)
    actual_number_of_columns <- ncol(out)

    # Carry out tests.
    expect_s3_class(out, "data.frame")
    expect_false(actual_number_of_rows == 0)
    expect_false(actual_number_of_columns == 0)
    expect_false(is.null(dim(out)))
    expect_equal(actual_number_of_columns, expected_number_of_columns)
  })

  test_that("the columns in the dataframe returned by load_MPAreports() correspond to the appropriate classes", {
    logger::log_threshold(logger::OFF)

    #  Get path to test directory.
    path_to_mpa <- get_test_mpa_report_dir("valid")

    # Load the MPA-style reports.
    out <- SPARKI::load_MPAreports(path_to_mpa)

    # Carry out tests.
    expect_type(out[[SPARKI:::COLNAME_MPA_SAMPLE]], "character")
    expect_type(out[[SPARKI:::COLNAME_MPA_TAXON_LEAF]], "character")
    expect_type(out[[SPARKI:::COLNAME_MPA_RANK]], "character")
    expect_type(out[[SPARKI:::COLNAME_MPA_N_FRAG_CLADE]], "integer")
    expect_type(out[[SPARKI:::COLNAME_MPA_DOMAIN]], "character")
    expect_type(out[[SPARKI:::COLNAME_MPA_KINGDOM]], "character")
    expect_type(out[[SPARKI:::COLNAME_MPA_PHYLUM]], "character")
    expect_type(out[[SPARKI:::COLNAME_MPA_CLASS]], "character")
    expect_type(out[[SPARKI:::COLNAME_MPA_ORDER]], "character")
    expect_type(out[[SPARKI:::COLNAME_MPA_FAMILY]], "character")
    expect_type(out[[SPARKI:::COLNAME_MPA_GENUS]], "character")
    expect_type(out[[SPARKI:::COLNAME_MPA_SPECIES]], "character")
  })

  test_that("the columns in the dataframe returned by load_MPAreports() have the right names", {
    logger::log_threshold(logger::OFF)

    #  Get path to test directory.
    path_to_mpa <- get_test_mpa_report_dir("valid")

    # Load the MPA-style reports and get actual column names.
    out <- SPARKI::load_MPAreports(path_to_mpa)
    actual_colnames <- colnames(out)

    # Get expected column names.
    expected_colnames <- c(
      SPARKI:::COLNAME_MPA_SAMPLE,
      SPARKI:::COLNAME_MPA_TAXON_LEAF,
      SPARKI:::COLNAME_MPA_RANK,
      SPARKI:::COLNAME_MPA_N_FRAG_CLADE,
      SPARKI:::COLNAME_MPA_DOMAIN,
      SPARKI:::COLNAME_MPA_KINGDOM,
      SPARKI:::COLNAME_MPA_PHYLUM,
      SPARKI:::COLNAME_MPA_CLASS,
      SPARKI:::COLNAME_MPA_ORDER,
      SPARKI:::COLNAME_MPA_FAMILY,
      SPARKI:::COLNAME_MPA_GENUS,
      SPARKI:::COLNAME_MPA_SPECIES
    )

    # Carry out tests.
    expect_equal(actual_colnames, expected_colnames)
  })

  test_that("load_MPAreports() does not fail if an MPA-style report is empty", {
    logger::log_threshold(logger::OFF)

    # Get path to test directory and define the name of the sample whose file
    # known to be empty.
    path_to_mpa <- get_test_mpa_report_dir("empty")
    empty_sample <- "sample2_empty"

    # First try to run load_MPAreports() to ensure no errors will be thrown.
    expect_no_error(SPARKI::load_MPAreports(path_to_mpa))

    # Now run load_MPAreports() again, storing the resulting dataframe for
    # further inspection.
    out <- SPARKI::load_MPAreports(path_to_mpa)

    # Assert that the sample whose file was empty is not present in the output.
    expect_false(empty_sample %in% out[[SPARKI:::COLNAME_MPA_SAMPLE]])
  })

  test_that("load_MPAreports() generates the rank columns (from domain to species) correctly", {
    logger::log_threshold(logger::OFF)

    # Get path to test directory.
    path_to_mpa <- get_test_mpa_report_dir("valid")

    # Load the MPA-style reports and obtain the actual numbers of rows and columns.
    out <- SPARKI::load_MPAreports(path_to_mpa)

    # Define columns that need to be checked.
    columns_to_check <- c(
      SPARKI:::COLNAME_MPA_DOMAIN,
      SPARKI:::COLNAME_MPA_KINGDOM,
      SPARKI:::COLNAME_MPA_PHYLUM,
      SPARKI:::COLNAME_MPA_CLASS,
      SPARKI:::COLNAME_MPA_ORDER,
      SPARKI:::COLNAME_MPA_FAMILY,
      SPARKI:::COLNAME_MPA_GENUS,
      SPARKI:::COLNAME_MPA_SPECIES
    )

    # For each column specified above...
    for (column in columns_to_check) {
      # Get the number of NAs in the column.
      n_NAs_in_column <- out[[column]] |> is.na() |> sum()

      # Get the number of elements in the column (i.e. column length).
      n_elements_in_column <- out[[column]] |> length()

      # The number of NAs should be less than the number of elements in the column.
      # If number of NAs = number of elements, then something has gone wrong.
      expect_lt(n_NAs_in_column, n_elements_in_column)
    }
  })

  test_that("the helper function check_mpa_lines() behaves as expected", {
    logger::log_threshold(logger::OFF)

    # Create mock MPA-style report with a few possibilities of taxon hierarachies.
    # Sample column.
    samples <- rep("sample", times = 7)
    # Taxon hierarchy column.
    taxon_hierarchies <- c(
      "d__Domain|k__Kingdom|p__Phylum|c__Class|o__Order|f__Family|g__Genus|s__Species", # Complete case.
      "d__Domain|p__Phylum|c__Class|o__Order|f__Family|g__Genus|s__Species", # Missing kingdom.
      "s__Species", # Missing everything but species.
      "d__Domain|p__Phylum|g__Genus|s__Species", # Missing kingdom, class, order, and family.
      "d__Domain|s__Species", # Missing everything but domain and species.
      "d__Domain|k__Kingdom|p__Phylum|c__Class|o__Order|f__Family", # Missing genus and species.
      "d__Domain" # Missing everything but domain.
    )
    # Column with clade-level fragments.
    n_frag_clade <- c(10, 2, 5, 10, 4, 12, 20)
    # Rank column.
    ranks <- c("S", "S", "S", "S", "S", "F", "D")
    # Create mock dataframe.
    mpa_reports <- data.frame(samples, taxon_hierarchies, n_frag_clade, ranks)
    colnames(mpa_reports) <- c(
      SPARKI:::COLNAME_MPA_SAMPLE,
      SPARKI:::COLNAME_MPA_TAXON_HIERARCHY,
      SPARKI:::COLNAME_MPA_N_FRAG_CLADE,
      SPARKI:::COLNAME_MPA_RANK
    )

    # Apply check_mpa_lines() to mock MPA-style report.
    checked_mpa_reports <- SPARKI:::check_mpa_lines(mpa_reports)

    # Define expected hierarchies.
    expected_hierachies <- c(
      "d__Domain|k__Kingdom|p__Phylum|c__Class|o__Order|f__Family|g__Genus|s__Species",
      "d__Domain|k__NA|p__Phylum|c__Class|o__Order|f__Family|g__Genus|s__Species",
      "d__NA|k__NA|p__NA|c__NA|o__NA|f__NA|g__NA|s__Species",
      "d__Domain|k__NA|p__Phylum|c__NA|o__NA|f__NA|g__Genus|s__Species",
      "d__Domain|k__NA|p__NA|c__NA|o__NA|f__NA|g__NA|s__Species",
      "d__Domain|k__Kingdom|p__Phylum|c__Class|o__Order|f__Family|g__NA|s__NA",
      "d__Domain|k__NA|p__NA|c__NA|o__NA|f__NA|g__NA|s__NA"
    )
    # Define actual hierarchies.
    actual_hierachies <- checked_mpa_reports[[SPARKI:::COLNAME_MPA_TAXON_HIERARCHY]]

    # Carry out test.
    expect_identical(actual_hierachies, expected_hierachies)
  })
})

describe("load_STDreports()", {
  test_that("load_STDreports() returns a dataframe with the appropriate dimensions", {
    logger::log_threshold(logger::OFF)

    # Define the expected number of columns.
    expected_number_of_columns <- 9

    #  Get path to test directory.
    path_to_std <- get_test_std_report_dir("valid")

    # Load the MPA-style reports and obtain the actual numbers of rows and columns.
    out <- SPARKI::load_STDreports(path_to_std)
    actual_number_of_rows <- nrow(out)
    actual_number_of_columns <- ncol(out)

    # Carry out tests.
    expect_s3_class(out, "data.frame")
    expect_false(actual_number_of_rows == 0)
    expect_false(actual_number_of_columns == 0)
    expect_false(is.null(dim(out)))
    expect_equal(actual_number_of_columns, expected_number_of_columns)
  })

  test_that("the columns in the dataframe returned by load_STDreports() correspond to the appropriate classes", {
    logger::log_threshold(logger::OFF)

    # Get path to test directory.
    path_to_std <- get_test_std_report_dir("valid")

    # Load the standard reports.
    out <- SPARKI::load_STDreports(path_to_std)

    # Carry out tests.
    expect_type(out[[SPARKI:::COLNAME_STD_SAMPLE]], "character")
    expect_type(out[[SPARKI:::COLNAME_STD_PCT_FRAG_CLADE]], "double")
    expect_type(out[[SPARKI:::COLNAME_STD_N_FRAG_CLADE]], "integer")
    expect_type(out[[SPARKI:::COLNAME_STD_N_FRAG_TAXON]], "integer")
    expect_type(out[[SPARKI:::COLNAME_STD_MINIMISERS]], "integer")
    expect_type(out[[SPARKI:::COLNAME_STD_UNIQ_MINIMISERS]], "integer")
    expect_type(out[[SPARKI:::COLNAME_STD_RANK]], "character")
    expect_type(out[[SPARKI:::COLNAME_STD_NCBI_ID]], "character")
    expect_type(out[[SPARKI:::COLNAME_STD_TAXON]], "character")
  })

  test_that("the columns in the dataframe returned by load_STDreports() have the right names", {
    logger::log_threshold(logger::OFF)

    # Get path to test directory.
    path_to_std <- get_test_std_report_dir("valid")

    # Load the standard reports and get actual column names.
    out <- SPARKI::load_STDreports(path_to_std)
    actual_colnames <- colnames(out)

    # Get expected column names.
    expected_colnames <- c(
      SPARKI:::COLNAME_STD_SAMPLE,
      SPARKI:::COLNAME_STD_PCT_FRAG_CLADE,
      SPARKI:::COLNAME_STD_N_FRAG_CLADE,
      SPARKI:::COLNAME_STD_N_FRAG_TAXON,
      SPARKI:::COLNAME_STD_MINIMISERS,
      SPARKI:::COLNAME_STD_UNIQ_MINIMISERS,
      SPARKI:::COLNAME_STD_RANK,
      SPARKI:::COLNAME_STD_NCBI_ID,
      SPARKI:::COLNAME_STD_TAXON
    )

    # Carry out tests.
    expect_equal(actual_colnames, expected_colnames)
  })

  test_that("load_STDreports() does not fail if a standard report is empty", {
    logger::log_threshold(logger::OFF)

    # Get path to test directory and define the name of the sample whose file
    # known to be empty.
    path_to_std <- get_test_std_report_dir("empty")
    empty_sample <- "sample2_empty"

    # First try to run load_STDreports() to ensure no errors will be thrown.
    expect_no_error(SPARKI::load_STDreports(path_to_std))

    # Now run load_STDreports() again, storing the resulting dataframe for
    # further inspection.
    out <- SPARKI::load_STDreports(path_to_std)

    # Assert that the sample whose file was empty is not present in the output.
    expect_false(empty_sample %in% out[[SPARKI:::COLNAME_STD_SAMPLE]])
  })

  test_that("load_STDreports() fails if a report does not contain 'D' as the domain identifier", {
    logger::log_threshold(logger::OFF)

    # Get path to test directory.
    path_to_std <- get_test_std_report_dir("db_incompatible")

    # Load the standard reports.
    expect_error(SPARKI::load_STDreports(path_to_std))
  })
})
