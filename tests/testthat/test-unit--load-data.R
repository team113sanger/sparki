describe("load_MPAreports()", {
  test_that("returns a dataframe with the appropriate dimensions", {
    #  Define the expected number of columns.
    expected_number_of_columns <- 12

    # Get path to test directory.
    path_to_mpa <- get_test_mpa_report_dir("usual")

    # Load the MPA-style reports and obtain the actual numbers of rows and columns.
    out <- SPARKI::load_MPAreports(path_to_mpa)
    actual_number_of_rows <- nrow(out)
    actual_number_of_columns <- ncol(out)

    # Carry out tests.
    expect_false(actual_number_of_rows == 0)
    expect_false(actual_number_of_columns == 0)
    expect_false(is.null(dim(out)))
    expect_equal(actual_number_of_columns, expected_number_of_columns)
  })

  test_that("the columns in the dataframe it returns correspond to the appropriate classes", {
    # Get path to test directory.
    path_to_mpa <- get_test_mpa_report_dir("usual")

    # Load the MPA-style reports.
    out <- SPARKI::load_MPAreports(path_to_mpa)

    # Carry out tests.
    expect_true(is.character(out[[COLNAME_MPA_SAMPLE]]))
    expect_true(is.character(out[[COLNAME_MPA_TAXON_LEAF]]))
    expect_true(is.character(out[[COLNAME_MPA_RANK]]))
    expect_true(is.double(out[[COLNAME_MPA_N_FRAG_CLADE]]))
    expect_true(is.character(out[[COLNAME_MPA_DOMAIN]]))
    expect_true(is.character(out[[COLNAME_MPA_KINGDOM]]))
    expect_true(is.character(out[[COLNAME_MPA_PHYLUM]]))
    expect_true(is.character(out[[COLNAME_MPA_CLASS]]))
    expect_true(is.character(out[[COLNAME_MPA_ORDER]]))
    expect_true(is.character(out[[COLNAME_MPA_FAMILY]]))
    expect_true(is.character(out[[COLNAME_MPA_GENUS]]))
    expect_true(is.character(out[[COLNAME_MPA_SPECIES]]))
  })

  test_that("the columns in the dataframe it returns have the right names", {
    # Get path to test directory.
    path_to_mpa <- get_test_mpa_report_dir("usual")

    # Load the MPA-style reports and get actual column names.
    out <- SPARKI::load_MPAreports(path_to_mpa)
    actual_colnames <- colnames(out)

    # Get expected column names.
    expected_colnames <- c(
      COLNAME_MPA_SAMPLE,
      COLNAME_MPA_TAXON_LEAF,
      COLNAME_MPA_RANK,
      COLNAME_MPA_N_FRAG_CLADE,
      COLNAME_MPA_DOMAIN,
      COLNAME_MPA_KINGDOM,
      COLNAME_MPA_PHYLUM,
      COLNAME_MPA_CLASS,
      COLNAME_MPA_ORDER,
      COLNAME_MPA_FAMILY,
      COLNAME_MPA_GENUS,
      COLNAME_MPA_SPECIES
    )

    # Carry out tests.
    expect_equal(actual_colnames, expected_colnames)
  })

  test_that("throws a warning if an MPA-style report is empty", {
    # Get path to test directory.
    path_to_mpa <- get_test_mpa_report_dir("empty")

    # Carry out test while trying to load an empty MPA-style report.
    expect_warning(SPARKI::load_MPAreports(path_to_mpa))
  })
})

describe("load_STDreports()", {
  test_that("returns a dataframe with the appropriate dimensions", {
    # Define the expected number of columns.
    expected_number_of_columns <- 9

    # Get path to test directory.
    path_to_std <- get_test_std_report_dir("usual")

    # Load the MPA-style reports and obtain the actual numbers of rows and columns.
    out <- SPARKI::load_STDreports(path_to_std)
    actual_number_of_rows <- nrow(out)
    actual_number_of_columns <- ncol(out)

    # Carry out tests.
    expect_false(actual_number_of_rows == 0)
    expect_false(actual_number_of_columns == 0)
    expect_false(is.null(dim(out)))
    expect_equal(actual_number_of_columns, expected_number_of_columns)
  })

  test_that("the columns in the dataframe it returns correspond to the appropriate classes", {
    # Get path to test directory.
    path_to_std <- get_test_std_report_dir("usual")

    # Load the MPA-style reports.
    out <- SPARKI::load_STDreports(path_to_std)

    # Carry out tests.
    expect_true(is.character(out[[COLNAME_STD_SAMPLE]]))
    expect_true(is.double(out[[COLNAME_STD_PCT_FRAG_CLADE]]))
    expect_true(is.double(out[[COLNAME_STD_N_FRAG_CLADE]]))
    expect_true(is.double(out[[COLNAME_STD_N_FRAG_TAXON]]))
    expect_true(is.double(out[[COLNAME_STD_MINIMISERS]]))
    expect_true(is.double(out[[COLNAME_STD_UNIQ_MINIMISERS]]))
    expect_true(is.character(out[[COLNAME_STD_RANK]]))
    expect_true(is.double(out[[COLNAME_STD_NCBI_ID]]))
    expect_true(is.character(out[[COLNAME_STD_TAXON]]))
  })

  test_that("the columns in the dataframe it returns have the right names", {
    # Get path to test directory.
    path_to_std <- get_test_std_report_dir("usual")

    # Load the MPA-style reports and get actual column names.
    out <- SPARKI::load_STDreports(path_to_std)
    actual_colnames <- colnames(out)

    # Get expected column names.
    expected_colnames <- c(
      COLNAME_STD_SAMPLE,
      COLNAME_STD_PCT_FRAG_CLADE,
      COLNAME_STD_N_FRAG_CLADE,
      COLNAME_STD_N_FRAG_TAXON,
      COLNAME_STD_MINIMISERS,
      COLNAME_STD_UNIQ_MINIMISERS,
      COLNAME_STD_RANK,
      COLNAME_STD_NCBI_ID,
      COLNAME_STD_TAXON
    )

    # Carry out tests.
    expect_equal(actual_colnames, expected_colnames)
  })

  test_that("throws a warning if a standard report is empty", {
    # Get path to test directory.
    path_to_std <- get_test_std_report_dir("empty")

    # Carry out test while trying to load an empty standard report.
    expect_warning(SPARKI::load_STDreports(path_to_std))
  })
})
