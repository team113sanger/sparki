describe("load_MPAreports()", {
  test_that("load_MPAreports() returns a dataframe with the appropriate dimensions", {
    #  Define the expected number of columns.
    expected_number_of_columns <- 12

    #  Get path to test directory.
    path_to_mpa <- get_test_mpa_report_dir("usual")

    # Load the MPA-style reports and obtain the actual numbers of rows and columns.
    out <- SPARKI::load_MPAreports(path_to_mpa, verbose = FALSE)
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
    #  Get path to test directory.
    path_to_mpa <- get_test_mpa_report_dir("usual")

    # Load the MPA-style reports.
    out <- SPARKI::load_MPAreports(path_to_mpa, verbose = FALSE)

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
    #  Get path to test directory.
    path_to_mpa <- get_test_mpa_report_dir("usual")

    # Load the MPA-style reports and get actual column names.
    out <- SPARKI::load_MPAreports(path_to_mpa, verbose = FALSE)
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
    #  Get path to test directory.
    path_to_mpa <- glue::single_quote(get_test_mpa_report_dir("empty"))

    # Build argument list.
    args <- c("-e", glue::glue("devtools::load_all(quiet = TRUE); SPARKI::load_MPAreports({path_to_mpa}, verbose = TRUE)"))

    # Run the command
    result <- processx::run(
      command = Sys.which("Rscript"),
      args = args,
      error_on_status = FALSE,
      echo = FALSE,
      stdout_line_callback = NULL,
      stderr_line_callback = NULL
    )

    #  Define the name of the sample whose file is empty.
    empty_sample <- "sample2_empty"

    #  Define the expected exit status and info/warning messages.
    expected_exit_code <- 0
    expected_info_in_stderr <- SPARKI:::LOAD_MPA_INFO_SUCCESS
    expected_warning_in_stderr <- SPARKI:::CHECK_EMPTY_FILE_WARNING

    # Get actual exit code and standard out/error.
    actual_exit_code <- result$status
    actual_stdout <- result$stdout
    actual_stderr <- result$stderr

    # Assert that the exit code is zero.
    expect_equal(
      actual_exit_code,
      expected_exit_code,
      info = sprintf("Command failed with status %s", actual_exit_code)
    )

    # Assert that stderr contains the expected log info/warning messages.
    expect_true(
      grepl(expected_info_in_stderr, actual_stderr, fixed = TRUE),
      info = sprintf("Standard out was:\n%s\nStandard error was:\n%s", actual_stdout, actual_stderr)
    )
    expect_true(
      grepl(expected_warning_in_stderr, actual_stderr, fixed = TRUE),
      info = sprintf("Standard out was: %s\nStandard error was: %s", actual_stdout, actual_stderr)
    )

    #  Assert that the sample whose file was empty is not present in the output.
    expect_false(
      grepl(empty_sample, actual_stdout, fixed = TRUE),
      info = sprintf("Standard out was:\n%s\nStandard error was:\n%s", actual_stdout, actual_stderr)
    )
  })

  test_that("load_MPAreports() generates the rank columns (from domain to species) correctly", {
    #  Get path to test directory.
    path_to_mpa <- get_test_mpa_report_dir("usual")

    # Load the MPA-style reports and obtain the actual numbers of rows and columns.
    out <- SPARKI::load_MPAreports(path_to_mpa, verbose = FALSE)

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
})

describe("load_STDreports()", {
  test_that("load_STDreports() returns a dataframe with the appropriate dimensions", {
    # Define the expected number of columns.
    expected_number_of_columns <- 9

    #  Get path to test directory.
    path_to_std <- get_test_std_report_dir("usual")

    # Load the MPA-style reports and obtain the actual numbers of rows and columns.
    out <- SPARKI::load_STDreports(path_to_std, verbose = FALSE)
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
    #  Get path to test directory.
    path_to_std <- get_test_std_report_dir("usual")

    # Load the MPA-style reports.
    out <- SPARKI::load_STDreports(path_to_std, verbose = FALSE)

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
    #  Get path to test directory.
    path_to_std <- get_test_std_report_dir("usual")

    # Load the MPA-style reports and get actual column names.
    out <- SPARKI::load_STDreports(path_to_std, verbose = FALSE)
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
    #  Get path to test directory.
    path_to_std <- get_test_std_report_dir("empty")
    path_to_std <- glue::single_quote(path_to_std)

    # Build argument list.
    args <- c("-e", glue::glue("devtools::load_all(quiet = TRUE); SPARKI::load_STDreports({path_to_std}, verbose = TRUE)"))

    # Run the command.
    result <- processx::run(
      command = Sys.which("Rscript"),
      args = args,
      error_on_status = FALSE,
      echo = FALSE,
      stdout_line_callback = NULL,
      stderr_line_callback = NULL
    )

    #  Define the name of the sample whose file is empty.
    empty_sample <- "sample2_empty"

    #  Define the expected exit status and info/warning messages.
    expected_exit_code <- 0
    expected_info_in_stderr <- SPARKI:::LOAD_STD_INFO_SUCCESS
    expected_warning_in_stderr <- SPARKI:::CHECK_EMPTY_FILE_WARNING

    # Get actual exit code and standard out/error.
    actual_exit_code <- result$status
    actual_stdout <- result$stdout
    actual_stderr <- result$stderr

    # Assert that the exit code is zero.
    expect_equal(
      actual_exit_code,
      expected_exit_code,
      info = sprintf("Command failed with status %s", actual_exit_code)
    )

    # Assert that stderr contains the expected log info/warning messages.
    expect_true(
      grepl(expected_info_in_stderr, actual_stderr, fixed = TRUE),
      info = sprintf("Standard out was:\n%s\nStandard error was:\n%s", actual_stdout, actual_stderr)
    )
    expect_true(
      grepl(expected_warning_in_stderr, actual_stderr, fixed = TRUE),
      info = sprintf("Standard out was: %s\nStandard error was: %s", actual_stdout, actual_stderr)
    )

    #  Assert that the sample whose file was empty is not present in the output.
    expect_false(
      grepl(empty_sample, actual_stdout, fixed = TRUE),
      info = sprintf("Standard out was:\n%s\nStandard error was:\n%s", actual_stdout, actual_stderr)
    )
  })
})
