describe("load_MPAreports()", {
  test_that("load_MPAreports() returns a dataframe with the appropriate dimensions", {
    #  Define the expected number of columns.
    expected_number_of_columns <- 12

    # Get path to test directory.
    path_to_mpa <- get_test_mpa_report_dir("usual")

    # Load the MPA-style reports and obtain the actual numbers of rows and columns.
    out <- SPARKI::load_MPAreports(path_to_mpa, verbose = FALSE)
    actual_number_of_rows <- nrow(out)
    actual_number_of_columns <- ncol(out)

    # Carry out tests.
    expect_false(actual_number_of_rows == 0)
    expect_false(actual_number_of_columns == 0)
    expect_false(is.null(dim(out)))
    expect_equal(actual_number_of_columns, expected_number_of_columns)
  })

  test_that("the columns in the dataframe returned by load_MPAreports() correspond to the appropriate classes", {
    # Get path to test directory.
    path_to_mpa <- get_test_mpa_report_dir("usual")

    # Load the MPA-style reports.
    out <- SPARKI::load_MPAreports(path_to_mpa, verbose = FALSE)

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

  test_that("the columns in the dataframe returned by load_MPAreports() have the right names", {
    # Get path to test directory.
    path_to_mpa <- get_test_mpa_report_dir("usual")

    # Load the MPA-style reports and get actual column names.
    out <- SPARKI::load_MPAreports(path_to_mpa, verbose = FALSE)
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

  test_that("load_MPAreports() does not fail if an MPA-style report is empty", {
    # Get path to test directory.
    path_to_mpa <- get_test_mpa_report_dir("empty")
    path_to_mpa <- paste0("'", path_to_mpa, "'")

    # Build argument list.
    args <- c(
      "-e",
      glue::glue("devtools::load_all(quiet = TRUE); SPARKI::load_MPAreports({path_to_mpa}, verbose = TRUE)")
    )

    # Run the command
    result <- processx::run(
      command = Sys.which("Rscript"),
      args = args,
      error_on_status = FALSE,
      echo = FALSE,
      stdout_line_callback = NULL,
      stderr_line_callback = NULL
    )

    # Define the name of the sample whose file is empty.
    empty_sample <- "sample2_empty"

    # Define the expected exit status and info/warning messages.
    expected_exit_code <- 0
    expected_info_in_stderr <- "MPA-style reports loaded successfully!"
    expected_warning_in_stderr <- "is empty, so this sample will not be included in the SPARKI analysis."

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

    # Assert that the sample whose file was empty is not present in the output.
    expect_false(
      grepl(empty_sample, actual_stdout, fixed = TRUE),
      info = sprintf("Standard out was:\n%s\nStandard error was:\n%s", actual_stdout, actual_stderr)
    )
  })
})

describe("load_STDreports()", {
  test_that("load_STDreports() returns a dataframe with the appropriate dimensions", {
    # Define the expected number of columns.
    expected_number_of_columns <- 9

    # Get path to test directory.
    path_to_std <- get_test_std_report_dir("usual")

    # Load the MPA-style reports and obtain the actual numbers of rows and columns.
    out <- SPARKI::load_STDreports(path_to_std, verbose = FALSE)
    actual_number_of_rows <- nrow(out)
    actual_number_of_columns <- ncol(out)

    # Carry out tests.
    expect_false(actual_number_of_rows == 0)
    expect_false(actual_number_of_columns == 0)
    expect_false(is.null(dim(out)))
    expect_equal(actual_number_of_columns, expected_number_of_columns)
  })

  test_that("the columns in the dataframe returned by load_STDreports() correspond to the appropriate classes", {
    # Get path to test directory.
    path_to_std <- get_test_std_report_dir("usual")

    # Load the MPA-style reports.
    out <- SPARKI::load_STDreports(path_to_std, verbose = FALSE)

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

  test_that("the columns in the dataframe returned by load_STDreports() have the right names", {
    # Get path to test directory.
    path_to_std <- get_test_std_report_dir("usual")

    # Load the MPA-style reports and get actual column names.
    out <- SPARKI::load_STDreports(path_to_std, verbose = FALSE)
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

  test_that("load_STDreports() does not fail if a standard report is empty", {
    # Get path to test directory.
    path_to_std <- get_test_std_report_dir("empty")
    path_to_std <- paste0("'", path_to_std, "'")

    # Build argument list.
    args <- c(
      "-e",
      glue::glue("devtools::load_all(quiet = TRUE); SPARKI::load_STDreports({path_to_std}, verbose = TRUE)")
    )

    # Run the command.
    result <- processx::run(
      command = Sys.which("Rscript"),
      args = args,
      error_on_status = FALSE,
      echo = FALSE,
      stdout_line_callback = NULL,
      stderr_line_callback = NULL
    )

    # Define the name of the sample whose file is empty.
    empty_sample <- "sample2_empty"

    # Define the expected exit status and info/warning messages.
    expected_exit_code <- 0
    expected_info_in_stderr <- "Standard reports loaded successfully!"
    expected_warning_in_stderr <- "is empty, so this sample will not be included in the SPARKI analysis."

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

    # Assert that the sample whose file was empty is not present in the output.
    expect_false(
      grepl(empty_sample, actual_stdout, fixed = TRUE),
      info = sprintf("Standard out was:\n%s\nStandard error was:\n%s", actual_stdout, actual_stderr)
    )
  })
})



# testthat the values under domain/kingdom/.../genus/species are not just a big block of NAs!
