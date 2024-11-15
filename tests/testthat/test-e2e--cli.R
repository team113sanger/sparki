describe("CLI integration", {
  test_that("version option works", {
    # Build the command and get the expected version
    rscript <- Sys.which("Rscript")
    args <- c(
      "-e", "devtools::load_all(quiet = TRUE); SPARKI::cli()",
      "--version"
    )
    pkg_version <- as.character(utils::packageVersion("SPARKI"))

    # Run the command
    result <- processx::run(
      command = rscript,
      args = args,
      error_on_status = FALSE,
      echo = FALSE,
      stdout_line_callback = NULL,
      stderr_line_callback = NULL
    )


    # Assert that the exit code is zero
    expect_equal(result$status, 0,
      info = sprintf("Command failed with status %s", result$status)
    )

    # Assert that the output contains the expected version
    expect_true(
      grepl(pkg_version, result$stdout, fixed = TRUE),
      info = sprintf(
        "Expected version %s in output, but got:\n%s",
        pkg_version,
        result$stdout
      )
    )
  })

  test_that("help option works", {
    # Build the command and get the expected version
    rscript <- Sys.which("Rscript")
    args <- c(
      "-e", "devtools::load_all(quiet = TRUE); SPARKI::cli()",
      "--help"
    )
    expected_pkg_description <- SPARKI:::CLI_DESCRIPTION
    expected_program_name <- SPARKI:::CLI_PROGRAM_NAME
    expected_exit_code <- 0

    # Run the command
    result <- processx::run(
      command = rscript,
      args = args,
      error_on_status = FALSE,
      echo = FALSE,
      stdout_line_callback = NULL,
      stderr_line_callback = NULL
    )
    actual_exit_code <- result$status
    actual_output <- result$stdout

    # Assert that the exit code is zero
    expect_equal(actual_exit_code, expected_exit_code,
      info = sprintf("Command failed with status %s", actual_exit_code)
    )

    # Assert that the output contains the expected package description.
    expect_true(grepl(expected_pkg_description, actual_output, fixed = TRUE))

    # Assert that the output contains the expected program name.
    expect_true(grepl(expected_program_name, actual_output, fixed = TRUE))
  })


test_that("main - all options", {
    # Build the command and get the expected version
    rscript <- Sys.which("Rscript")
    args <- c(
      "-e", "devtools::load_all(quiet = TRUE); SPARKI::cli()",
      "--help"
    )
    expected_exit_code <- 0
    expect_true(FALSE) # Deliberately fail the test

    # # Run the command
    # result <- processx::run(
    #   command = rscript,
    #   args = args,
    #   error_on_status = FALSE,
    #   echo = FALSE,
    #   stdout_line_callback = NULL,
    #   stderr_line_callback = NULL
    # )
    # actual_exit_code <- result$status
    # actual_output <- result$stdout

    # # Assert that the exit code is zero
    # expect_equal(actual_exit_code, expected_exit_code,
    #   info = sprintf("Command failed with status %s", actual_exit_code)
    # )
  })

test_that("main - required options", {
    # Setup
    std_directory <- get_test_std_report_dir("usual")
    option_std_directory <- glue::glue("--std-reports '{std_directory}'")

    mpa_directory <- get_test_mpa_report_dir("usual")
    option_mpa_directory <- glue::glue("--mpa-reports '{mpa_directory}'")

    organism <- "Homo sapiens"
    option_organism <- glue::glue("--organism '{organism}'")

    option_prefix <- "--prefix 'foo'"

    reference <- get_test_reference()
    option_reference <- glue::glue("--reference '{reference}xyz'")

    domain <- "Viruses,Bacteria"
    option_domain <- glue::glue("--domain '{domain}'")

    outdir <- get_local_tmp_dir()
    option_outdir <- glue::glue("--outdir '{outdir}'")


    # Build the command and get the expected version
    rscript <- Sys.which("Rscript")
    args <- c(
      "-e", "'devtools::load_all(quiet = TRUE); SPARKI::cli()'",
      option_std_directory,
      option_mpa_directory,
      option_reference,
      option_prefix,
      option_organism,
      option_domain,
      option_outdir
    )
    expected_exit_code <- 512

    # Run the command
    result <- processx::run(
      command = rscript,
      args = args,
      error_on_status = FALSE,
      echo = FALSE,
      stdout_line_callback = NULL,
      stderr_line_callback = NULL
    )
    actual_exit_code <- result$status


    # Assert that the exit code is zero
    expect_equal(actual_exit_code, expected_exit_code,
      info = sprintf("Command failed with status %s and with args:\n%s\nand stdout:\n%s \nand stderr:\n%s", actual_exit_code, paste(args, collapse = " "), result$stdout, result$stderr)
    )
  })
})
