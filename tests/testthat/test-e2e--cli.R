describe("CLI integration", {
  test_that("version command outputs correct version and exits successfully", {
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
    expected_pkg_description <- CLI_DESCRIPTION
    expected_program_name <- CLI_PROGRAM_NAME
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
})
