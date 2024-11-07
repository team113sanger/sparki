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
})
