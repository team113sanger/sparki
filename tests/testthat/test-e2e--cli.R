describe("CLI integration", {
  test_that("the --version option in the CLI works", {
    # Define the arguments that will be used with SPARKI's CLI.
    args <- c("-e", "devtools::load_all(quiet = TRUE); SPARKI::cli()", "--version")

    # Get the expected version.
    expected_version <- as.character(utils::packageVersion("SPARKI"))

    # Run the command
    result <- processx::run(
      command = Sys.which("Rscript"),
      args = args,
      error_on_status = FALSE,
      echo = FALSE,
      stdout_line_callback = NULL,
      stderr_line_callback = NULL
    )

    #  Define the expected exit code.
    expected_exit_code <- 0

    #  Get actual outputs of the run.
    actual_output <- result$stdout
    actual_exit_code <- result$status

    # Assert that the output contains the expected version
    expect_true(
      grepl(expected_version, actual_output, fixed = TRUE),
      info = sprintf("Expected version is %s, but got:\n%s", expected_version, actual_output)
    )

    # Assert that the exit code is zero
    expect_equal(
      actual_exit_code,
      expected_exit_code,
      info = sprintf("Command failed with status %s", actual_exit_code)
    )
  })

  test_that("the --help option in the CLI works", {
    # Define the arguments that will be used with SPARKI's CLI.
    args <- c("-e", "devtools::load_all(quiet = TRUE); SPARKI::cli()", "--help")

    # Run SPARKI's CLI with the options that were specified.
    result <- processx::run(
      command = Sys.which("Rscript"),
      args = args,
      error_on_status = FALSE,
      echo = FALSE,
      stdout_line_callback = NULL,
      stderr_line_callback = NULL
    )

    #  Define the expected program name and description, as well as
    # the expected exit status code.
    expected_program_name <- SPARKI:::CLI_PROGRAM_NAME
    expected_program_description <- SPARKI:::CLI_DESCRIPTION
    expected_exit_code <- 0

    #  Get actual outputs of the run.
    actual_output <- result$stdout
    actual_exit_code <- result$status

    # Assert that the output contains the expected program name.
    expect_true(grepl(expected_program_name, actual_output, fixed = TRUE))

    # Assert that the output contains the expected program description.
    expect_true(grepl(expected_program_description, actual_output, fixed = TRUE))

    # Assert that the exit code is zero.
    expect_equal(
      actual_exit_code,
      expected_exit_code,
      info = sprintf("Command failed with status %s", actual_exit_code)
    )
  })


  # test_that("CLI works when all arguments, both required and optional, are provided", {
  #   #  Define input arguments for SPARKI.
  #   arg_std_directory <- glue::single_quote(get_test_std_report_dir("usual"))
  #   arg_mpa_directory <- glue::single_quote(get_test_mpa_report_dir("usual"))
  #   arg_organism <- glue::single_quote("Homo sapiens")
  #   arg_reference <- glue::single_quote(get_test_reference())
  #   arg_domain <- "Viruses"
  #   arg_outdir <- glue::single_quote(get_local_tmp_dir())
  #   arg_metadata <- glue::single_quote(get_test_metadata())
  #   arg_metadata_sample_col <- "sample"
  #   arg_metadata_columns <- "type,date,status"
  #   arg_prefix <- "test"
  #   arg_samples_to_remove <- glue::single_quote(get_test_samples_to_remove())

  #   # Build argument list.
  #   args <- c(
  #     "-e", "devtools::load_all(quiet = TRUE); SPARKI::cli()",
  #     glue::glue("--std-reports {arg_std_directory}"),
  #     glue::glue("--mpa-reports {arg_mpa_directory}"),
  #     glue::glue("--organism {arg_organism}"),
  #     glue::glue("--reference {arg_reference}"),
  #     glue::glue("--domain {arg_domain}"),
  #     glue::glue("--outdir {arg_outdir}"),
  #     glue::glue("--metadata {arg_metadata}"),
  #     glue::glue("--sample-col {arg_metadata_sample_col}"),
  #     glue::glue("--columns {arg_metadata_sample_col}"),
  #     glue::glue("--prefix {arg_prefix}"),
  #     glue::glue("--samples-to-remove {arg_samples_to_remove}"),
  #     glue::glue("--verbose"),
  #     glue::glue("--include-eukaryotes"),
  #     glue::glue("--include-sample-names")
  #   )

  #   # Ensure outdir is empty before the program is run.
  #   initial_outdir_size <- length(list.files(arg_outdir))
  #   initial_outdir_contents <- list.files(arg_outdir)
  #   expect_equal(initial_outdir_size, 0)

  #   # Run the command.
  #   result <- processx::run(
  #     command = Sys.which("Rscript"),
  #     args = args,
  #     error_on_status = FALSE,
  #     echo = TRUE,
  #     stdout_line_callback = NULL,
  #     stderr_line_callback = NULL
  #   )

  #   # Define expected exit status.
  #   expected_exit_code <- 0

  #   # Get actual exit status.
  #   actual_exit_code <- result$status

  #   # Assert that the exit code is zero.
  #   expect_equal(
  #     actual_exit_code,
  #     expected_exit_code,
  #     info = sprintf("Stdout\n%s\nStderr\n%s", result$stdout, result$stderr)
  #   )

  #   # Assert that the outdir size is non-zero after the program is run.
  #   final_outdir_size <- length(list.files(arg_outdir))
  #   final_outdir_contents <- list.files(arg_outdir)
  #   #expect_gt(final_outdir_size, 0)

  #   # a <- sprintf("Initial size was %s, final size was %s. The initial contents were %s, the final contents were %s", initial_outdir_size, final_outdir_size, initial_outdir_contents, final_outdir_contents)
  #   # write.csv(a, "/lustre/scratch126/casm/team113da/users/jb62/projects/sparki/testing.csv")
  # })

  # test_that("CLI works when all required arguments are provided", {
  #   #  Define input arguments for SPARKI.
  #   std_directory <- get_test_std_report_dir("usual")
  #   mpa_directory <- get_test_mpa_report_dir("usual")
  #   organism <- "'Homo sapiens'"
  #   reference <- get_test_reference()
  #   domain <- "Viruses"
  #   outdir <- get_local_tmp_dir()

  #   # Ensure outdir is empty before the program is run.
  #   initial_outdir_size <- length(list.files(outdir))
  #   expect_equal(initial_outdir_size, 0)

  #   # Build argument list.
  #   args <- c(
  #     "-e", "'devtools::load_all(quiet = TRUE); SPARKI::cli()'",
  #     glue::glue("--std-reports {std_directory}"),
  #     glue::glue("--mpa-reports {mpa_directory}"),
  #     glue::glue("--organism {organism}"),
  #     glue::glue("--reference {reference}"),
  #     glue::glue("--domain {domain}"),
  #     glue::glue("--outdir {outdir}")
  #   )

  #   # Run the command
  #   result <- processx::run(
  #     command = Sys.which("Rscript"),
  #     args = args,
  #     error_on_status = FALSE,
  #     echo = FALSE,
  #     stdout_line_callback = NULL,
  #     stderr_line_callback = NULL
  #   )

  #   #  Define expected exit status.
  #   expected_exit_code <- 0

  #   # Get actual exit status.
  #   actual_exit_code <- result$status

  #   # Assert that the exit code is zero
  #   expect_equal(
  #     actual_exit_code,
  #     expected_exit_code,
  #     info = sprintf("Command failed with status %s", actual_exit_code)
  #   )

  #   # Assert that the outdir size is non-zero after the program is run.
  #   final_outdir_size <- length(list.files(outdir))
  #   #expect_gt(final_outdir_size, initial_outdir_size)
  # })

  test_that("CLI fails when the required argument --std-reports is not provided", {
    #  Define input arguments for SPARKI.
    std_directory <- NA
    mpa_directory <- get_test_mpa_report_dir("usual")
    organism <- "'Homo sapiens'"
    reference <- get_test_reference()
    domain <- "Viruses"
    outdir <- get_local_tmp_dir()

    # Build argument list.
    args <- c(
      "-e", "devtools::load_all(quiet = TRUE); SPARKI::cli()",
      glue::glue("--std-reports {std_directory}"),
      glue::glue("--mpa-reports {mpa_directory}"),
      glue::glue("--organism {organism}"),
      glue::glue("--reference {reference}"),
      glue::glue("--domain {domain}"),
      glue::glue("--outdir {outdir}")
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

    # Get actual exit status.
    actual_exit_code <- result$status

    # Assert that the exit code is not zero.
    expect_gt(actual_exit_code, 0)
  })

  test_that("CLI fails when the required argument --mpa-reports is not provided", {
    #  Define input arguments for SPARKI.
    std_directory <- get_test_std_report_dir("usual")
    mpa_directory <- NA
    organism <- "'Homo sapiens'"
    reference <- get_test_reference()
    domain <- "Viruses"
    outdir <- get_local_tmp_dir()

    # Build argument list.
    args <- c(
      "-e", "devtools::load_all(quiet = TRUE); SPARKI::cli()",
      glue::glue("--std-reports {std_directory}"),
      glue::glue("--mpa-reports {mpa_directory}"),
      glue::glue("--organism {organism}"),
      glue::glue("--reference {reference}"),
      glue::glue("--domain {domain}"),
      glue::glue("--outdir {outdir}")
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

    # Get actual exit status.
    actual_exit_code <- result$status

    # Assert that the exit code is not zero.
    expect_gt(actual_exit_code, 0)
  })

  test_that("CLI fails when the required argument --organism is not provided", {
    #  Define input arguments for SPARKI.
    std_directory <- get_test_std_report_dir("usual")
    mpa_directory <- get_test_mpa_report_dir("usual")
    organism <- NA
    reference <- get_test_reference()
    domain <- "Viruses"
    outdir <- get_local_tmp_dir()

    # Build argument list.
    args <- c(
      "-e", "devtools::load_all(quiet = TRUE); SPARKI::cli()",
      glue::glue("--std-reports {std_directory}"),
      glue::glue("--mpa-reports {mpa_directory}"),
      glue::glue("--organism {organism}"),
      glue::glue("--reference {reference}"),
      glue::glue("--domain {domain}"),
      glue::glue("--outdir {outdir}")
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

    # Get actual exit status.
    actual_exit_code <- result$status

    # Assert that the exit code is not zero.
    expect_gt(actual_exit_code, 0)
  })

  test_that("CLI fails when the required argument --reference is not provided", {
    #  Define input arguments for SPARKI.
    std_directory <- get_test_std_report_dir("usual")
    mpa_directory <- get_test_mpa_report_dir("usual")
    organism <- "'Homo sapiens'"
    reference <- NA
    domain <- "Viruses"
    outdir <- get_local_tmp_dir()

    # Build argument list.
    args <- c(
      "-e", "devtools::load_all(quiet = TRUE); SPARKI::cli()",
      glue::glue("--std-reports {std_directory}"),
      glue::glue("--mpa-reports {mpa_directory}"),
      glue::glue("--organism {organism}"),
      glue::glue("--reference {reference}"),
      glue::glue("--domain {domain}"),
      glue::glue("--outdir {outdir}")
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

    # Get actual exit status.
    actual_exit_code <- result$status

    # Assert that the exit code is not zero.
    expect_gt(actual_exit_code, 0)
  })

  test_that("CLI fails when the required argument --domain is not provided", {
    #  Define input arguments for SPARKI.
    std_directory <- get_test_std_report_dir("usual")
    mpa_directory <- get_test_mpa_report_dir("usual")
    organism <- "'Homo sapiens"
    reference <- get_test_reference()
    domain <- NA
    outdir <- get_local_tmp_dir()

    # Build argument list.
    args <- c(
      "-e", "devtools::load_all(quiet = TRUE); SPARKI::cli()",
      glue::glue("--std-reports {std_directory}"),
      glue::glue("--mpa-reports {mpa_directory}"),
      glue::glue("--organism {organism}"),
      glue::glue("--reference {reference}"),
      glue::glue("--domain {domain}"),
      glue::glue("--outdir {outdir}")
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

    # Get actual exit status.
    actual_exit_code <- result$status

    # Assert that the exit code is not zero.
    expect_gt(actual_exit_code, 0)
  })

  test_that("CLI fails when the required argument --outdir is not provided", {
    #  Define input arguments for SPARKI.
    std_directory <- get_test_std_report_dir("usual")
    mpa_directory <- get_test_mpa_report_dir("usual")
    organism <- "'Homo sapiens"
    reference <- get_test_reference()
    domain <- "Viruses"
    outdir <- NA

    # Build argument list.
    args <- c(
      "-e", "devtools::load_all(quiet = TRUE); SPARKI::cli()",
      glue::glue("--std-reports {std_directory}"),
      glue::glue("--mpa-reports {mpa_directory}"),
      glue::glue("--organism {organism}"),
      glue::glue("--reference {reference}"),
      glue::glue("--domain {domain}"),
      glue::glue("--outdir {outdir}")
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

    # Get actual exit status.
    actual_exit_code <- result$status

    # Assert that the exit code is not zero.
    expect_gt(actual_exit_code, 0)
  })
})
