test_that("It is possible to retrieve the package version through cli()", {
    logger::log_threshold(logger::OFF)

    rscript <- Sys.which("Rscript")

    # Build CLI command.
    cli_cmd <- paste(
        rscript, "-e", glue::single_quote("devtools::load_all(quiet = TRUE); SPARKI::cli()"),
        "--version", sep = " "
    )

    # Launch the CLI and check that the exit code is zero (i.e. success!).
    exit_code <- system(
        command = cli_cmd,
        ignore.stdout = TRUE, # Ignoring this as we do not want log messages to be printed out during the test.
        ignore.stderr = TRUE, # Ignoring this as we do not want log messages to be printed out during the test.
        intern = FALSE
    )
    expect_equal(exit_code, 0)

    # Launch the CLI again, now capturing the output.
    out <- system(
        command = cli_cmd,
        ignore.stdout = FALSE,
        ignore.stderr = FALSE,
        intern = TRUE
    )
    expect_true(grepl(packageVersion("SPARKI"), out))
})

test_that("It is possible to get help through cli()", {
    logger::log_threshold(logger::OFF)

    rscript <- Sys.which("Rscript")

    # Build CLI command.
    cli_cmd <- paste(
        rscript, "-e", glue::single_quote("devtools::load_all(quiet = TRUE); SPARKI::cli()"),
        "--help", sep = " "
    )

    # Launch the CLI and check that the exit code is zero (i.e. success!).
    exit_code <- system(
        command = cli_cmd,
        ignore.stdout = TRUE, # Ignoring this as we do not want log messages to be printed out during the test.
        ignore.stderr = TRUE, # Ignoring this as we do not want log messages to be printed out during the test.
        intern = FALSE
    )
    expect_equal(exit_code, 0)
})


test_that("cli() works when only required arguments are specified", {
    logger::log_threshold(logger::OFF)

    rscript <- Sys.which("Rscript")

    # Build CLI command.
    cli_cmd <- paste(
        rscript, "-e", glue::single_quote("devtools::load_all(quiet = TRUE); SPARKI::cli()"),
        glue::glue("--std-reports", {get_test_std_report_dir("valid")}, .sep = " "),
        glue::glue("--mpa-reports", {get_test_mpa_report_dir("valid")}, .sep = " "),
        glue::glue("--reference", {get_test_reference()}, .sep = " "),
        glue::glue("--outdir", {get_local_tmp_dir()}, .sep = " "),
        "--organism 'Homo sapiens'",
        "--domain Viruses",
        "--verbosity trace", # We want to test at the trace level to ensure logging will be produced with no errors;
        sep = " "            # the log messages will exist but they will not be printed out (as defined below).
    )

    # Launch the CLI and check that the exit code is zero (i.e. success!).
    exit_code <- system(
        command = cli_cmd,
        ignore.stdout = TRUE, # Ignoring this as we do not want log messages to be printed out during the test.
        ignore.stderr = TRUE, # Ignoring this as we do not want log messages to be printed out during the test.
        intern = FALSE
    )
    expect_equal(exit_code, 0)
})

test_that("cli() works when the report contains white space-delimited species names but the user-defined species is underscore-delimited", {
    logger::log_threshold(logger::OFF)

    rscript <- Sys.which("Rscript")

    # Build CLI command.
    cli_cmd <- paste(
        rscript, "-e", glue::single_quote("devtools::load_all(quiet = TRUE); SPARKI::cli()"),
        glue::glue("--std-reports", {get_test_std_report_dir("valid")}, .sep = " "),
        glue::glue("--mpa-reports", {get_test_mpa_report_dir("valid")}, .sep = " "),
        glue::glue("--reference", {get_test_reference()}, .sep = " "),
        glue::glue("--outdir", {get_local_tmp_dir()}, .sep = " "),
        "--organism 'Homo_sapiens'", # This species is present in the report with a white space as the delimiter!
        "--domain Viruses",
        "--verbosity trace", # We want to test at the trace level to ensure logging will be produced with no errors;
        sep = " "            # the log messages will exist but they will not be printed out (as defined below).
    )

    # Launch the CLI and check that the exit code is zero (i.e. success!).
    exit_code <- system(
        command = cli_cmd,
        ignore.stdout = TRUE, # Ignoring this as we do not want log messages to be printed out during the test.
        ignore.stderr = TRUE, # Ignoring this as we do not want log messages to be printed out during the test.
        intern = FALSE
    )
    expect_equal(exit_code, 0)
})

test_that("cli() fails when the user-defined species's delimiter is neither a white space or an underscore", {
    logger::log_threshold(logger::OFF)

    rscript <- Sys.which("Rscript")

    # Build CLI command.
    cli_cmd <- paste(
        rscript, "-e", glue::single_quote("devtools::load_all(quiet = TRUE); SPARKI::cli()"),
        glue::glue("--std-reports", {get_test_std_report_dir("valid")}, .sep = " "),
        glue::glue("--mpa-reports", {get_test_mpa_report_dir("valid")}, .sep = " "),
        glue::glue("--reference", {get_test_reference()}, .sep = " "),
        glue::glue("--outdir", {get_local_tmp_dir()}, .sep = " "),
        "--organism 'Homo-sapiens'", # This species is present in the report with a white space as the delimiter!
        "--domain Viruses",
        "--verbosity trace", # We want to test at the trace level to ensure logging will be produced with no errors;
        sep = " "            # the log messages will exist but they will not be printed out (as defined below).
    )

    # Launch the CLI and check that the exit code is not zero (i.e. failure).
    exit_code <- system(
        command = cli_cmd,
        ignore.stdout = TRUE, # Ignoring this as we do not want log messages to be printed out during the test.
        ignore.stderr = TRUE, # Ignoring this as we do not want log messages to be printed out during the test.
        intern = FALSE
    )
    expect_gt(exit_code, 0)
})
