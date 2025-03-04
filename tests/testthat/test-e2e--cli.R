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
        "--verbosity off", # Not a required argument but we do not want log messages to be printed out during the test.
        sep = " "
    )

    exit_code <- system(
        command = cli_cmd,
        ignore.stdout = TRUE, # Ignoring this as we do not want log messages to be printed out during the test.
        ignore.stderr = FALSE,
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
        "--verbosity off", # Not a required argument but we do not want log messages to be printed out during the test.
        sep = " "
    )

    exit_code <- system(
        command = cli_cmd,
        ignore.stdout = TRUE, # Ignoring this as we do not want log messages to be printed out during the test.
        ignore.stderr = FALSE,
        intern = FALSE
    )

    expect_equal(exit_code, 0)
})

test_that("cli() fails when the user-defined species' delimiter is neither a white space or an underscore", {
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
        "--verbosity off", # Not a required argument but we do not want log messages to be printed out during the test.
        sep = " "
    )

    exit_code <- system(
        command = cli_cmd,
        ignore.stdout = TRUE, # Ignoring this as we do not want log messages to be printed out during the test.
        ignore.stderr = TRUE, # In this case we know this will fail so we will ignore the standard error as well.
        intern = FALSE
    )

    expect_gt(exit_code, 0)
})
