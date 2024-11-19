
test_that("run_sparki - required options", {
    # Setup
    std_reports_path <- get_test_std_report_dir("usual")
    mpa_reports_path <- get_test_mpa_report_dir("usual")
    organism <- "Homo sapiens"
    reference_path <- get_test_reference()
    metadata_path <- NA
    metadata_sample_col <- NA
    metadata_columns <- NA
    outdir_path <- get_local_tmp_dir()
    prefix <- NA
    verbose <- TRUE
    include_eukaryotes <- FALSE
    include_sample_names <- FALSE
    domain <- "Viruses,Bacteria"
    remove <- NA

    # Check the outdir is empty
    initial_outdir_size <- length(list.files(outdir_path))
    expect_equal(initial_outdir_size, 0)


    # Build the command and get the expected version
    SPARKI::run_sparki(
        std_reports_path,
        mpa_reports_path,
        organism,
        reference_path,
        metadata_path,
        metadata_sample_col,
        metadata_columns,
        outdir_path,
        prefix,
        verbose,
        include_eukaryotes,
        include_sample_names,
        domain,
        remove
    )

    # Assert that the outdir size is non-zero
    final_outdir_size <- length(list.files(outdir_path))
    expect_gt(final_outdir_size, 0)
    }
)
