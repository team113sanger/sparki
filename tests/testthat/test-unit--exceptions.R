test_that("check_directory() returns an error if a directory does not exist", {
    logger::log_threshold(logger::OFF)

    # Define path to a directory that does not exist, and ensure it really
    # does not exist.
    dirpath <- "directory_that_does_not_exist/"
    expect_false(dir.exists(dirpath))

    # Carry out test.
    expect_error(SPARKI:::check_directory(dirpath))
})

test_that("check_directory() handles empty directories appropriately depending on the expectation", {
    logger::log_threshold(logger::OFF)

    # Define path to an empty directory.
    dirpath <- get_test_empty_dir()

    # Carry out tests.
    expect_error(SPARKI:::check_directory(dirpath, expectation = "not_empty"))
    expect_no_error(SPARKI:::check_directory(dirpath, expectation = "empty"))
})

test_that("check_directory() handles non-empty directories appropriately depending on the expectation", {
    logger::log_threshold(logger::OFF)

    # Define path to a non-empty directory.
    dirpath <- get_test_not_empty_dir()

    # Carry out tests.
    expect_error(SPARKI:::check_directory(dirpath, expectation = "empty"))
    expect_no_error(SPARKI:::check_directory(dirpath, expectation = "not_empty"))
})

test_that("check_directory() returns a valid path after adding a slash symbol to the end of the path", {
    logger::log_threshold(logger::OFF)

    # Define path to a non-empty directory.
    dirpath <- get_test_not_empty_dir()

    # Carry out tests.
    post_check_dirpath <- SPARKI:::check_directory(dirpath, expectation = "not_empty")
    expect_true(substr(post_check_dirpath, nchar(post_check_dirpath), nchar(post_check_dirpath)) == "/")
    expect_true(dir.exists(post_check_dirpath))
})
