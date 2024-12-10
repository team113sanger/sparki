test_that("helpers.R get_test_mpa_report_dir() returns a valid path for each case type", {
  for (case_type in c("usual", "empty")) {
    # Assert that the directory exists
    expect_true(fs::dir_exists(get_test_mpa_report_dir(case_type)))
    # Assert that directory contains files with the expected extension .kraken.mpa
    expect_true(length(fs::dir_ls(get_test_mpa_report_dir(case_type), glob = "*.kraken.mpa$")) > 0)
  }
})

test_that("helpers.R get_test_std_report_dir returns a valid path for each case type", {
  for (case_type in c("usual", "empty")) {
    # Assert that the directory exists
    expect_true(fs::dir_exists(get_test_std_report_dir(case_type)))
    # Assert that directory contains files with the expected extension .kraken.mpa
    expect_true(length(fs::dir_ls(get_test_std_report_dir(case_type), glob = "*.kraken$")) > 0)
  }
})

test_that("helpers.R get_test_reference returns valid inspect.txt path", {
  ref_path <- get_test_reference()

  # Check path exists
  expect_true(file.exists(ref_path))

  # Check file content matches expected format
  first_lines <- readLines(ref_path, n = 3)
  expect_match(first_lines[1], "^# Database options:")
  expect_match(first_lines[2], "^# Spaced mask =")
  expect_match(first_lines[3], "^# Toggle mask =")

  # Verify file extension
  expect_equal(tools::file_ext(ref_path), "txt")

  # Verify directory structure
  expect_equal(basename(dirname(ref_path)), "testdata")
})

test_that("helpers.R get_local_tmp_dir creates a directory", {
  # Test directory creation
  dir_path <- get_local_tmp_dir()
  expect_true(dir.exists(dir_path))

  # Test file creation
  test_file <- file.path(dir_path, "test.txt")
  writeLines("test", test_file)
  expect_true(file.exists(test_file))
})

test_that("helpers.R get_local_tmp_dir creates unique directories", {
  dir1 <- get_local_tmp_dir()
  dir2 <- get_local_tmp_dir()

  expect_false(dir1 == dir2)
})

test_that("helpers.R get_test_metadata returns valid inspect.txt path", {
  metadata_path <- get_test_metadata()

  # Check path exists
  expect_true(file.exists(metadata_path))

  # Verify file extension
  expect_equal(tools::file_ext(metadata_path), "txt")

  # Verify directory structure
  expect_equal(basename(dirname(metadata_path)), "testdata")
})
