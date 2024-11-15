test_that("helpers.R get_test_mpa_report_dir() returns a valid path", {
  # Assert that the directory exists
  expect_true(fs::dir_exists(get_test_mpa_report_dir()))
  # Assert that directory contains files with the expected extension .kraken.mpa
  expect_true(length(fs::dir_ls(get_test_mpa_report_dir(), glob = "*.kraken.mpa")) > 0)
})
