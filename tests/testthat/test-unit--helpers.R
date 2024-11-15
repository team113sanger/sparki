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
