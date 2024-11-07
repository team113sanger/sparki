#################################
## Tests for load_MPAreports() ##
#################################

test_that("load_MPAreports() returns a dataframe with the appropriate dimensions", {
  #  Define the expected number of columns.
  mpa_reports_dir <- get_test_mpa_report_dir()
  expected_number_of_columns <- 12

  # Load the MPA-style reports and obtain the actual numbers of rows and columns.
  out <- SPARKI::load_MPAreports(mpa_reports_dir)
  actual_number_of_rows <- nrow(out)
  actual_number_of_columns <- ncol(out)

  #  Carry out tests.
  expect_false(actual_number_of_rows == 0)
  expect_false(actual_number_of_columns == 0)
  expect_false(is.null(dim(out)))
  expect_equal(actual_number_of_columns, expected_number_of_columns)
})

test_that("the columns in the dataframe returned by load_MPAreports() correspond to the appropriate classes", {
  #  Define the expected number of columns.
  mpa_reports_dir <- get_test_mpa_report_dir()
  expected_number_of_columns <- 12

  # Load the MPA-style reports and obtain the actual numbers of rows and columns.
  out <- SPARKI::load_MPAreports(mpa_reports_dir)
  actual_number_of_rows <- nrow(out)
  actual_number_of_columns <- ncol(out)

  #  Carry out tests.
  expect_false(actual_number_of_rows == 0)
  expect_false(actual_number_of_columns == 0)
  expect_false(is.null(dim(out)))
  expect_equal(actual_number_of_columns, expected_number_of_columns)
})




#################################
## Tests for load_STDreports() ##
#################################

test_that("load_STDreports() returns a dataframe with the appropriate dimensions", {
  #  Define the expected number of columns.
  mpa_reports_dir <- get_test_mpa_report_dir()
  expected_number_of_columns <- 9

  # Load the MPA-style reports and obtain the actual numbers of rows and columns.
  out <- SPARKI::load_STDreports(mpa_reports_dir)
  actual_number_of_rows <- nrow(out)
  actual_number_of_columns <- ncol(out)

  #  Carry out tests.
  expect_false(actual_number_of_rows == 0)
  expect_false(actual_number_of_columns == 0)
  expect_false(is.null(dim(out)))
  expect_equal(actual_number_of_columns, expected_number_of_columns)
})
