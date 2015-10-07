context("reading")

test_that("get kallisto path", {
  dir_name <- "small_test_data"

  # the standard case
  file_name <- file.path(dir_name, "abundance.h5")
  result <- get_kallisto_path(dir_name)
  expect_equal(result$ext, "h5")
  expect_equal(result$path, file_name)

  # user gives the full path
  result <- get_kallisto_path(file_name)
  expect_equal(result$ext, "h5")
  expect_equal(result$path, file_name)

  # user gives plain text file (full path)
  file_name <- file.path(dir_name, "abundance.tsv")
  result <- get_kallisto_path(file_name)
  expect_equal(result$ext, "tsv")
  expect_equal(result$path, file_name)

  missing_file_name <- file.path(dir_name, "missing.h5")
  expect_error(get_kallisto_path(missing_file_name))

  expect_error(get_kallisto_path("missing_directory"))

})