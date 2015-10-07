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

test_that("both read types", {
  dir_name <- "small_test_data"

  h5_file_name <- file.path(dir_name, "abundance.h5")
  kal_h5 <- read_kallisto_h5(h5_file_name, read_bootstrap = FALSE)
  kal_h5$abundance <- dplyr::arrange(kal_h5$abundance, target_id)

  tsv_file_name <- file.path(dir_name, "abundance.tsv")
  kal_tsv <- read_kallisto_tsv(tsv_file_name)

  expect_equal(kal_h5$abundance$target_id, kal_tsv$abundance$target_id)
  expect_equal(ncol(kal_h5$abundance), ncol(kal_tsv$abundance))
})