context("preparing a data set")


test_that("ensure that samples and directories are included", {
  expect_error(sleuth_prep(study_formula, study_mapping))
  expect_error(sleuth_prep(3, study_mapping))
  expect_error(sleuth_prep(study_mapping, 3))


})
