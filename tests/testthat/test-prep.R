context("preparing a data set")


test_that("ensure that samples and directories are included", {
  expect_error(sleuth_prep(study_formula, study_mapping))
  expect_error(sleuth_prep(3, study_mapping))
  expect_error(sleuth_prep(study_mapping, 3))
})

test_that("normalization factors", {
  expect_error(sleuth_prep(study_mapping, ~condition, norm_fun_counts = 3))
  all_ones <- function(x) {
    p <- ncol(x)
    sf <- rep.int(1, p)
    names(sf) <- colnames(x)

    sf
  }
  result <- sleuth_prep(study_mapping, ~condition, norm_fun_counts = all_ones)
  temp_result <- rep.int(1, nrow(study_mapping))
  expect_equivalent(result$est_counts_sf, temp_result)
})
