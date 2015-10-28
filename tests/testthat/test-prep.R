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

test_that("give a design matrix", {
  design_matrix <- matrix(c(0, 0, 0, 1, 1, 1), ncol = 1)
  result <- sleuth_prep(study_mapping, design_matrix)

  expect_equal(result$design_matrix, design_matrix)
  expect_equal(result$design_matrix, result$full_formula)

  expect_error(sleuth_prep(study_mapping, design_matrix[1:5,]))
})
