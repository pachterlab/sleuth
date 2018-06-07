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

  expect_error(result <- sleuth_prep(study_mapping, design_matrix))
  expect_error(sleuth_prep(study_mapping, design_matrix[1:5,]))

  colnames(design_matrix) <- c('Intercept')
  rownames(design_matrix) <- study_mapping$sample

  result <- sleuth_prep(study_mapping, design_matrix)

  expect_equal(result$design_matrix, design_matrix)
  expect_equal(result$design_matrix, result$full_formula)

  expect_equal(result$bs_summary, trans_test_data$bs_summary)
  expect_equal(result$obs_norm, trans_test_data$obs_norm)
  expect_equal(result$obs_norm_filt, trans_test_data$obs_norm_filt)
})

test_that("gene level", {
  expect_error(result <- sleuth_prep(study_mapping, study_formula,
                                     aggregation_column = "gene_name"))
  result <- sleuth_prep(study_mapping, study_formula,
                        target_mapping = target_mapping,
                        aggregation_column = "gene_name", gene_mode = TRUE)

  expect_equal(result$bs_summary, gene_test_data$bs_summary)
  expect_equal(result$obs_norm, gene_test_data$obs_norm)
  expect_equal(result$obs_norm_filt, gene_test_data$obs_norm_filt)

  expect_warning(result_incomplete <- sleuth_prep(study_mapping, study_formula,
                        target_mapping = incomplete_mapping,
                        aggregation_column = "gene_name", gene_mode = TRUE))
})

test_that(".N target mappings", {
  expect_warning(result.N <- sleuth_prep(small_study_map,
                          target_mapping = small_target_mapping))
  expect_warning(result.N <- sleuth_prep(small_study_map,
                          target_mapping = small_target_mapping,
                          aggregation_column = "gene_name", gene_mode = TRUE))
})

test_that("duplicated IDS", {
  test_mapping <- small_target_mapping
  test_mapping$target_id[2] <- test_mapping$target_id[1]
  expect_warning(result_dup <- sleuth_prep(small_study_map,
                            target_mapping = test_mapping))
  expect_error(result_dup <- sleuth_prep(small_study_map,
                          target_mapping = test_mapping,
                          aggregation_column = "gene_name", gene_mode = TRUE))
})

test_that("directly changing sleuth object", {
  test_data <- sleuth_prep(study_mapping, study_formula, read_bootstrap_tpm = TRUE)
  expect_error(transform_fun_counts(test_data) <- function(x) 1)
  expect_error(transform_fun_tpm(test_data) <- function(x) 1)

  test_data <- sleuth_fit(test_data)
  test_data <- sleuth_fit(test_data, fit_name = "tpm", which_var = "obs_tpm")
  expect_error(transform_fun_counts(test_data) <- 2)
  expect_error(transform_fun_tpm(test_data) <- 2)
  expect_warning(transform_fun_counts(test_data) <- function(x) 3)
  expect_error(transform_fun_counts(test_data) <- function(x) 4)
  expect_warning(transform_fun_tpm(test_data) <- function(x) 3)
  expect_error(transform_fun_tpm(test_data) <- function(x) 4)

  expect_error(p <- plot_mean_var(test_data))
  expect_error(p <- plot_vars(test_data))
  expect_error(test_data <- sleuth_wt(test_data))
  expect_error(test_data <- sleuth_lrt(test, "tpm", "full"))
  expect_error(test_data <- sleuth_fit(test_data, ~condition+0, fit_name = "new"))

  expect_error(test_data$norm_fun_counts <- function(x) 1)
  expect_error(test_data$norm_fun_tpm <- function(x) 1)
})
