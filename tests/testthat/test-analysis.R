context("running full sleuth pipeline")

test_that("transcript level analysis", {
  # transcript level results
  transcript_result <- sleuth_prep(study_mapping, study_formula)
  transcript_result <- sleuth_fit(transcript_result)

  expect_equal(transcript_result$fits$full$models, trans_test_data$fits$full$models)
  expect_equal(transcript_result$fits$full$beta_covars, trans_test_data$fits$full$beta_covars)

  transcript_result <- sleuth_wt(transcript_result, 'conditionwildtype')
  test_table <- sleuth_results(transcript_result, 'conditionwildtype')
  comparison <- sleuth_results(trans_test_data, 'conditionwildtype')

  expect_equal(test_table, comparison)
})

test_that("gene level analysis", {
  # gene level results
  gene_result <- sleuth_prep(study_mapping, study_formula,
                        target_mapping = target_mapping,
                        aggregation_column = "gene_name")
  gene_result <- sleuth_fit(gene_result)

  expect_equal(gene_result$fits$full$models, gene_test_data$fits$full$models)
  expect_equal(gene_result$fits$full$beta_covars, gene_test_data$fits$full$beta_covars)

  gene_result <- sleuth_wt(gene_result, 'conditionwildtype')
  test_table <- sleuth_results(gene_result, 'conditionwildtype')
  comparison <- sleuth_results(gene_test_data, 'conditionwildtype')

  expect_equal(test_table, comparison)
})
