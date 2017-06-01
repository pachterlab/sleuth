# prepare a data set to be used below

data_path <- 'ellahi'

sample_ids <- grep('^SR', dir(data_path), value = TRUE)
sample_ids <- rev(sample_ids)

target_mapping <- read.table(file.path(data_path, 'target_mappings.txt'), header = TRUE,
  stringsAsFactors = FALSE, sep="\t", quote="")
incomplete_mapping <- read.table(file.path(data_path, 'target_mappings_incomplete.txt'), header = TRUE,
  stringsAsFactors = FALSE, sep="\t", quote="")

small_study_map <- data.frame(sample = "small_sample", condition = "test",
                                path = "small_test_data/kallisto.N",
                                stringsAsFactors = F)
small_target_mapping <- read.table('small_test_data/target_mapping.txt', header = TRUE,
  stringsAsFactors = FALSE, sep="\t", quote="")

study_mapping <- read.table(file.path(data_path, 'study_design.txt'), header = TRUE,
  stringsAsFactors = FALSE)
study_mapping <- dplyr::select(study_mapping, sample = run, condition)

stopifnot(sample_ids == study_mapping$sample)

result_paths <- file.path(data_path, sample_ids, 'kallisto')
study_mapping <- dplyr::mutate(study_mapping, path = result_paths)

study_formula <- ~condition

trans_test_data <- sleuth_load(file.path(data_path, 'ellahi_transcript.rda'))
gene_test_data <- sleuth_load(file.path(data_path, 'ellahi_gene.rda'))
