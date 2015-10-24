# prepare a data set to be used below

data_path <- 'ellahi'

sample_ids <- grep('^SR', dir(data_path), value = TRUE)
sample_ids <- rev(sample_ids)

study_mapping <- read.table(file.path(data_path, 'study_design.txt'), header = TRUE,
  stringsAsFactors = FALSE)
study_mapping <- dplyr::select(study_mapping, sample = run, condition)

stopifnot(sample_ids == study_mapping$sample)

result_paths <- file.path(data_path, sample_ids, 'kallisto')
study_mapping <- dplyr::mutate(study_mapping, path = result_paths)

study_formula <- ~condition
