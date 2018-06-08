This data set is a subset of the data from the Ellahi et. al. paper.
The actual analysis can be found at: https://github.com/pachterlab/bears_analyses

The target mappings were collected using biomaRt, Ensembl v82 (Sep. 2015), and
extracting out the attributes "ensembl_transcript_id" and "external_gene_name."
Some Yeast transcripts do not have an external gene name, so target_mappings.txt
uses the transcript ID as the gene name, and target_mappings_incomplete.txt has
those entries with an "NA" (and thus excluded from gene-level tests).

The RDS files were computed using sleuth version 0.30.0. Here was the code to
generate those files, assuming R is in this directory:

```
library(sleuth)

data_path <- getwd()
sample_ids <- grep('^SR', dir(data_path), value = TRUE)
sample_ids <- rev(sample_ids)

target_mapping <- read.table(file.path(data_path, 'target_mappings.txt'), header = TRUE,
  stringsAsFactors = FALSE, sep="\t", quote="")

study_mapping <- read.table(file.path(data_path, 'study_design.txt'), header = TRUE,
  stringsAsFactors = FALSE)
study_mapping <- dplyr::select(study_mapping, sample = run, condition)
result_paths <- file.path(data_path, sample_ids, 'kallisto')
study_mapping <- dplyr::mutate(study_mapping, path = result_paths)

study_formula <- ~condition

trans_so <- sleuth_prep(study_mapping, study_formula)
trans_so <- sleuth_fit(trans_so)
trans_so <- sleuth_wt(trans_so, 'conditionwildtype')
sleuth_save(trans_so, 'ellahi_transcript.rds')

gene_so <- sleuth_prep(study_mapping, study_formula, target_mapping = target_mapping,
                       aggregation_column = 'gene_name', gene_mode = TRUE)
gene_so <- sleuth_fit(gene_so)
gene_so <- sleuth_wt(gene_so, 'conditionwildtype')
sleuth_save(gene_so, 'ellahi_gene.rds')
```

The data set is small enough that you can use it as a test while you develop.
