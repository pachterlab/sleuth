## ----eval=FALSE----------------------------------------------------------
#  source("http://bioconductor.org/biocLite.R")
#  biocLite("rhdf5")

## ----eval=FALSE----------------------------------------------------------
#  install.packages("devtools")

## ----eval=FALSE----------------------------------------------------------
#  devtools::install_github("pachterlab/sleuth")

## ------------------------------------------------------------------------
library("sleuth")

## ------------------------------------------------------------------------
base_dir <- "~/Downloads/cuffdiff2_data_kallisto_results"

## ------------------------------------------------------------------------
sample_id <- dir(file.path(base_dir,"results"))

## ------------------------------------------------------------------------
sample_id

## ------------------------------------------------------------------------
kal_dirs <- sapply(sample_id, function(id) file.path(base_dir, "results", id, "kallisto"))
kal_dirs

## ------------------------------------------------------------------------
s2c <- read.table(file.path(base_dir,"hiseq_info.txt"), header = TRUE, stringsAsFactors=FALSE)
s2c <- dplyr::select(s2c, sample = run_accession, condition)
s2c

## ----eval=TRUE-----------------------------------------------------------
so <- sleuth_prep(kal_dirs, s2c, ~ condition)

## ----eval=TRUE-----------------------------------------------------------
so <- sleuth_fit(so)

## ----eval=TRUE-----------------------------------------------------------
so <- sleuth_test(so, which_beta = 'conditionscramble')

## ----eval=TRUE-----------------------------------------------------------
models(so)

## ------------------------------------------------------------------------
results_table <- sleuth_results(so, 'conditionscramble') 

