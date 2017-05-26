This data set is a subset of the data from the Ellahi et. al. paper.
The actual analysis can be found at: https://github.com/pachterlab/bears_analyses

The target mappings were collected using biomaRt, Ensembl v82 (Sep. 2015), and
extracting out the attributes "ensembl_transcript_id" and "external_gene_name."
Some Yeast transcripts do not have an external gene name, so target_mappings.txt
uses the transcript ID as the gene name, and target_mappings_incomplete.txt has
those entries with an "NA" (and thus excluded from gene-level tests).

The rda files were computed using sleuth version 0.28.1.
(Warren's forked repository at github.com/warrenmcg/sleuth, commit 325c267).
Both were run using read_bootstrap_tpm and extra_bootstrap_summary, and
the gene-level results used target_mappings.txt (full mappings).

The data set is small enough that you can use it as a test while you develop.
