# version 0.30.0

This version integrates [p-value aggregation](https://github.com/pachterlab/sleuth/pull/148) as described in [Yi et al.](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1419-z).
The behavior of gene-level differential expression testing now follows this procedure:

1. Isoform-level testing.
2. P-value aggregation at the gene level (using `target_mapping`) by the lancaster method.

Thank you to [Lynn Yi](https://github.com/lynnyi) for implementing p-value aggregation.
Please see [pull request #148](https://github.com/pachterlab/sleuth/pull/148) for details.

The API has also slightly changed. Particularly, for `sleuth_prep`, several options have been moved to optional arguments via `...`. See [pull request #168](https://github.com/pachterlab/sleuth/pull/168) for more information or `?sleuth_prep` in R.

A fair amount of speed up and bug fixes have also been implemented.

- [Patch: bugs in sleuth_results & other miscellaneous fixes](https://github.com/pachterlab/sleuth/pull/163)
- [Fix behavior of sleuth_results when gene_mode is TRUE (and error reporting)](https://github.com/pachterlab/sleuth/pull/160)
- [Shiny and Plot Fixes / Enhancements](https://github.com/pachterlab/sleuth/pull/159)
- [Quick Patch: UseMethod typo](https://github.com/pachterlab/sleuth/pull/157)
- [Update `write_kallisto_hdf5` function and add ability ot subset kallisto object (address #131)](https://github.com/pachterlab/sleuth/pull/150)
- [extend sleuth to model TPMs](https://github.com/pachterlab/sleuth/pull/145)
- [Fixes to various miscellaneous issues (#73, #84, #97, #122, #135, #142)](https://github.com/pachterlab/sleuth/pull/144)
- [Improvements to shiny and plot functions (solving several open issues)](https://github.com/pachterlab/sleuth/pull/143)
- [Possible solution to NAs in sleuth_lrt, addressing #68](https://github.com/pachterlab/sleuth/pull/118)
- [bug fix patches](https://github.com/pachterlab/sleuth/pull/117)
- [address #113 - patch bug where TPM bootstrap summary target_ids are moved](https://github.com/pachterlab/sleuth/pull/116)
- [New tests for ".N" target mappings](https://github.com/pachterlab/sleuth/pull/115)
- [Misc bug fixes + Allow sleuth_prep to process just one sample](https://github.com/pachterlab/sleuth/pull/114)

A major things to [Warren McGee](https://github.com/warrenmcg) for doing the majority of the heavy lifting on all of the bug fixes.


# version 0.29.0

This version has numerous bug fixes and several performance upgrades.
Most notably, memory usage has been decreased greatly by no longer storing the bootstraps in memory.
Additionally, speed has been improved in numerous areas — particularly `sleuth_prep` — by changing several of the computations as well as changing the order of the parallelization (special thanks to [Warren McGee](https://github.com/warrenmcg) for his contributions to this).

Below is an incomplete list of new features:

- The full model no longer has to be specified in `sleuth_prep`.
- A new function `extract_model` allow users to extract the effect sizes for a model in a tidy format similar to [broom](https://cran.r-project.org/web/packages/broom/vignettes/broom.html).
- An arbitrary transformation can be specified/used in `sleuth_prep` (see argument `transformation_function`).

A big thanks to our users for fixing and reporting bugs.
A special thanks to [Warren McGee](https://github.com/warrenmcg) for making several of the performance improvements as well as fixing several bugs.
Below is a partial list of many of the upgrades and the pull requests by the community.

- [Memory overhaul to reduce overall usage](https://github.com/pachterlab/sleuth/pull/63) (@psturmfels)
- [Bugfix to drop unused factors](https://github.com/pachterlab/sleuth/pull/71) (@roryk and @warrenmcg)
- [Reduce memory footprint and improve parallelization](https://github.com/pachterlab/sleuth/pull/94) (@warrenmcg)
- [Add gene annotations when using `sleuth_results`](https://github.com/pachterlab/sleuth/pull/95) (@warrenmcg)
- [Improve sample name handling](https://github.com/pachterlab/sleuth/pull/96) (@warrenmcg)
- [Reconcile memory overhaul and gene aggregation and allow arbitrary transformations](https://github.com/pachterlab/sleuth/pull/99) (@warrenmcg)
- [Do not parallelize when in RStudio](https://github.com/pachterlab/sleuth/pull/108) (@warrenmcg)
- [Remove warning in `sliding_window_grouping`](https://github.com/pachterlab/sleuth/pull/106) (@warrenmcg)
- [Bug fix in `sleuth_live` in gene mode](https://github.com/pachterlab/sleuth/pull/107) (@warrenmcg)
