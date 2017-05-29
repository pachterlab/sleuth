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
