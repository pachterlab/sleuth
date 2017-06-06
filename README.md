# sleuth

__sleuth__ is a program for differential analysis of RNA-Seq data. It makes use of quantification uncertainty estimates obtained via [kallisto](https://github.com/pimentel/kallisto) for accurate differential analysis of isoforms or genes, allows testing in the context of experiments with complex designs, and supports interactive exploratory data analysis via __sleuth live__. The sleuth methods are described in

H Pimentel,	NL Bray,	S Puente,	P Melsted	and Lior Pachter, [Differential analysis of RNA-seq incorporating quantification uncertainty](http://www.nature.com/nmeth/journal/vaop/ncurrent/full/nmeth.4324.html), Nature Methods (2017), advanced access.

Scripts reproducing all the results of the paper are available [here](https://github.com/pachterlab/sleuth_paper_analysis).

# Installation

The easiest way to install is using the `devtools` package through Bioconductor.

```{r}
source("http://bioconductor.org/biocLite.R")
biocLite("devtools")    # only if devtools not yet installed
biocLite("pachterlab/sleuth")
```

These commands will install `sleuth` along with all of its dependencies. You
can then load `sleuth` like any other R package:

```{r}
library('sleuth')
```

## Installation via `conda`

If you have [`conda`](http://conda.pydata.org/docs/), a cross-platform package manager installed, you can install `sleuth` via the [`bioconda`](https://bioconda.github.io/) channel.

```
conda install --channel bioconda r-sleuth
```

# Documentation

We recommend starting with the vignette:

```{r}
vignette('intro', package = 'sleuth')
```

Detailed documentation can be retrieved within R using the `help()` command:

```{r}
help(package = 'sleuth')
```

Specific function documentation can also be accessed using `?` as you would for
any other function in R:

```{r}
?sleuth_prep
```

# Conventions

- All sleuth "core" functionality is prefixed by `sleuth_` (e.g.
`sleuth_prep()`).
- All sleuth plots are prefixed with `plot_` (e.g. `plot_ma()`)


# Further help

Please visit the [sleuth website](https://pachterlab.github.io/sleuth) for ways to get help.
We have several new [walk-throughs](https://pachterlab.github.io/sleuth/walkthroughs) at the main sleuth website.
In particular, you might find the [kallisto-sleuth users](https://groups.google.com/forum/#!forum/kallisto-sleuth-users) Google group helpful.

Please post bug reports [on GitHub](https://github.com/pachterlab/sleuth/issues).

# Copyright

Copyright (C) 2017 Harold Pimentel, Nicolas Bray, Pall Melsted, Lior Pachter
