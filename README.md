# sleuth: inspect your RNA-Seq

Investigate RNA-Seq transcript abundance from
[kallisto](https://github.com/pimentel/kallisto) and perform differential expression
analysis.

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

Please visit the [sleuth website](https://pachterlab.github.io/sleuth) for
ways to get help. In particular, you might find the [kallisto-sleuth
users](https://groups.google.com/forum/#!forum/kallisto-sleuth-users) Google
group helpful.

Please post bug reports [on GitHub](https://github.com/pachterlab/sleuth/issues).

# Copyright

Copyright (C) 2015  Harold Pimentel, Nicolas Bray, Pall Melsted, Lior Pachter
