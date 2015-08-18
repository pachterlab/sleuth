# sleuth: inspect your RNA-Seq

Investigate transcript abundance from
[kallisto](https://github.com/pimentel/kallisto) and differential expression
analysis from RNA-Seq data.

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

# Documentation

We recommend going through the vignette, first though:

```{r}
vignette('intro', package = 'sleuth')
```

Detailed documentation can be accessed inside of R using the `help()` command:

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
users](https://groups.google.com/forum/#!forum/sleuth-sleuth-users) Google
group helpful.

# Copyright

Copyright (C) 2015  Harold Pimentel, Nicolas Bray, Pall Melsted, Lior Pachter
