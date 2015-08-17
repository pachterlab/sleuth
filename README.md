# sleuth: inspect your RNA-Seq

Investigate transcript abundance from
[kallisto](https://github.com/pimentel/kallisto) and differential expression
analysis from RNA-Seq data.

# Dependencies

All depencies except `rhd5` are available on CRAN so they should get pulled in
when you install `sleuth`. Here is how to `rhdf5`:

```{r}
source("http://bioconductor.org/biocLite.R")
biocLite("rhdf5")
```

# Installation

The easiest way to install is using the `devtools` package. If you don't have
the devtools package, you can get it from CRAN:

```
install.packages('devtools')
```

Then you can install sleuth:

```
devtools::install_github('pachterlab/github')
```

# Copyright

Copyright (C) 2015  Harold Pimentel, Nicolas Bray, Pall Melsted, Lior Pachter
