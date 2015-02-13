# sleuth: inspect your RNA-Seq

Investigate transcript abundance from
[kallisto](https://github.com/pimentel/kallisto) and differential expression
analysis from RNA-Seq data.


# Installation

The easiest way to install is using the `devtools` package. The package is
changing quite a bit in this stage, so I suggest you clone the repo, then
install using "devtools":

```
install("~/path/to/sleuth")
```

# Development workflow

Clone this repository and make sure you have `devtools` and `roxygen2`
installed in `R`. Then move to the repo:

```{r}
install.packages(c("devtools", "roxygen2"))
library("devtools")
library("roxygen2")
setwd("~/path/to/sleuth")
```

## The build process

Typically as you write new functions, you'll need to update the NAMESPACE. This
can be done with `roxygen2` documentation in the file (e.g. the `@export` tag).
Then you can run:

```{r}
document()
install()
```

`document()` updates the NAMESPACE file, `install()` installed the package.
