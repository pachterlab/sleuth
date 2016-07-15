---
layout: page
title: "Download"
group: navigation
---

{% include JB/setup %} 

#### Installation

To install __sleuth__ start [R](https://www.r-project.org) and first install `rhdf5` by typing: 

~~~
source("http://bioconductor.org/biocLite.R")
biocLite("rhdf5")
~~~

Then install devtools by typing

~~~
install.packages("devtools")
~~~

and install __sleuth__ by typing

~~~
devtools::install_github("pachterlab/sleuth")
~~~

If you have [__conda__](http://conda.pydata.org/docs/), a cross-platform package manager installed, you can install __sleuth__ via the [__bioconda__](https://bioconda.github.io/) channel.

~~~
conda install --channel bioconda r-sleuth
~~~

Next load __sleuth__ with

~~~
library("sleuth")
~~~

#### Repository

The __sleuth__ GitHub repository is [here](http://github.com/pachterlab/sleuth).

