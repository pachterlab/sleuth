---
layout: page
title: "Download"
group: navigation
---

{% include JB/setup %} 

#### Repository

The __sleuth__ GitHub repository is [here](http://github.com/pachterlab/sleuth).

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

Next load __sleuth__ with

~~~
library("sleuth")
~~~

