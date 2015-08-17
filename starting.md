---
layout: page
title: "Getting Started"
description: ""
group: navigation
---
{% include JB/setup %}

<p>To explain how to use <strong>sleuth</strong> we provide an example based on the data in the “Cuffdiff2 paper”:</p>
<ul>
<li><a href="http://www.nature.com/nbt/journal/v31/n1/full/nbt.2450.html">Differential analysis of gene regulation at transcript resolution with RNA-seq</a> by Cole Trapnell, David G Henderickson, Martin Savageau, Loyal Goff, John L Rinn and Lior Pachter, Nature Biotechnology <strong>31</strong>, 46–53 (2013).</li>
</ul>
<p>The human fibroblast RNA-Seq data for the paper is available on GEO at accession <a href="http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE37704">GSE37704</a>. The samples to be analyzed are the six samples LFB_scramble_hiseq_repA, LFB_scramble_hiseq_repB, LFB_scramble_hiseq_repC, LFB_HOXA1KD_hiseq_repA, LFB_HOXA1KD_hiseq_repA, and LFB_HOXA1KD_hiseq_repC. These are three biological replicates in each of two conditions (scramble and HoxA1 knockdown) that will be compared with <strong>sleuth</strong>.</p>
<p>To analyze the data, first download the raw reads, install <strong>kallisto</strong> and then quantify the data with boostraps as described <a href="http://pachterlab.github.io/kallisto/starting.html">here</a>. This step can be skipped for the purposes of the vignette, by downloading the <strong>kallisto</strong> processsed data directly by clicking <a href="http://bio.math.berkeley.edu/sleuth/cuffdiff2/cuffdiff2_data_kallisto_results.zip">here</a>.</p>

<p>
Next, install sleuth in R following the directions <a href="http://pachterlab.github.io/sleuth/download.html">here</a>.
</p>
<p>The first step in a <strong>sleuth</strong> analysis is to specify where the <strong>kallisto</strong> results are stored. Begin by storing the base directory of the results in a variable,</p>
<pre class="sourceCode r"><code class="sourceCode r">base_dir &lt;-<span class="st"> &quot;~/Downloads/cuffdiff2_data_kallisto_results&quot;</span></code></pre>
<p>Next get the list of sample IDs with</p>
<pre class="sourceCode r"><code class="sourceCode r">sample_id &lt;-<span class="st"> </span><span class="kw">dir</span>(<span class="kw">file.path</span>(base_dir,<span class="st">&quot;results&quot;</span>))</code></pre>
<p>The result can be displayed by typing</p>
<pre class="sourceCode r"><code class="sourceCode r">sample_id</code></pre>
<pre><code>## [1] &quot;SRR493366&quot; &quot;SRR493367&quot; &quot;SRR493368&quot; &quot;SRR493369&quot; &quot;SRR493370&quot; &quot;SRR493371&quot;</code></pre>
<p>In the box above, lines beginning with ## show the output of the command (in what follows we include the output that should appear with each command).</p>
<p>A list of paths to the <strong>kallisto</strong> results indexed by the sample IDs is collated with</p>
<pre class="sourceCode r"><code class="sourceCode r">kal_dirs &lt;-<span class="st"> </span><span class="kw">sapply</span>(sample_id, function(id) <span class="kw">file.path</span>(base_dir, <span class="st">&quot;results&quot;</span>, id, <span class="st">&quot;kallisto&quot;</span>))
kal_dirs</code></pre>
<pre><code>##                                                                SRR493366 
## &quot;~/Downloads/cuffdiff2_data_kallisto_results/results/SRR493366/kallisto&quot; 
##                                                                SRR493367 
## &quot;~/Downloads/cuffdiff2_data_kallisto_results/results/SRR493367/kallisto&quot; 
##                                                                SRR493368 
## &quot;~/Downloads/cuffdiff2_data_kallisto_results/results/SRR493368/kallisto&quot; 
##                                                                SRR493369 
## &quot;~/Downloads/cuffdiff2_data_kallisto_results/results/SRR493369/kallisto&quot; 
##                                                                SRR493370 
## &quot;~/Downloads/cuffdiff2_data_kallisto_results/results/SRR493370/kallisto&quot; 
##                                                                SRR493371 
## &quot;~/Downloads/cuffdiff2_data_kallisto_results/results/SRR493371/kallisto&quot;</code></pre>
<p>The next step is to load an auxillary table that describes the experimental design and the relationship between the kallisto directories and the samples:</p>
<pre class="sourceCode r"><code class="sourceCode r">s2c &lt;-<span class="st"> </span><span class="kw">read.table</span>(<span class="kw">file.path</span>(base_dir,<span class="st">&quot;hiseq_info.txt&quot;</span>), <span class="dt">header =</span> <span class="ot">TRUE</span>, <span class="dt">stringsAsFactors=</span><span class="ot">FALSE</span>)
s2c &lt;-<span class="st"> </span>dplyr::<span class="kw">select</span>(s2c, <span class="dt">sample =</span> run_accession, condition)
s2c</code></pre>
<pre><code>##      sample condition
## 1 SRR493366  scramble
## 2 SRR493367  scramble
## 3 SRR493368  scramble
## 4 SRR493369   HOXA1KD
## 5 SRR493370   HOXA1KD
## 6 SRR493371   HOXA1KD</code></pre>
<p>The next stage is to build a “sleuth object”. This requires three commands that (1) load the kallisto processed data into the object (2) estimate parameters for the <strong>sleuth</strong> response error measurement model and (3) perform differential analyis (testing). On a laptop the three steps should take about 2 minutes altogether.</p>
<p>First type</p>
<pre class="sourceCode r"><code class="sourceCode r">so &lt;-<span class="st"> </span><span class="kw">sleuth_prep</span>(kal_dirs, s2c, ~<span class="st"> </span>condition)</code></pre>
<pre><code>## Reading in kallisto results
## ......
## Normalizing 'est_counts'
## 37624 targets passed the filter.
## Normalizing bootstrap samples</code></pre>
<p>then</p>
<pre class="sourceCode r"><code class="sourceCode r">so &lt;-<span class="st"> </span><span class="kw">sleuth_fit</span>(so)</code></pre>
<pre><code>## Summarizing bootstraps
## Fitting measurement error models
## Shrinkage estimation
## Computing variance of betas</code></pre>
<p>and finally</p>
<pre class="sourceCode r"><code class="sourceCode r">so &lt;-<span class="st"> </span><span class="kw">sleuth_test</span>(so, <span class="dt">which_beta =</span> <span class="st">'conditionscramble'</span>)</code></pre>
<p>In general, one can see the possible tests that could be performed using the <code>which_beta</code> parameter in <code>sleuth_test</code> and examining the coefficients:</p>
<pre class="sourceCode r"><code class="sourceCode r"><span class="kw">models</span>(so)</code></pre>
<pre><code>## [  full  ]
## formula:  ~condition 
## coefficients:
##  (Intercept)
##      conditionscramble
## tests:
##  conditionscramble</code></pre>
<p>At this point the sleuth object constructed from the kallisto runs has information about the data, the experimental design, the <strong>kallisto</strong> estimates, the model fit, and the testing. In other words it contains the entire analysis of the data.</p>
<p>The best way to view the results is to generate the Shiny webpage that allows for exploratory data analysis:</p>
<pre><code>sleuth_live(so)</code></pre>
<p>To generate a table of results type</p>
<pre class="sourceCode r"><code class="sourceCode r">results_table &lt;-<span class="st"> </span><span class="kw">sleuth_results</span>(so, <span class="st">'conditionscramble'</span>) </code></pre>



<!-- dynamically load mathjax for compatibility with self-contained -->
<script>
  (function () {
    var script = document.createElement("script");
    script.type = "text/javascript";
    script.src  = "https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML";
    document.getElementsByTagName("head")[0].appendChild(script);
  })();
</script>



