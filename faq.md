---
layout: page
title: "FAQ"
group: navigation
---

{% include JB/setup %}

- I'm having trouble with __sleuth__. Can I get help?
  - Yes. If you think you have discovered a bug that needs to be fixed please
    file a report on the GitHub page. If you have a question about installing
    or running the program, please ask on the [kallisto-and-applications Google user group](https://groups.google.com/forum/#!forum/kallisto-and-applications).

- __sleuth__ is spiffy but is it as accurate as other differential expression tools?
  - Yes. __sleuth__ takes advantage of the boostraps of kallisto, thereby effectively leveraging technical replicates in the determination of differential expression. In the tests performed on both real and simulated data in the [sleuth paper](http://biorxiv.org/content/early/2016/06/10/058164), we find that __sleuth__ is _more_ accurate than currently popular differential expression tools such as Cuffdiff2, DESeq2 and edgeR.

- Can __sleuth__  be used for quantifying abundances of transcripts or genes from RNA-Seq data?
  - No. __sleuth__ is used _after_ transcripts have been quantified using [__kallisto__](http://pachterlab.github.io/kallisto/).

- Is __sleuth__ usable with both single-end and paired-end reads?
  - Yes.

- Does __sleuth__ require reads to be of the same length?
  - No.

- Can __sleuth__ be used to analyze single-cell RNA-Seq data?
  - Yes. However there are unique challenges in single-cell analysis that are not currently addressed by __sleuth__ (but will be in the near future).

- I've already mapped all my reads and counted the number of alignments to genes. Can I use those mappings with __sleuth__?
  - No. The mappings are not relevant and therefore cannot be used with __sleuth__.


- Are you distributing transcriptomes for use with __sleuth__?
  - No. This is unnecessary because __sleuth__ obtains transcript names from the kallisto quantification output. 

- My RNA-Seq was prepared with a stranded library. Is there a special option I need to use with sleuth?
  - No.

- Is there a reason you picked the name sleuth for your program?
  - Yes.
