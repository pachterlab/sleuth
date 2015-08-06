---
layout: page
title: "FAQ"
group: navigation
---

{% include JB/setup %}

- I'm having trouble with __sleuth__. Can I get help?
  - Yes. If you think you have discovered a bug that needs to be fixed please
    file a report on the GitHub page. If you have a question about installing
    or running the program please ask on the [kallisto-sleuth-users Google user
    group](https://groups.google.com/forum/#!forum/sleuth-sleuth-users).

- Where can I get announcements about new releases?
  - You can get announcements via the [kallisto-sleuth-announcements Google
    group](https://groups.google.com/forum/#!forum/sleuth-sleuth-announcements).
    This is a read-only, low traffic mailing list that only sends an email when
    a major version is released.

- __sleuth__ is spiffy but is it as accurate as other differential expression tools?
  - Yes. __sleuth__ takes advantage of the boostraps of kallisto, thereby effectively leveraging technical replicates in the determination of differential expression. In our tests on both real and simulated data we find that __sleuth__ is _more_ accurate than currently popular differential expression tools such as Cuffdiff2 and DESeq2. 

- Can __sleuth__  be used for quantifying abundances of transcripts or genes from RNA-Seq data?
  - No. __sleuth__ can be used _after_ transcript have been quantified using kallisto.

- Is __sleuth__ usable with both single-end and paired-end reads?
  - Yes.

- Does __sleuth__ require reads to be of the same length?
  - No.

- Can __sleuth__ be used to analyze single-cell RNA-Seq data?
  - Yes. However there are unique challenges in single-cell analysis that are not currently addressed by __sleuth__ (but will be in the near future).

- I've already mapped all my reads and counted the number of alignments to genes. Can I use those mappings with __sleuth__?
  - No. The mappings are not relevant and not needed for __sleuth__.


- Are you distributing for use with __sleuth__?
  - No. __sleuth__ relies on the transcriptome used to quantify the RNA-Seq data being analyzed.

- My RNA-Seq was prepared with a stranded library. Is there a special option I need to use with sleuth?
  - No.

- Is there a reason you picked the name sleuth for your program?
  - Yes.
