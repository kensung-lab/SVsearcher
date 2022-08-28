# SVsearcher
---	
## Introduction
Structural variations (SVs) represent genomic rearrangements such as deletions, insertions, inversions, duplications, and translocations whose sizes are larger than 50bp. A number of long read SV callers have been proposed to call SVs and they perform well. However, the long reads generated by Oxford Nanopore (ONT) have high error rate, which affect the correctness of the long read alignment. Existing long read SV callers do not perform well. We propose a novel method, SVsearcher, to resolve these issues. Compared with existing methods, SVsearcher has highest recall, precision and F1-score.
---	
## Dependence
  1. python3
	2. pysam
	3. Biopython
	4. cigar
	5. numpy
	6. pyvcf
