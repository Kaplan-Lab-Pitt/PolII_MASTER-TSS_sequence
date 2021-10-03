# Among 45 samples, 27 samples have 2 batches of TSS-seq data.

# Preliminary analysis for 2 batches were processed seperately.
# After investigating reproducibility, data from 2 batches were merged after DNA-seq & TSS-seq merging step.

# Analysis for both batches are the same. Scipts here are the version for batch-1.
# Note two batches have:
     # different file names and number of fastq files;
     # different file names of intermediate files ("b1" vs "b2");
     # different read lengthes, which therefore affect the value of "max_len" parameter in <7-MkMockSam-UmiBarCXed-TssOnUmiBar.py>. Details are documented in scripts.