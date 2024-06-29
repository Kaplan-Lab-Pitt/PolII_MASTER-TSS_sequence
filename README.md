# Pol II MASTER data analysis

This repository includes scripts used for analyzing "Pol II MASTER" (Pol II MAssively Systematic Transcript End Readout) libraies studying how sequence and Pol II activity determine transcription start site selection in budding yeast *Saccharomyces cerevisiae*.
 - Yunye Zhu, Irina O. Vvedenskaya, Sing-Hoi Sze, Bryce E. Nickels, and Craig D. Kaplan. ["Quantitative analysis of transcription start site selection reveals control by DNA sequence, RNA polymerase II activity and NTP levels."](https://rdcu.be/du9F1) *Nature Structural & Molecular Biology* 31.1 (2024): 190-202.

Three major parts: (1) Pol II MASTER libraries analysis, including DNA-seq analysis, TSS-seq analysis, TSS sequence preferece analysis; (2) genomic TSS-seq analysis; (3) modeling analysis.

## Pol II MASTER libraries analysis
The purpose for DNA-seq analysis is to assign each 9 nt randomized TSS sequence to a corresponding 20N barcode. Scripts in folder [1-DNAseq] performed following analyses:
1. paired-end reads merging using PEAR;
2. TSS sequence variant and corresponding barcode extraction;
3. TSS and barcode error correction using UMI-tools;
4. TSS-barcode linkages pool generation.
The final output is tab delimited files containing TSS-Barcode linkages and coresponding DNA-seq read counts in each *E.coli* or yeast samples.

The purpose for TSS-seq analysis is to link RNA products to barcodes, therefore assign TSS usage to corresponding DNA templates. Scripts in folder [2-TSSseq] performed following analyses:
1. 5’-UMI, "RNA 5'-end", and barcode extraction;
2. 5’-UMI, "RNA 5'-end", and barcode error correction using UMI-tools;
3. UMIs-based deduplication;
The final output is tab delimited files containing "RNA 5'-end"-Barcode linkages and coresponding deduplicated TSS-seq read counts.

Subsequently, DNA template and TSS usage information of each TSS promoter variant were merged based on barcode. Then, reads with “RNA 5’-end” sequence perfectly matched to corresponding template sequence were kept and used for downstream analysis. A TSS-seq count table containing TSS usage distribution of each TSS promoter variant was generated. Scripts in folder [3-DTmerge] also include downstream analysis and visualization.

## Genomic TSS-seq analysis
The purposed for genomic TSS-seq analysis is to investigate TSS sequence preference of natural promoters. Scripts in folder [4-genome] performed following analyses:
1. DNA context extraction for all positions within known promoter windows;
2. sequence preference analysis for "median" TSSs;
3. data preparation for later model prediction;

## Modeling analysis
The purposed for modeling analysis is to quantitatively identify key features (sequences and interactions) for TSS efficiency and evaluate what extent does DNA sequence around a TSS contribute to TSS efficiency in genomic promoters. Scripts in folder [5-modeling] performed following analyses:
1. regression modeling using a forward stepwise strategy with a 5-fold Cross-Validation;
2. prediction on genomic TSSs


