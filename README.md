# fahrner_kabuki

Code and associated files for analysis of KS2 RNAseq data.

The `RNAseq/final_analysis` directory contains scripts and data objects for Fig 4 and Supp Fig 6 and 8.
  - `RNAseq_final.R`: Imports the quant files from Salmon (available on GEO accession number GSE262539), summarizes to transcript level using `tximeta` and performs differential analysis with `DESeq2`
  - `Fig4.Rmd`: Generates plots displayed in Fig 4 and Supp Fig 6 and 8
  - Several data objects `MGI.*.rda`: Gene lists obtained from MGI used in the Wilcoxon rank sum analyses
  - Several other data objects `permutation.rank.*.rda`: Rank permutations generated in the course of the Wilcoxon rank sum analyses

The subdirectory `RNAseq/final_analysis/230510_KS1-KS2_NGS` contains the script and sample sheet for Fig 5 and Supp Fig 9. This is the combined RNAseq experiment with both KS2 and KS1 cell lines, at Days 7 and 14 of chondrogenic differentiation.
