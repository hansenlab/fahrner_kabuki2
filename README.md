# fahrner_kabuki

Code and associated files for analysis of KS2 RNAseq data.
Raw and processed data files available at GEO accession number GSE262539.

The `RNAseq/final_analysis` directory contains scripts and data objects for Fig 4 and Supp Fig 6 and 8.
  - `RNAseq_final.R`: Imports the quant files, summarizes to transcript level using `tximeta` and performs differential analysis with `DESeq2`
  - `Fig4.Rmd`: Generates plots displayed in Fig 4 and Supp Fig 6 and 8
  - Several data objects `MGI.*.rda`: Gene lists obtained from MGI used in the Wilcoxon rank sum analyses
  - Several other data objects `permutation.rank.*.rda`: Rank permutations generated in the course of the Wilcoxon rank sum analyses
  - `quants_final/`: Salmon output files for the experiment Kabuki 2 Mut vs Wt cell lines at Day 14
  - `quants_final_KS1/`: Salmon output files for a previously published experiment Kabuki 1 Mut vs Wt at Day 7 (Fahrner et al., JCI Insight 2019), processed identically
  - `RNAseq_final_KS1.R`: Similar to `RNAseq_final.R` but for the previously published Kabuki 1 experiment

The subdirectory `RNAseq/final_analysis/230510_KS1-KS2_NGS` contains the script and sample sheet for Fig 5 and Supp Fig 9. This is the combined RNAseq experiment with both KS2 and KS1 cell lines, at Days 7 and 14 of chondrogenic differentiation.
  - `quants_final_KS1-KS2/`: Salmon output files for the combined RNAseq experiment
