Run the scripts in the following order

```
index_build.sh
bowtie2_ChIPseq.sh
picard_dup.sh
merge-reheader.sh
merge_sort_index.sh
macs2-merge.sh
```

ALTERNATIVE, with comments

1. `index_build.sh`: builds a combined bowtie index for mouse and drosophila
2. `bowtie2_ChIPseq.sh`: aligns reads, sort and index, separate into mouse and drosophila `.bam` files
3. `picard_dup.sh`: marks and removes read duplicates (requires `picard.jar`)
4. `merge-reheader`: merge processed `.bam` files by genotype and rename chromosomes to chr format
5. `merge_sort_index.sh`: sort and index the merged `.bam` files
6. `macs2-merge.sh`: call peaks on merged IP/input files by genotype
