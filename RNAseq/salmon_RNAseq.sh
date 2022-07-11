#!/bin/bash

echo "Start time: $(date)"

echo "Building index"
salmon index -t Mus_musculus.GRCm38.cdna.all.fa.gz -i GRCm38_cdna_index

echo "Mapping reads"
for FILE in `ls *.fastq.gz | sed 's/_L001_R[12]_001.fastq.gz//g' | sort -u`
do
	salmon quant -i GRCm38_cdna_index -l A \
		-1 ${FILE}_L001_R1_001.fastq.gz \
		-2 ${FILE}_L001_R2_001.fastq.gz \
		-p 8 --validateMappings --gcBias \
		-o quants/${FILE}_quant
done

echo "End time: $(date)"
