#$ -l mem_free=5G,h_vmem=5G
#$ -l h_fsize=50G
#$ -cwd
#$ -V
#$ -m e
#$ -M cgao13@jhmi.edu

echo "Start time: $(date)"

bowtie2 --version

for FILE in `ls *.fastq.gz | sed 's/_L002_R[12]_001.fastq.gz//g' | sort -u`
do
echo "Aligning ${FILE}"
bowtie2 -q -p 16 -t --phred33 --sensitive \
	-x combined-index \
	-1 ${FILE}_L002_R1_001.fastq.gz \
	-2 ${FILE}_L002_R2_001.fastq.gz \
	-S ${FILE}_aligned.sam
done 

echo FINISHED ALL ALIGNMENTS

samtools --version

for ALIGNED in `ls *.sam | sed 's/_aligned.sam//g' | sort -u`
do
	echo ${ALIGNED}_aligned.sam
        samtools view -@ 16 -S -b ${ALIGNED}_aligned.sam > ${ALIGNED}_aligned.bam
        samtools sort -o ${ALIGNED}_aligned_sorted.bam \
                ${ALIGNED}_aligned.bam
        samtools index ${ALIGNED}_aligned_sorted.bam
        samtools view -@ 16 -b -h ${ALIGNED}_aligned_sorted.bam \
                -L dm6.bed > ${ALIGNED}_aligned_sorted_dm6.bam
        samtools view -@ 16 -b -h ${ALIGNED}_aligned_sorted.bam \
                -L GRCm38.bed > ${ALIGNED}_aligned_sorted_GRCm38.bam

done

echo "End time: $(date)"

