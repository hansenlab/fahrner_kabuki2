#$ -l mem_free=5G,h_vmem=5G
#$ -l h_fsize=50G
#$ -cwd
#$ -V
#$ -m e
#$ -M cgao13@jhmi.edu

samtools --version

echo "starting merge 1, $(date)"

samtools merge -o merge-WT_IP_GRCM38.bam \
	10-3_IP_S32_aligned_sorted_GRCm38_filtered.bam \
	10-29_IP_S29_aligned_sorted_GRCm38_filtered.bam \
	96-13_IP_S20_aligned_sorted_GRCm38_filtered.bam \
	96-14_IP_S30_aligned_sorted_GRCm38_filtered.bam \
	96-19_IP_S24_aligned_sorted_GRCm38_filtered.bam

echo "starting merge 2, $(date)"

samtools merge -o merge-WT_Input_GRCM38.bam \
	10-3_Input_S27_aligned_sorted_GRCm38_filtered.bam \
        10-29_Input_S35_aligned_sorted_GRCm38_filtered.bam \
        96-13_Input_S25_aligned_sorted_GRCm38_filtered.bam \
        96-14_Input_S22_aligned_sorted_GRCm38_filtered.bam \
        96-19_Input_S23_aligned_sorted_GRCm38_filtered.bam 

echo "starting merge 3, $(date)"

samtools merge -o merge-Mut_IP_GRCM38.bam \
	10-37_IP_S37_aligned_sorted_GRCm38_filtered.bam \
        10-39_IP_S33_aligned_sorted_GRCm38_filtered.bam \
        96-15_IP_S21_aligned_sorted_GRCm38_filtered.bam \
        96-20_IP_S39_aligned_sorted_GRCm38_filtered.bam \
        96-21_IP_S28_aligned_sorted_GRCm38_filtered.bam

echo "starting merge 4, $(date)"

samtools merge -o merge-Mut_Input_GRCM38.bam \
	10-37_Input_S26_aligned_sorted_GRCm38_filtered.bam \
        10-39_Input_S36_aligned_sorted_GRCm38_filtered.bam \
        96-15_Input_S31_aligned_sorted_GRCm38_filtered.bam \
        96-20_Input_S34_aligned_sorted_GRCm38_filtered.bam \
        96-21_Input_S38_aligned_sorted_GRCm38_filtered.bam

echo "done with merges, $(date)"

for FILE in `ls merge-*.bam | sed 's/.bam//g'`; do
	samtools view -H ${FILE}.bam | \
	sed -e 's/SN:\([0-9XY]\)/SN:chr\1/' -e 's/SN:MT/SN:chrM/' | \
	samtools reheader - ${FILE}.bam > ${FILE}_chr.bam
done
