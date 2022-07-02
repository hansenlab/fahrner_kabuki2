#$ -l mem_free=5G,h_vmem=5G
#$ -l h_fsize=50G
#$ -cwd
#$ -V
#$ -m e
#$ -M cgao13@jhmi.edu

samtools --version

for FILE in `ls merge*chr.bam | sed 's/.bam//g'`
do
	samtools sort -o ${FILE}_sorted.bam ${FILE}.bam
	samtools index ${FILE}_sorted.bam
done
