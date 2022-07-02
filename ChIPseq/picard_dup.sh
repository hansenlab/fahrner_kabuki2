#$ -l mem_free=10G,h_vmem=10G
#$ -l h_fsize=50G
#$ -cwd
#$ -V
#$ -m e
#$ -M cgao13@jhmi.edu

echo "Start time: $(date)"

java -jar picard.jar -h

for FILE in `ls *{GRCm38,dm6}.bam | sed 's/.bam//g' | sort -u`

do

java -jar picard.jar MarkDuplicates I=${FILE}.bam O=${FILE}_filtered.bam M=${FILE}_dup_metrics.txt REMOVE_DUPLICATES=TRUE
	
done

echo "End time: $(date)"
