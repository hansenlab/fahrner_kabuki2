#$ -l mem_free=5G,h_vmem=5G
#$ -l h_fsize=50G
#$ -cwd
#$ -V
#$ -m e
#$ -M cgao13@jhmi.edu

echo "Start time: $(date)"

conda activate macs2
echo $CONDA_DEFAULT_ENV

macs2 --version

echo "Now calling sample merge-WT, time is $(date)..."

macs2 callpeak \
        -t merge-WT_IP_GRCM38_chr_sorted.bam \
        -c merge-WT_Input_GRCM38_chr_sorted.bam \
        -f BAM \
        -g mm \
        --outdir MACS2-merge \
        -n merge-WT_GRCM38 \
        --broad

echo "Now calling sample merge-Mut, time is $(date)..."

macs2 callpeak \
        -t merge-Mut_IP_GRCM38_chr_sorted.bam \
        -c merge-Mut_Input_GRCM38_chr_sorted.bam \
        -f BAM \
        -g mm \
        --outdir MACS2-merge \
        -n merge-Mut_GRCM38 \
        --broad

echo "End time: $(date)" 

