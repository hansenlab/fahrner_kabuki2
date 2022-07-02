#$ -l mem_free=50G,h_vmem=50G
#$ -l h_fsize=50G
#$ -cwd
#$ -V
#$ -m e
#$ -M cgao13@jhmi.edu

echo "Start time: $(date)"

bowtie2-build --version
bowtie2-build -f --threads 8 combined-ref.fa combined-index

echo "End time: $(date)"
