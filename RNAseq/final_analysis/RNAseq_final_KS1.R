library(TxDb.Mmusculus.UCSC.mm10.ensGene)
library(tximport)
library(DESeq2)
library(sva)
library(biomaRt)
library(ggplot2)
library(ggrepel)
library(tximeta)
library(SummarizedExperiment)
library(edgeR)
library(limma)
library(statmod)
library(tidyverse)


### BEGIN IMPORT
samples <- read.csv('RNAseq-samples_KS1.csv')
samples$genotype <- ifelse(str_detect(samples$Sample.ID, 'DH'), 'Mut', 'Wt')
samples$Sample.ID <- str_replace_all(samples$Sample.ID, c('Kmt2d '='', ' D7_[12]'=''))

dir <- '~/Desktop/temp/Fahrner lab/Data/KS2/220506_NGS/fahrner_kabuki2/RNAseq/final_analysis/quants_final_KS1'
files <- file.path(dir, paste(samples$fasta, '_quant', sep=''), 'quant.sf')

file.exists(files)
coldata <- data.frame(files, names=samples$Sample.ID, condition=samples$genotype, stringsAsFactors = FALSE)

# import transcript abundances
se <- tximeta(coldata)

# summarize abundances to gene level
gse <- summarizeToGene(se)

# load SummarizedExperiment into DESeqDataSet
dds <- DESeqDataSet(gse, design=~condition)
### END IMPORT


### BEGIN ANALYSIS 
plot <- plotPCA(vst(dds))

# setEPS()
# postscript('221016_KS2_PCA-uncollapsed.eps')
plot + geom_label_repel(aes(label=gsub("_.*", "", colnames(dds))), show.legend = FALSE)
# dev.off()

# collapse technical replicates
dds$sample <- colnames(dds)
dds.coll <- collapseReplicates(dds, dds$sample)
stopifnot(all(rowSums(counts(dds[,which(dds$sample=="Con 10-5")])) == counts(dds.coll[,1])))
dds.coll.counts <- counts(dds.coll)

# filter genes with low median counts
keep <- rowMedians(counts(dds.coll)) > 10
dds.coll.filtered <- dds.coll[keep,]

#dds.coll.filtered$condition <- relevel(dds.coll.filtered$condition, ref="WT")
plot <- plotPCA(vst(dds.coll.filtered), returnData=TRUE)
# plot + geom_label_repel(aes(label=colnames(dds.coll.filtered)), show.legend = FALSE)
plot$group <- relevel(plot$group, ref='Wt')
group.colors <- c('deepskyblue','red')
setEPS()
postscript('PCA_KS1.eps', width=9, height=5)
ggplot(data = plot, aes(x = PC1, y = PC2, color=group)) + 
  geom_point(size = 4) + 
  scale_color_manual(labels=c('Kmt2d +/+','Kmt2d -/-'), values=group.colors) + 
  xlab('PC1: 59% variance') +
  ylab('PC2: 25% variance') +
  theme_bw(base_size=18) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
dev.off()

# find surrogate variables for batch effects using sva
dds.coll.filtered$condition <- relevel(dds.coll.filtered$condition, ref="Wt")
dds.coll.filtered.counts <- counts(dds.coll.filtered)
mod <- model.matrix(~condition, colData(dds.coll.filtered))
mod0 <- model.matrix(~1, colData(dds.coll.filtered))
svseq <- svaseq(dds.coll.filtered.counts, mod, mod0)

# append surrogate variables
ddssva <- dds.coll.filtered
ddssva$SV1 <- svseq$sv[,1]
design(ddssva) <- formula(~SV1 + condition)

# DESeq2 call
ddssva <- DESeq(ddssva)
dds <- ddssva
save(dds, file='BigStuff/dds_KS1.rda')
res <- results(dds)

diffexp.subset <- as.data.frame(res[which(res$padj <0.1),])
### END ANALYSIS

# gene annotation
ensembl.id <- rownames(diffexp.subset)

ensembl102.mmusculus <- useEnsembl(biomart='genes', dataset='mmusculus_gene_ensembl', version=102)
mgi <- getBM(attributes=c('ensembl_gene_id','mgi_symbol','mgi_description'),
             filters='ensembl_gene_id',
             values=ensembl.id,
             mart=ensembl102.mmusculus)
diffexp.subset$ensembl_gene_id <- ensembl.id
diffexp.subset <- merge(diffexp.subset, mgi, by='ensembl_gene_id')
diffexp.subset <- diffexp.subset[order(diffexp.subset$log2FoldChange, decreasing=TRUE),]
# write.csv(diffexp.subset, file="diffexp_final_no-sva.csv")