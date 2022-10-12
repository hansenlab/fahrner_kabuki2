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
samples <- read.csv('RNAseq-samples.csv')

dir <- '~/Desktop/temp/Fahrner lab/Data/KS2/220506_NGS/fahrner_kabuki2/RNAseq/final_analysis/quants_final'
files <- file.path(dir, paste(samples$Sample.ID, '_quant', sep=''), 'quant.sf')

file.exists(files)
coldata <- data.frame(files, names=samples$Sample.ID, condition=samples$Genotype, stringsAsFactors = FALSE)

# import transcript abundances
se <- tximeta(coldata)

# summarize abundances to gene level
gse <- summarizeToGene(se)

# load SummarizedExperiment into DESeqDataSet
dds <- DESeqDataSet(gse, design=~condition)
### END IMPORT


### BEGIN ANALYSIS 
dds$condition <- relevel(dds$condition, ref="WT")
plot <- plotPCA(vst(dds))
plot + geom_label_repel(aes(label=gsub("_.*", "", colnames(dds))), show.legend = FALSE)

# collapse technical replicates
dds$sample <- factor(gsub("_.*", "", colnames(dds)))
dds.coll <- collapseReplicates(dds, dds$sample)
stopifnot(all(rowSums(counts(dds[,which(dds$sample=="10-29")])) == counts(dds.coll[,1])))
dds.coll.counts <- counts(dds.coll)

# filter genes with low median counts
keep <- rowMedians(counts(dds.coll)) > 10
dds.coll.filtered <- dds.coll[keep,]

#dds.coll.filtered$condition <- relevel(dds.coll.filtered$condition, ref="WT")
plot <- plotPCA(vst(dds.coll.filtered))
plot + geom_label_repel(aes(label=gsub("_.*", "", colnames(dds.coll.filtered))), show.legend = FALSE)

# # find surrogate variables for batch effects using sva
# dds.coll.filtered.counts <- counts(dds.coll.filtered)
# mod <- model.matrix(~condition, colData(dds.coll.filtered))
# mod0 <- model.matrix(~1, colData(dds.coll.filtered))
# svseq <- svaseq(dds.coll.filtered.counts, mod, mod0)
# 
# # append surrogate variables
# ddssva <- dds.coll.filtered
# ddssva$SV1 <- svseq$sv[,1]
# ddssva$SV2 <- svseq$sv[,2]
# ddssva$SV3 <- svseq$sv[,3]
# design(ddssva) <- formula(~SV1 + SV2 + SV3 + condition)
# #design(ddssva) <- formula(~SV1 + SV2 + condition)

# DESeq2 call
#ddssva <- DESeq(ddssva)
dds <- DESeq(dds.coll.filtered)
#dds <- ddssva
save(dds, file='dds_final_no-sva.rda')
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
write.csv(diffexp.subset, file="diffexp_final.csv")


### TEST SECTION
load(file='dds_final.rda')
dds_sva <- dds
load(file='dds_final_no-sva.rda')
dds_no.sva <- dds
rm(dds)

res_sva <- results(dds_sva)
res_no.sva <- results(dds_no.sva)

comp <- merge(as.data.frame(res_sva), as.data.frame(res_no.sva), by="row.names", suffixes=c('.sva','.no.sva'))

setEPS()
postscript('sva_comparison.eps')
plot(comp$log2FoldChange.sva, comp$log2FoldChange.no.sva,
     pch=19, cex=0.5,
     #xlim=c(-8.5,3), ylim=c(-8.5,3),
     xlab='log2FC, sva', ylab='log2FC, no sva', main='sva vs no sva, all genes')
abline(a=0, b=1, col='red', lty=2)
abline(h=0)
abline(v=0)
dev.off()


