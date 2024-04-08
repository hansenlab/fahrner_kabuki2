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
plot <- plotPCA(vst(dds))

# setEPS()
postscript('221016_KS2_PCA-uncollapsed.eps')
plot + geom_label_repel(aes(label=gsub("_.*", "", colnames(dds))), show.legend = FALSE)
# dev.off()

# collapse technical replicates
dds$sample <- factor(gsub("_.*", "", colnames(dds)))
dds.coll <- collapseReplicates(dds, dds$sample)
stopifnot(all(rowSums(counts(dds[,which(dds$sample=="10-29")])) == counts(dds.coll[,1])))
dds.coll.counts <- counts(dds.coll)

# filter genes with low median counts
keep <- rowMedians(counts(dds.coll)) > 10
dds.coll.filtered <- dds.coll[keep,]

#dds.coll.filtered$condition <- relevel(dds.coll.filtered$condition, ref="WT")
plot <- plotPCA(vst(dds.coll.filtered), returnData=TRUE)

# setEPS()
postscript('221016_KS2_PCA-collapsed.eps', width=7, height=4)
plot + geom_label_repel(aes(label=gsub("_.*", "", colnames(dds.coll.filtered))), show.legend = FALSE)
# dev.off()

# find surrogate variables for batch effects using sva
dds.coll.filtered$condition <- relevel(dds.coll.filtered$condition, ref="WT")
dds.coll.filtered.counts <- counts(dds.coll.filtered)
# mod <- model.matrix(~condition, colData(dds.coll.filtered))
# mod0 <- model.matrix(~1, colData(dds.coll.filtered))
# svseq <- svaseq(dds.coll.filtered.counts, mod, mod0)
# svseq2 <- svaseq(dds.coll.filtered.counts, mod, mod0)
# svseq <- svaseq(dds.coll.filtered.counts, mod, mod0, n.sv=1)

# append surrogate variables
# ddssva <- dds.coll.filtered
# ddssva$SV1 <- svseq$sv[,1]
# ddssva$SV2 <- svseq$sv[,2]
# ddssva$SV3 <- svseq$sv[,3]
# design(ddssva) <- formula(~SV1 + SV2 + SV3 + condition)
# design(ddssva) <- formula(~SV1 + SV2 + condition)
# design(ddssva) <- formula(~SV1 + condition)

# DESeq2 call
# ddssva <- DESeq(ddssva)
dds <- DESeq(dds.coll.filtered)
# dds <- ddssva
# save(dds, file='dds_final_no-sva.rda')
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
load(file='MGI.chon.rda')
diffexp.subset$chondrocyte <- 
  ifelse(toupper(diffexp.subset$mgi_symbol) %in% MGI.chon$gene_name, 'yes','no')
write.csv(diffexp.subset, file="diffexp_final_no-sva_chon.csv")


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

### PLOTS
# Wilcoxon rank-sum, chondrogenesis genes
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)
edb_mouse <- EnsDb.Mmusculus.v79
proms_mouse <- promoters(edb_mouse, filter = TxBiotypeFilter("protein_coding"), upstream = 2000, downstream = 2000, columns = c("gene_name", "tx_id", "gene_id"))
genome(seqinfo(proms_mouse)) <- "mm10"
seqlevelsStyle(proms_mouse) <- "ucsc"
proms_mouse <- proms_mouse[which(seqnames(proms_mouse) %in% seqnames(Mmusculus)[1:21])]
proms_mouse$gene_name <- toupper(proms_mouse$gene_name)

genes_chon <- read_csv('~/Desktop/temp/Fahrner lab/Data/KS2/220506_NGS/fahrner_kabuki2/RNAseq/final_analysis/chondrogenesis_genes.csv', col_names = FALSE)
genes_chon <- genes_chon$X1
genes_chon <- toupper(genes_chon)
# genes_chon <- genes_chon[-which(duplicated(genes_chon))]
gene_ids <- unique(proms_mouse$gene_id[which(proms_mouse$gene_name %in% genes_chon)])

rank_KS2 <- wilcox.test(res$pvalue[which(rownames(res) %in% gene_ids)], 
                       res$pvalue[-which(rownames(res) %in% gene_ids)])$statistic

length1 <- length(which(rownames(res) %in% gene_ids))

permutation_rank_chon <- replicate(10000, {
  indices <- sample(1:length(rownames(res)), length1)
  wilcox.test(res$pvalue[indices], res$pvalue[-indices])$statistic
})

save(permutation_rank_chon, file='permutation_rank_chon.rda')
prop.table(table(permutation_rank_chon<rank_KS2))

setEPS()
postscript('221016_chon.eps', width=7, height=6)
par(mar=c(5,5,6,5))
hist(permutation_rank_chon, col = "cornflowerblue", lty = 0, 
     breaks = 50, freq = FALSE, xlab = "Rank-sum statistic", cex.lab = 1.5, yaxt = 'n',
     main = "Chondrogenesis genes", cex.main = 1.5, font.main = 1, xlim = c(rank_KS2-20000, max(permutation_rank_chon)+0.05), xaxt = 'n')
axis(1, at = c(round(quantile(permutation_rank_chon, c(0.01, 0.99))),rank_KS2), cex.axis = 1.5)
axis(2, at = c(0, 0.000011), cex.axis = 1.5)
abline(v = rank_KS2, col = "red", lwd = 2.5)
legend <- legend("topright", legend = c("random", "observed"), 
                 col = c("cornflowerblue", "red"), bty = 'n', 
                 cex = 1.5, lty = "solid", lwd = 2.5)
dev.off()

genes_ost <- read.csv(file='osteogenesis_genes.csv')
gene_ids <- unique(genes_ost$Gene.ID)
rank_ost <- wilcox.test(res$pvalue[which(rownames(res) %in% gene_ids)], 
                        res$pvalue[-which(rownames(res) %in% gene_ids)])$statistic

length1 <- length(which(rownames(res) %in% gene_ids))

permutation_rank_ost <- replicate(10000, {
  indices <- sample(1:length(rownames(res)), length1)
  wilcox.test(res$pvalue[indices], res$pvalue[-indices])$statistic
})
save(permutation_rank_ost, file='permutation_rank_ost.rda')


makeVolcanoPlot <- function(DEgenes_df, lfc_cutoff, gene_ids, color){
  tab = data.frame(logFC = DEgenes_df$log2FoldChange, negLogPval = -log10(DEgenes_df$pvalue))
  rownames(tab) <- rownames(DEgenes_df)
  par(mar=c(5,5,5,5))
  plot(tab[-which(rownames(tab) %in% gene_ids), ], pch = 16, cex = 0.75, 
       xlab = expression(log[2]~fold~change),
       ylab = expression(-log[10]~pvalue), 
       cex.axis=1.5,
       cex.lab=1.5,
       col = "gray60", bty = 'l', 
       xlim = c(min(tab$logFC), max(tab$logFC)), 
       ylim = c(0, max(-log10(res$pvalue)[which(is.na(-log10(res$pvalue))==FALSE)])))
  points(tab[gene_ids, ], pch = 19, cex = 0.75, col = color)
  abline(h = -log10(0.00551414), col= "red", lty=2, lwd=2)
  # legend <- legend('topright', legend=c('FDR = 0.1','chondrogenesis gene'),
  #                  col=c('red','#04C6F7'), lty=c(2, NA), pch=c(NA,19), lwd=c(2,NA), cex=1.25, bty='n')
  # legend <- legend('topright', legend=c('FDR = 0.1'),
  #                  col=c('cornflowerblue'), lty=c(2), pch=c(NA), lwd=c(2), cex=1.25, bty='n')
  if (lfc_cutoff != "none"){
    abline(v = c(-lfc_cutoff, lfc_cutoff), col = rgb(0,0,0,0.75), lty = "longdash", lwd = 1) 
  }
}

setEPS()
postscript('221016_KS2_volcano-gray.eps', width=7, height=6)
makeVolcanoPlot(res, "none", gene_ids, "red") #gene_ids have been set to the chondrogenesis gene ids
dev.off()

fdr <- res[which(res$padj>0.1 & is.na(res$padj)==FALSE),]
fdr[which(fdr$padj==min(fdr$padj)),]

plot <- plotPCA(vst(dds.coll.filtered), returnData=TRUE)
plot$group <- relevel(plot$group, ref='WT')
group.colors <- c('cornflowerblue','red')
setEPS()
postscript('221023_KS2_PCA-collapsed.eps', width=9, height=5)
ggplot(data = plot, aes(x = PC1, y = PC2, color=group)) + 
  geom_point(size = 4) + 
  scale_color_manual(labels=c('+/+','-/-'), values=group.colors) + 
  xlab('PC1: 44% variance') +
  ylab('PC2: 18% variance') +
  theme_bw(base_size=18)
dev.off()

