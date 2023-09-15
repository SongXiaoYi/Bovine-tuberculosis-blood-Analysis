####### Bovine-tuberculosis Blood
library(dplyr)
library(stringr)
library(patchwork)
library(knitr)
library(msigdbr) 
library(fgsea)
library(ggplot2)
library(GSVA)
library(GSEABase)
library(limma)
library(clusterProfiler)
library(DESeq2)
library(ggplotify)
library(cowplot)
library(edgeR)
library(tinyarray)
library(gridExtra)
library(Seurat)
######################################
setwd('F:\\Bull Figure')
############################################### Matrix
pbmc <- read.csv('./Normalized_data.csv')
pbmc <- as.matrix(pbmc)
rownames(pbmc) <- pbmc[,1]
Rownames <- pbmc[,1]
pbmc <- pbmc[,-1]
pbmc <- apply(pbmc,2,as.numeric)
rownames(pbmc) <- Rownames
pbmc <- avereps(pbmc)
pbmc <- pbmc[rowMeans(pbmc)>0,]

############################################ meta data
pbmc_meta <- readRDS('./datTraits.rds')
pbmc_meta$id <- paste0('X',pbmc_meta$id)
pbmc_meta$id <- gsub('-','.',pbmc_meta$id)
rownames(pbmc_meta) <- pbmc_meta$id
pbmc_meta$Infect <- c(rep('Uninfect',6),rep('Infect',6),rep('Uninfect',6),rep('Infect',6))
pbmc_meta$Infect_time <- paste(pbmc_meta$Infect,pbmc_meta$type,sep = "-")
pbmc_meta$Infect <- factor(pbmc_meta$Infect, levels = c('Uninfect','Infect'))
############################################################## Uninfect-Infect
############################################ DEGeq2
#colData <- data.frame(row.names = colnames(pbmc),
#                     condition = pbmc_meta$Infect)

#DEG_matrix <- apply(pbmc,2,as.integer)
#rownames(DEG_matrix) <- rownames(pbmc)

#dds <- DESeqDataSetFromMatrix(
#  countData = DEG_matrix,
#  colData = colData,
#  design = ~ condition
#)

#dds <- DESeq(dds)

#res <- results(dds, contrast = c("condition", rev(levels(pbmc_meta$Infect))))
#resOrdered <- res[order(res$pvalue),]
#DEG <- as.data.frame(resOrdered)
#DEG <- na.omit(DEG)

#logFC_cutoff <- with(DEG,mean(abs(log2FoldChange)) + 2*sd(abs(log2FoldChange)))
#k1 <- (DEG$pvalue < 0.05) & (DEG$log2FoldChange < -logFC_cutoff)
#k2 <- (DEG$pvalue < 0.05) & (DEG$log2FoldChange > logFC_cutoff)
#DEG$change <- ifelse(k1,"DOWN",ifelse(k2,"UP","NOT"))
#table(DEG$change)

#DESeq2_DEG <- DEG
############################################ edgeR
dge <- DGEList(counts = pbmc, group = pbmc_meta$Infect)
dge$samples$lib.size <- colSums(dge$counts)
#######################
dge <- calcNormFactors(dge) 

design <- model.matrix(~0 + pbmc_meta$Infect)
rownames(design) <- colnames(dge)
colnames(design) <- levels(pbmc_meta$Infect)

dge <- estimateGLMCommonDisp(dge, design)
dge <- estimateGLMTrendedDisp(dge, design)
dge <- estimateGLMTagwiseDisp(dge, design)

fit <- glmFit(dge, design)
fit2 <- glmLRT(fit, contrast = c(-1,1))

DEG <- topTags(fit2, n = nrow(pbmc))
DEG <- as.data.frame(DEG) 

logFC_cutoff <- with(DEG,mean(abs(logFC)) + 2*sd(abs(logFC)))
k1 = (DEG$PValue < 0.05)&(DEG$logFC < -logFC_cutoff)
k2 = (DEG$PValue < 0.05)&(DEG$logFC > logFC_cutoff)
DEG$change <- ifelse(k1, "DOWN", ifelse(k2, "UP","NOT"))
table(DEG$change)

edgeR_DEG <- DEG

############################################# limma
pbmc_meta$Infect <- factor(pbmc_meta$Infect,levels = c('Infect','Uninfect'))
design <- model.matrix(~0 + pbmc_meta$Infect)
colnames(design) <- levels(pbmc_meta$Infect)
rownames(design) <- colnames(pbmc)

dge <- DGEList(counts = pbmc)
dge <- calcNormFactors(dge)

v <- voom(dge, design, normalize = 'quantile')
fit <- lmFit(v, design)

constrasts <- paste(rev(levels(pbmc_meta$Infect)), collapse = "-")
cont.matrix <- makeContrasts(contrasts = constrasts, levels = design)
fit2 = contrasts.fit(fit, cont.matrix)
fit2 = eBayes(fit2)

DEG = topTable(fit2, coef = constrasts, n = Inf)
DEG = na.omit(DEG)
##################################################################### cutoff
logFC_cutoff <- with(DEG,mean(abs(logFC)) + 2*sd(abs(logFC)))
k1 <- (DEG$P.Value < 0.05) & (DEG$logFC < -logFC_cutoff)
k2 <- (DEG$P.Value < 0.05) & (DEG$logFC > logFC_cutoff)
DEG$change <- ifelse(k1, "DOWN", ifelse(k2, "UP", "NOT"))
table(DEG$change)

limma_DEG <- DEG

##################################################################### PCA
########################## Early-Late
pbmc_meta$Infect <- factor(pbmc_meta$Infect,levels = c('Uninfect','Infect'))
pca.plot <- draw_pca(pbmc,pbmc_meta$Infect)
p1 <- pca.plot + theme_bw() + NoLegend()

setwd('F:\\Bull Figure\\Figure1')
pdf(file = 'PCA-plot_p1_legend.pdf', width = 3.5, height = 3.55)
pca.plot + theme_bw()
dev.off()

########################## 
pbmc_meta$Infect_time <- factor(pbmc_meta$Infect_time,levels = c('Uninfect-Early','Infect-Early','Uninfect-Late','Infect-Late'))
pca.plot <- draw_pca(pbmc,pbmc_meta$Infect_time,color = c("#E78AC3","#A6D854","#868686","#66C2A5"))
p2 <- pca.plot + theme_bw() + NoLegend()

setwd('F:\\Bull Figure\\Figure1')
pdf(file = 'PCA-plot_p2_legend.pdf', width = 3.5, height = 3.55)
pca.plot + theme_bw()
dev.off()

setwd('F:\\Bull Figure\\Figure1')
pdf(file = 'PCA-plot.pdf', width = 7, height = 3.55)
grid.arrange(p1,p2,ncol=2)
dev.off()


##################################################################### Intersection
UP_function <- function(df){
  rownames(df)[df$change == 'UP']
}

DOWN_function <- function(df){
  rownames(df)[df$change == 'DOWN']
}

Up_all <- intersect(UP_function(limma_DEG),UP_function(edgeR_DEG))
DOWN_all <- intersect(DOWN_function(edgeR_DEG_EL),DOWN_function(limma_DEG_EL))


