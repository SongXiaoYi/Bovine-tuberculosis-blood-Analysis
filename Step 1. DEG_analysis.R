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

############################################ edgeR
dge <- DGEList(counts = pbmc, group = pbmc_meta$type)
dge$samples$lib.size <- colSums(dge$counts)
#######################
dge <- calcNormFactors(dge) 

design <- model.matrix(~0 + pbmc_meta$type)
rownames(design) <- colnames(dge)
colnames(design) <- levels(pbmc_meta$type)

dge <- estimateGLMCommonDisp(dge, design)
dge <- estimateGLMTrendedDisp(dge, design)
dge <- estimateGLMTagwiseDisp(dge, design)

fit <- glmFit(dge, design)
fit2 <- glmLRT(fit, contrast = c(-1,1))

DEG <- topTags(fit2, n = nrow(pbmc))
DEG <- as.data.frame(DEG) 

logFC_cutoff <- with(DEG, mean(abs(logFC)) + 2*sd(abs(logFC)))
k1 = (DEG$PValue < 0.05)&(DEG$logFC < -logFC_cutoff)
k2 = (DEG$PValue < 0.05)&(DEG$logFC > logFC_cutoff)
DEG$change <- ifelse(k1, "DOWN", ifelse(k2, "UP","NOT"))

edgeR_DEG <- DEG









