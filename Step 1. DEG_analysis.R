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

logFC_cutoff <- 1
k1 = (DEG$FDR < 0.05)&(DEG$logFC < -logFC_cutoff)
k2 = (DEG$FDR < 0.05)&(DEG$logFC > logFC_cutoff)
DEG$change <- ifelse(k1, "DOWN", ifelse(k2, "UP","NOT"))

edgeR_DEG <- DEG

############################################# limma
pbmc_meta$type <- factor(pbmc_meta$type,levels = c('Late','Early'))
design <- model.matrix(~0 + pbmc_meta$type)
colnames(design) <- levels(pbmc_meta$type)
rownames(design) <- colnames(pbmc)

dge <- DGEList(counts = pbmc)
dge <- calcNormFactors(dge)

v <- voom(dge, design, normalize = 'quantile')
fit <- lmFit(v, design)

constrasts <- paste(rev(levels(pbmc_meta$type)), collapse = "-")
cont.matrix <- makeContrasts(contrasts = constrasts, levels = design)
fit2 = contrasts.fit(fit, cont.matrix)
fit2 = eBayes(fit2)

DEG = topTable(fit2, coef = constrasts, n = Inf)
DEG = na.omit(DEG)
##################################################################### cutoff
logFC_cutoff <- 1
k1 <- (DEG$adj.P.Val < 0.05) & (DEG$logFC < -logFC_cutoff)
k2 <- (DEG$adj.P.Val < 0.05) & (DEG$logFC > logFC_cutoff)
DEG$change <- ifelse(k1, "DOWN", ifelse(k2, "UP", "NOT"))
table(DEG$change)

limma_DEG <- DEG
##################################################################### PCA
########################## Early-Late
pbmc_meta$type <- factor(pbmc_meta$type,levels = c('Early','Late'))
pca.plot <- draw_pca(pbmc,pbmc_meta$type)
p1 <- pca.plot + theme_bw()

########################## 
pbmc_meta$Infect <- factor(pbmc_meta$Infect,levels = c('Uninfect','Infect'))
pca.plot <- draw_pca(pbmc,pbmc_meta$Infect,color = c("#E78AC3","#A6D854"))
p2 <- pca.plot + theme_bw()

setwd('F:\\Bull Figure\\Figure1')
pdf(file = 'PCA-plot.pdf', width = 8.7, height = 3.4)
grid.arrange(p1,p2,ncol=2)
dev.off()

##################################################################### heatmap













