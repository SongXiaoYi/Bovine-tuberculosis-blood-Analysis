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

############################################################## Early-Late
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
k1 = (DEG$PValue < 0.05)&(DEG$logFC < -logFC_cutoff)
k2 = (DEG$PValue < 0.05)&(DEG$logFC > logFC_cutoff)
DEG$change <- ifelse(k1, "DOWN", ifelse(k2, "UP","NOT"))
table(DEG$change)

edgeR_DEG_EL <- DEG

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
k1 <- (DEG$P.Value < 0.05) & (DEG$logFC < -logFC_cutoff)
k2 <- (DEG$P.Value < 0.05) & (DEG$logFC > logFC_cutoff)
DEG$change <- ifelse(k1, "DOWN", ifelse(k2, "UP", "NOT"))
table(DEG$change)

limma_DEG_EL <- DEG


############################################################## Uninfect-Infect
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

logFC_cutoff <- 1
k1 = (DEG$PValue < 0.05)&(DEG$logFC < -logFC_cutoff)
k2 = (DEG$PValue < 0.05)&(DEG$logFC > logFC_cutoff)
DEG$change <- ifelse(k1, "DOWN", ifelse(k2, "UP","NOT"))
table(DEG$change)

edgeR_DEG_Infect <- DEG

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
logFC_cutoff <- 1
k1 <- (DEG$P.Value < 0.05) & (DEG$logFC < -logFC_cutoff)
k2 <- (DEG$P.Value < 0.05) & (DEG$logFC > logFC_cutoff)
DEG$change <- ifelse(k1, "DOWN", ifelse(k2, "UP", "NOT"))
table(DEG$change)

limma_DEG_Infect <- DEG
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

##################################################################### Heatmap
######################################## Early-Late
heat_data <- as.data.frame(pbmc)
cg1_EL <- rownames(edgeR_DEG_EL)[edgeR_DEG_EL$change != "NOT"]
cg2_EL <- rownames(limma_DEG_EL)[limma_DEG_EL$change != "NOT"]

h1_EL <- draw_heatmap(heat_data[cg1_EL,],pbmc_meta$type,n_cutoff = 2)
h2_EL <- draw_heatmap(heat_data[cg2_EL,],pbmc_meta$type,n_cutoff = 2)

######################################## Infect-Uninfect
heat_data <- as.data.frame(pbmc)
cg1_Infect <- rownames(edgeR_DEG_Infect)[edgeR_DEG_Infect$change != "NOT"]
cg2_Infect <- rownames(limma_DEG_Infect)[limma_DEG_Infect$change != "NOT"]

h1_Infect <- draw_heatmap(heat_data[cg1_Infect,],pbmc_meta$Infect,n_cutoff = 2)
h2_Infect <- draw_heatmap(heat_data[cg2_Infect,],pbmc_meta$Infect,n_cutoff = 2)

##################################################################### Intersection
UP_function <- function(df){
  rownames(df)[df$change == 'UP']
}

DOWN_function <- function(df){
  rownames(df)[df$change == 'DOWN']
}

#################################################### Early-Late
Up_EL <- intersect(UP_function(edgeR_DEG_EL),UP_function(limma_DEG_EL))
DOWN_EL <- intersect(DOWN_function(edgeR_DEG_EL),DOWN_function(limma_DEG_EL))

#################################################### Infect-Uninfect
Up_Infect <- intersect(UP_function(edgeR_DEG_Infect),UP_function(limma_DEG_Infect))
DOWN_Infect <- intersect(DOWN_function(edgeR_DEG_Infect),DOWN_function(limma_DEG_Infect))

#################################################### Heatmap
hp <- draw_heatmap(heat_data[c(Up_Infect,DOWN_Infect),],pbmc_meta$Infect,n_cutoff = 2,annotation_legend = TRUE)

############################################### Venn
Up_genes <- list(edgeR = UP_function(edgeR_DEG_Infect),
                limma = UP_function(limma_DEG_Infect))

DOWN_genes <- list(edgeR = DOWN_function(edgeR_DEG_Infect),
                limma = DOWN_function(limma_DEG_Infect))

up.plot <- draw_venn(Up_genes, 'UPgene')
down.plot <- draw_venn(DOWN_genes, 'DOWNgene')







