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

#group_list <- pbmc_meta$Infect
#batch <- pbmc_meta$type

#g=factor(group_list)
#g=relevel(g,'Uninfect')
#design=model.matrix(~g) 

#ex_b_limma <- removeBatchEffect(pbmc,
#                                batch = batch,
#                                design = design)

#pbmc <- ex_b_limma
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


