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

########################################### Early infection vs Other
pbmc_meta$InEvsOther <- pbmc_meta$Infect_time
pbmc_meta[which(pbmc_meta$Infect_time != 'Infect-Early'),'InEvsOther'] <- 'Control'
pbmc_meta[which(pbmc_meta$Infect_time == 'Infect-Early'),'InEvsOther'] <- 'Case'

description <- factor(pbmc_meta$InEvsOther, levels = c('Case','Control'))
matrix <- as.data.frame(pbmc)
design <- model.matrix(~description + 0, matrix)
colnames(design) <- c('Case','Control')

dge <- DGEList(counts = pbmc)
v <- voom(dge, design, normalize = 'quantile')
fit <- lmFit(v, design)

cont.matrix <- makeContrasts(Case-Control, levels = design)
fit2 = contrasts.fit(fit, cont.matrix)
fit2 = eBayes(fit2)

DEG = topTable(fit2, adjust.method = "fdr", sort.by = "B", number = nrow(matrix))
DEG = na.omit(DEG)

DEG <- data.frame(Gene = row.names(DEG),DEG)
diff <- subset(DEG,DEG$P.Val < 0.05 & abs(DEG$logFC) > 1.2)
diff$Trend <- "up"
diff$Trend[diff$logFC < 0] <- "down" 
#####################

####################################################### 火山图
library(EnhancedVolcano)
library(magrittr)

DEG$P.Val <- as.numeric(DEG$P.Val)

p1 <- EnhancedVolcano(DEG,
              lab = rownames(DEG),
              x = "logFC",
              y = "P.Val",
              xlim = c(-3,3),
              ylim = c(0,6),                
              title = NULL,
              subtitle = NULL,
              pCutoff = 0.05,
              pointSize = 4.0,
              labSize = 6,
              #shape = 8,
              colAlpha = 1,
              #drawConnectors = TRUE,
              widthConnectors = 0.75)

setwd('F:\\Bull Figure\\Figure2')
pdf(file = 'InfectEarly_vs_Other.pdf', width = 6.5, height = 6.5)
p1
dev.off()

##################################################  Early infection vs Late infection
pbmc <- pbmc[,1:12]
pbmc_meta <- pbmc_meta[1:12,]

description <- factor(pbmc_meta$InEvsOther, levels = c('Case','Control'))
matrix <- as.data.frame(pbmc)
design <- model.matrix(~description + 0, matrix)
colnames(design) <- c('Case','Control')

dge <- DGEList(counts = pbmc)
v <- voom(dge, design, normalize = 'quantile')
fit <- lmFit(v, design)

cont.matrix <- makeContrasts(Case-Control, levels = design)
fit2 = contrasts.fit(fit, cont.matrix)
fit2 = eBayes(fit2)

DEG = topTable(fit2, adjust.method = "fdr", sort.by = "B", number = nrow(matrix))
DEG = na.omit(DEG)

DEG <- data.frame(Gene = row.names(DEG),DEG)
diff <- subset(DEG,DEG$P.Val < 0.05 & abs(DEG$logFC) > 1.2)
diff$Trend <- "up"
diff$Trend[diff$logFC < 0] <- "down" 







