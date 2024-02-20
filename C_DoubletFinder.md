## Remove doublets using DoubletFinder

### CS
```
## load modules before starting R session
# module purge && module load shared slurm python/3.7.x-anaconda
# module load gcc/8.3.0
# module load hdf5_18/1.8.17
# module load R/4.0.2-gccmkl


##------------------------------------
## CB+DF FILTERED COUNT TABLES FOR CS
## load libraries
rm(list = ls())
library(dplyr)
library(patchwork)
library(ggplot2)
library(reshape2)
library(clustree)
library(Matrix)
library(Seurat)
library(SeuratObject)
set.seed(10)

## reference gene symbol list
ref <- read.table("../gencode.vM17.protein_coding_gene_id_symbol.txt.gz", header = TRUE, sep = "\t")

## load cellbender retained cellbarcodes list
cs_a1_cb <- scan("/work/psychiatry/rgree2/SNRNASEQ/03_CELLBENDER/CS_A1/CellBender_Out_cell_barcodes.csv", what = "", sep = "\n") # 19998
cs_a2_cb <- scan("/work/psychiatry/rgree2/SNRNASEQ/03_CELLBENDER/CS_A2/CellBender_Out_cell_barcodes.csv", what = "", sep = "\n") # 11160
cs_a5_cb <- scan("/work/psychiatry/rgree2/SNRNASEQ/03_CELLBENDER/CS_A5/CellBender_Out_cell_barcodes.csv", what = "", sep = "\n") # 12757
cs_a6_cb <- scan("/work/psychiatry/rgree2/SNRNASEQ/03_CELLBENDER/CS_A6/CellBender_Out_cell_barcodes.csv", what = "", sep = "\n") # 10441
cs_x2_cb <- scan("/work/psychiatry/rgree2/SNRNASEQ/03_CELLBENDER/CS_X2/CellBender_Out_cell_barcodes.csv", what = "", sep = "\n") # 20000
cs_x3_cb <- scan("/work/psychiatry/rgree2/SNRNASEQ/03_CELLBENDER/CS_X3/CellBender_Out_cell_barcodes.csv", what = "", sep = "\n") # 20000


## load dataset from 10X run
cs_a1_cr <- Read10X(data.dir = "/work/psychiatry/rgree2/SNRNASEQ/02_CELLRANGER/CS_A1/A1_filtered_feature_bc_matrix/")
cs_a2_cr <- Read10X(data.dir = "/work/psychiatry/rgree2/SNRNASEQ/02_CELLRANGER/CS_A2/A2_filtered_feature_bc_matrix/")
cs_a5_cr <- Read10X(data.dir = "/work/psychiatry/rgree2/SNRNASEQ/02_CELLRANGER/CS_A5/A5_filtered_feature_bc_matrix/")
cs_a6_cr <- Read10X(data.dir = "/work/psychiatry/rgree2/SNRNASEQ/02_CELLRANGER/CS_A6/A6_filtered_feature_bc_matrix/")
cs_x2_cr <- Read10X(data.dir = "/work/psychiatry/rgree2/SNRNASEQ/02_CELLRANGER/CS_X2/X2_filtered_feature_bc_matrix/")
cs_x3_cr <- Read10X(data.dir = "/work/psychiatry/rgree2/SNRNASEQ/02_CELLRANGER/CS_X3/X3_filtered_feature_bc_matrix/")


## common cell-barcodes
cs_a1.common <- intersect(colnames(cs_a1_cr), cs_a1_cb) # 10110
cs_a2.common <- intersect(colnames(cs_a2_cr), cs_a2_cb) # 11153
cs_a5.common <- intersect(colnames(cs_a5_cr), cs_a5_cb) # 12641
cs_a6.common <- intersect(colnames(cs_a6_cr), cs_a6_cb) # 10423
cs_x2.common <- intersect(colnames(cs_x2_cr), cs_x2_cb) # 20000
cs_x3.common <- intersect(colnames(cs_x3_cr), cs_x3_cb) # 16683


## filtered count tables
cs_a1 <- cs_a1_cr[,colnames(cs_a1_cr) %in% cs_a1.common]    # 53715 10110
cs_a2 <- cs_a2_cr[,colnames(cs_a2_cr) %in% cs_a2.common]    # 53715 11153
cs_a5 <- cs_a5_cr[,colnames(cs_a5_cr) %in% cs_a5.common]    # 53715 12641
cs_a6 <- cs_a6_cr[,colnames(cs_a6_cr) %in% cs_a6.common]    # 53715 10423
cs_x2 <- cs_x2_cr[,colnames(cs_x2_cr) %in% cs_x2.common]    # 53715 20000
cs_x3 <- cs_x3_cr[,colnames(cs_x3_cr) %in% cs_x3.common]    # 53715 16683


## adding sample & condition info to column names
colnames(cs_a1) <- paste("CS_A1", colnames(cs_a1), sep = "_")
colnames(cs_a2) <- paste("CS_A2", colnames(cs_a2), sep = "_")
colnames(cs_a5) <- paste("CS_A5", colnames(cs_a5), sep = "_")
colnames(cs_a6) <- paste("CS_A6", colnames(cs_a6), sep = "_")
colnames(cs_x2) <- paste("CS_X2", colnames(cs_x2), sep = "_")
colnames(cs_x3) <- paste("CS_X3", colnames(cs_x3), sep = "_")


## fetch genes or rows corresponding to gencode protein coding gene symbols
cs_a1.temp <- cs_a1[row.names(cs_a1) %in% ref$GeneSymbol,]
cs_a2.temp <- cs_a2[row.names(cs_a2) %in% ref$GeneSymbol,]
cs_a5.temp <- cs_a5[row.names(cs_a5) %in% ref$GeneSymbol,]
cs_a6.temp <- cs_a6[row.names(cs_a6) %in% ref$GeneSymbol,]
cs_x2.temp <- cs_x2[row.names(cs_x2) %in% ref$GeneSymbol,]
cs_x3.temp <- cs_x3[row.names(cs_x3) %in% ref$GeneSymbol,]


## make a data from from matrix
CS_A1_COUNTS <- as.data.frame(as.matrix(cs_a1.temp))
CS_A2_COUNTS <- as.data.frame(as.matrix(cs_a2.temp))
CS_A5_COUNTS <- as.data.frame(as.matrix(cs_a5.temp))
CS_A6_COUNTS <- as.data.frame(as.matrix(cs_a6.temp))
CS_X2_COUNTS <- as.data.frame(as.matrix(cs_x2.temp))
CS_X3_COUNTS <- as.data.frame(as.matrix(cs_x3.temp))


## add gene symbol as a column
CS_A1_COUNTS$Genes <- row.names(CS_A1_COUNTS)
CS_A2_COUNTS$Genes <- row.names(CS_A2_COUNTS)
CS_A5_COUNTS$Genes <- row.names(CS_A5_COUNTS)
CS_A6_COUNTS$Genes <- row.names(CS_A6_COUNTS)
CS_X2_COUNTS$Genes <- row.names(CS_X2_COUNTS)
CS_X3_COUNTS$Genes <- row.names(CS_X3_COUNTS)


## combine individual tables into a giant data frame
dataCombined <- list("CS_A1" = CS_A1_COUNTS, "CS_A2" = CS_A2_COUNTS, "CS_A5" = CS_A5_COUNTS, "CS_A6" = CS_A6_COUNTS, "CS_X2" = CS_X2_COUNTS, "CS_X3" = CS_X3_COUNTS)

combinedData <- Reduce(function(x, y) { merge(x, y, all = TRUE, by = "Genes") } , dataCombined)
row.names(combinedData) <- combinedData$Genes
combinedData$Genes <- NULL

cs.data <- combinedData[row.names(combinedData) %in% ref$GeneSymbol,]

cs.data[is.na(cs.data)] <- 0
print(dim(cs.data)) ## 21967 81010

save(cs.data, file = "SLEEPSEQ_COUNTS_CS.RData")



##------------------------------------
## SEURAT ANALYSIS FOR CS
rm(list = ls())

## load combined count table
load("SLEEPSEQ_COUNTS_CS.RData")
# cs.data

## Initialize the Seurat object with the raw (non-normalized data).
seuObj <- CreateSeuratObject(counts = cs.data, project = "SLEEPSEQ_CS")
print(dim(seuObj)) ## 21967 81010

## Add and update meta data
metaTemp <- as.data.frame(seuObj@meta.data)
metaTemp2 <- as.data.frame(matrix(unlist(strsplit(row.names(metaTemp), "_")), ncol = 3, byrow = TRUE))
row.names(metaTemp2) <- row.names(metaTemp)
metaData <- merge(metaTemp, metaTemp2, by = "row.names")
row.names(metaData) <- metaData$Row.names
metaData$Row.names <- NULL
colnames(metaData) <- c("Ident", "nUMI", "nGenes", "Condition", "Sample", "CellBarcode")

metaSample <- metaData$Sample
names(metaSample) <- row.names(metaData)

metaCond <- metaData$Condition
names(metaCond) <- row.names(metaData)

seuObj$Condition <- seuObj$orig.ident
seuObj$Sample <- metaSample
seuObj$ConditionSample <- paste(seuObj$Condition, seuObj$Sample, sep = "_")

## Calculate Percent Mito
seuObj[["pMito_RNA"]] <- PercentageFeatureSet(seuObj, pattern = "^mt-")

## Set default identities to GenoAgeSample
Idents(seuObj) <- "ConditionSample"

table(seuObj@active.ident)
# CS_A1 CS_A2 CS_A5 CS_A6 CS_X2 CS_X3 
# 10110 11153 12641 10423 20000 16683

## Visualize Data QC
p1 <- VlnPlot(seuObj, features = c("nFeature_RNA", "nCount_RNA", "pMito_RNA"), ncol = 3, pt.size = 0)
p2 <- FeatureScatter(seuObj, feature1 = "nCount_RNA", feature2 = "pMito_RNA")
p3 <- FeatureScatter(seuObj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

ggsave(filename = "SEURAT_SLEEPSEQ_QC_1.pdf", plot = p1, width = 8, height = 4, units = "in", dpi = 150)
ggsave(filename = "SEURAT_SLEEPSEQ_QC_2.pdf", plot = p2, width = 5, height = 4, units = "in", dpi = 150)
ggsave(filename = "SEURAT_SLEEPSEQ_QC_3.pdf", plot = p3, width = 5, height = 4, units = "in", dpi = 150)

## Data Normalization
seuObj <- NormalizeData(seuObj)

## Identify top 2000 highly variable genes
seuObj <- FindVariableFeatures(seuObj, selection.method = "vst", nfeatures = 2000)

## Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seuObj), 10)

## plot variable features with and without labels
p4 <- VariableFeaturePlot(seuObj)
p5 <- LabelPoints(plot = p4, points = top10, repel = TRUE)
ggsave(filename = "SEURAT_SLEEPSEQ_VARGENES_1.pdf", plot = p4, width = 8, height = 4, units = "in", dpi = 150)
ggsave(filename = "SEURAT_SLEEPSEQ_VARGENES_2.pdf", plot = p5, width = 8, height = 4, units = "in", dpi = 150)

## Data Scaling
all.genes <- rownames(seuObj)
seuObj <- ScaleData(seuObj, features = all.genes, vars.to.regress = c("nCount_RNA", "pMito_RNA"))

## PCA & Jackstraw
## Compute PCA
seuObj <- RunPCA(seuObj, features = VariableFeatures(object = seuObj), npcs = 100)

## Compute and Score Jackstraw
seuObj <- JackStraw(seuObj, num.replicate = 100, dims = 100)
seuObj <- ScoreJackStraw(seuObj, dims = 1:100)

## Plot PCA loadings and Jackstraw scores
p6 <- ElbowPlot(seuObj, ndims = 100)
p7 <- JackStrawPlot(seuObj, dims = 1:100)
ggsave(filename = "SEURAT_SLEEPSEQ_PCA_1.pdf", plot = p6, width = 6, height = 4, units = "in", dpi = 150)
ggsave(filename = "SEURAT_SLEEPSEQ_PCA_2.pdf", plot = p7, width = 12, height = 6, units = "in", dpi = 150)

## Data Clustering
## IDENTIFY PCS
pcScores <- seuObj@reductions$pca@jackstraw$overall.p.values
pcScoresNoSign <- pcScores[pcScores[,2] > 0.05,]
selpcs <- min(pcScoresNoSign[,1]) - 1
print(selpcs)
# 73

seuObj <- FindNeighbors(seuObj, dims = 1:selpcs)
seuObj <- FindClusters(seuObj, resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 0.8, 1.4, 1.6, 1.8, 2.0))
seuObj <- RunUMAP(seuObj, dims = 1:selpcs)

## clustree
seutree <- clustree(seuObj, prefix = "RNA_snn_res.", node_colour = "sc3_stability") # + scale_color_brewer(palette = "Set1") + scale_edge_color_continuous(low = "blue", high = "red")
ggsave(filename = "SLEEPSEQ_SEURAT_CS_CLUSTREE.pdf", plot = seutree, width = 24, height = 24, units = "in", dpi = 150)

## UMAP plots and nUMI violin plots
## for all resolutions
ResolutionList <- grep("RNA_snn_res", colnames(seuObj@meta.data), value = TRUE)

for (Resolution in ResolutionList)
    {
    print(paste("====> RESOLUTION ", Resolution, sep = ""))

    pdf(paste0("SLEEPSEQ_SEURAT_CS_UMAP_RES_", Resolution, ".pdf"), width = 7, height = 6)
    g <- DimPlot(object = seuObj, label = TRUE, reduction = "umap", group.by = Resolution)
    print(g)
    dev.off()

    pdf(paste0("SLEEPSEQ_SEURAT_CS_UMAP_VIOLIN_nUMI_RES_", Resolution, ".pdf"), width = 7, height = 3)
    v <- VlnPlot(object = seuObj, features = "nCount_RNA", ncol = 1, pt.size = 0, group.by = Resolution)
    print(v)
    dev.off()
    }

## Save RData
save(seuObj, file = "SLEEPSEQ_SEURAT_CS.RData")


##------------------------------------
## DOUBLET FINDER FOR CS
## load libraries
rm(list = ls())
library(dplyr)
library(patchwork)
library(ggplot2)
library(reshape2)
library(clustree)
library(Seurat)
library(DoubletFinder)
set.seed(10)

## load seurat object
load("SLEEPSEQ_SEURAT_CS.RData")
# seuObj

pcScores <- seuObj@reductions$pca@jackstraw$overall.p.values
pcScoresNoSign <- pcScores[pcScores[,2] > 0.05,]
selpcs <- min(pcScoresNoSign[,1]) - 1
print(selpcs)
# 73

## pK Identification
sweep.res.list_str <- paramSweep_v3(seuObj, PCs = 1:selpcs, sct = FALSE)
sweep.stats_str <- summarizeSweep(sweep.res.list_str, GT = FALSE)
bcmvn_str <- find.pK(sweep.stats_str)
selpk <- 0.005

## Homotypic Doublet Proportion Estimate
homotypic.prop <- modelHomotypic(seuObj@meta.data$RNA_snn_res.0.8) 
nExp_poi <- round(0.075*nrow(seuObj@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies
seuObjDF <- doubletFinder_v3(seuObj, PCs = 1:selpcs, pN = 0.25, pK = selpk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
seuObjDF <- doubletFinder_v3(seuObjDF, PCs = 1:selpcs, pN = 0.25, pK = selpk, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_6076", sct = FALSE)

table(seuObjDF$DF.classifications_0.25_0.005_5770)
# Doublet Singlet 
#    5770   75240

## Save RData
save(seuObjDF, file = "SLEEPSEQ_SEURAT_CS_DOUBLETFINDER.RData")


```

