## Clustering, Annotation and DEG Analysis

### Seurat analysis for CS
```
# module purge && module load shared slurm python/3.7.x-anaconda
# module load gcc/8.3.0
# module load hdf5_18/1.8.17
# module load R/4.0.2-gccmkl

##------------------------------------
## Libraries
rm(list = ls())
library(dplyr)
library(patchwork)
library(ggplot2)
library(reshape2)
library(clustree)
# library(Matrix)
library(Seurat)
# library(SeuratObject)
set.seed(10)


##------------------------------------
## reference gene symbol list
ref <- read.table("/work/psychiatry/rgree2/SNRNASEQ/04_SEURAT_DOUBLETFINDER/gencode.vM17.protein_coding_gene_id_symbol.txt.gz", header = TRUE, sep = "\t")

##------------------------------------
## load DoubletFinder retained count table
load("/work/psychiatry/rgree2/SNRNASEQ/04_SEURAT_DOUBLETFINDER/01_CS/SLEEPSEQ_SEURAT_CS_DOUBLETFINDER.RData")
# seuObjDF

seuObjDFremoved <- subset(seuObjDF, subset = DF.classifications_0.25_0.005_5770 == "Singlet")

cs.data <- GetAssayData(seuObjDFremoved, slot = "counts") # 21967 75240

save(cs.data, file = "SLEEPSEQ_COUNTS_CS_FILTERED.RData")


##------------------------------------
rm(list = ls())

## load combined count table
load("SLEEPSEQ_COUNTS_CS_FILTERED.RData")
# cs.data

## Initialize the Seurat object with the raw (non-normalized data).
seuObj <- CreateSeuratObject(counts = cs.data, project = "SLEEPSEQ_CS")
print(dim(seuObj)) ## 21967 75240


##------------------------------------
## Add and update meta data
metaTemp <- as.data.frame(seuObj@meta.data)
metaTemp2 <- as.data.frame(matrix(unlist(strsplit(row.names(metaTemp), "_")), ncol = 3, byrow = TRUE))
row.names(metaTemp2) <- row.names(metaTemp)
metaData <- merge(metaTemp, metaTemp2, by = "row.names")
row.names(metaData) <- metaData$Row.names
metaData$Row.names <- NULL
colnames(metaData) <- c("Ident", "nUMI", "nGenes", "Condition", "Sample", "CellBarcode")

metaData$Batch <- metaData$Sample
metaData$Batch <- gsub("A1", "1", metaData$Batch)
metaData$Batch <- gsub("A2", "1", metaData$Batch)
metaData$Batch <- gsub("A5", "1", metaData$Batch)
metaData$Batch <- gsub("A6", "1", metaData$Batch)
metaData$Batch <- gsub("X2", "2", metaData$Batch)
metaData$Batch <- gsub("X3", "2", metaData$Batch)


metaSample <- metaData$Sample
names(metaSample) <- row.names(metaData)

metaCond <- metaData$Condition
names(metaCond) <- row.names(metaData)

metaBatch <- metaData$Batch
names(metaBatch) <- row.names(metaData)


seuObj$Condition <- seuObj$orig.ident
seuObj$Sample <- metaSample
seuObj$ConditionSample <- paste(seuObj$Condition, seuObj$Sample, sep = "_")
seuObj$Batch <- metaBatch

##------------------------------------
## Calculate Percent Mito
seuObj[["pMito_RNA"]] <- PercentageFeatureSet(seuObj, pattern = "^mt-")

## Set default identities to GenoAgeSample
Idents(seuObj) <- "ConditionSample"

table(seuObj@active.ident)
# CS_A1 CS_A2 CS_A5 CS_A6 CS_X2 CS_X3 
#  9462  9774 12334 10084 18460 15126

seuObj <- subset(seuObj, idents = c("CS_X2", "CS_X3"), invert = TRUE)

table(seuObj@active.ident)
# CS_A1 CS_A2 CS_A5 CS_A6 
#  9462  9774 12334 10084


## Visualize Data QC
p1 <- VlnPlot(seuObj, features = c("nFeature_RNA", "nCount_RNA", "pMito_RNA"), ncol = 3, pt.size = 0)
p2 <- FeatureScatter(seuObj, feature1 = "nCount_RNA", feature2 = "pMito_RNA")
p3 <- FeatureScatter(seuObj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

ggsave(filename = "SEURAT_SLEEPSEQ_QC_1.pdf", plot = p1, width = 8, height = 4, units = "in", dpi = 150)
ggsave(filename = "SEURAT_SLEEPSEQ_QC_2.pdf", plot = p2, width = 5, height = 4, units = "in", dpi = 150)
ggsave(filename = "SEURAT_SLEEPSEQ_QC_3.pdf", plot = p3, width = 5, height = 4, units = "in", dpi = 150)


##------------------------------------
## Data Filtering
seuObj <- subset(seuObj, subset = nCount_RNA < 20000 & pMito_RNA < 0.5)
print(dim(seuObj)) ## 21967  41398

##------------------------------------
## Data Normalization
seuObj <- NormalizeData(seuObj)

##------------------------------------
## Identify top 2000 highly variable genes
seuObj <- FindVariableFeatures(seuObj, selection.method = "vst", nfeatures = 2000)

## Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seuObj), 10)

## plot variable features with and without labels
p4 <- VariableFeaturePlot(seuObj)
p5 <- LabelPoints(plot = p4, points = top10, repel = TRUE)

ggsave(filename = "SEURAT_SLEEPSEQ_VARGENES_1.pdf", plot = p4, width = 8, height = 4, units = "in", dpi = 150)
ggsave(filename = "SEURAT_SLEEPSEQ_VARGENES_2.pdf", plot = p5, width = 8, height = 4, units = "in", dpi = 150)


##------------------------------------
## Data Scaling
all.genes <- rownames(seuObj)
seuObj <- ScaleData(seuObj, features = all.genes, vars.to.regress = c("nCount_RNA", "pMito_RNA"))

##------------------------------------
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


# save(seuObj, file = "SLEEPSEQ_SEURAT_CS_NORM_SCALE_PCA.RData")


##------------------------------------
# rm(list = ls())

# load("SLEEPSEQ_SEURAT_CS_NORM_SCALE_PCA.RData")
# # seuObj

## Data Clustering
## IDENTIFY PCS
pcScores <- seuObj@reductions$pca@jackstraw$overall.p.values
pcScoresNoSign <- pcScores[pcScores[,2] > 0.05,]
selpcs <- min(pcScoresNoSign[,1]) - 1
print(selpcs)


seuObj <- FindNeighbors(seuObj, dims = 1:selpcs)
seuObj <- FindClusters(seuObj, resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 0.8, 1.4, 1.6, 1.8, 2.0), n.iter = 1000)
seuObj <- RunUMAP(seuObj, dims = 1:selpcs, n.epochs = 200, min.dist = 0.001, spread = 5)


##------------------------------------
## clustree
seutree <- clustree(seuObj, prefix = "RNA_snn_res.", node_colour = "sc3_stability") # + scale_color_brewer(palette = "Set1") + scale_edge_color_continuous(low = "blue", high = "red")
ggsave(filename = "SLEEPSEQ_SEURAT_CS_CLUSTREE.pdf", plot = seutree, width = 24, height = 24, units = "in", dpi = 150)


##------------------------------------
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


##------------------------------------
## UMAP/Feature plots for marker genes
# mygenes <- c("Mki67", "Sox4", "Sox11", "Foxp1", "Foxp2", "Drd", "Drd2", "Tac1", "Penk", "Casz1", "Aqp4", "Olig1", "Cx3cr1", "Flt1", "Slc17a7", "Dlx2", "Npy", "Ascl1", "Sp9", "Ppp1r1b")

mygenes <- c("Aqp4", "Aldoc", "Tac1", "Clock", "Thy1", "Slc17a6", "Slc17a7", "Gad1", "Foxp1", "Foxp2", "Olig1", "Olig2", "Mobp", "Mbp", "Mag", "Cldn5", "Cx3cr1", "Vtn", "Mrc1")

fpl1 <- FeaturePlot(object = seuObj, features = mygenes, cols = c("gray90", "blue"), reduction = "umap", ncol = 5, pt.size = 0.1, min.cutoff = 0, max.cutoff = 3, order = TRUE, raster = TRUE)
ggsave(filename = "SLEEPSEQ_SEURAT_PLOT_FeaturePlot_orderT.pdf", plot = fpl1, width = 36, height = 24, units = "in", dpi = 300)

fpl2 <- FeaturePlot(object = seuObj, features = mygenes, cols = c("gray90", "blue"), reduction = "umap", ncol = 5, pt.size = 0.1, min.cutoff = 0, max.cutoff = 3, order = FALSE, raster = TRUE)
ggsave(filename = "SLEEPSEQ_SEURAT_PLOT_FeaturePlot_orderF.pdf", plot = fpl2, width = 36, height = 24, units = "in", dpi = 300)

umap1 <- DimPlot(object = seuObj, group.by = "RNA_snn_res.1.4", pt.size = 0.1, label = TRUE, reduction = "umap", raster = TRUE) + NoLegend()
ggsave(filename = "SLEEPSEQ_SEURAT_PLOT_DimPlot_Res_1.4.pdf", plot = umap1, width = 6, height = 6, units = "in", dpi = 300)


##------------------------------------
## Save RData
save(seuObj, file = "SLEEPSEQ_SEURAT_CS.RData")




# ##------------------------------------
# ## Find Cluster Markers
# rm(list = ls())
# load("SLEEPSEQ_SEURAT_CS.RData")
# # seuObj

# Idents(seuObj) <- "RNA_snn_res.1.4"

# seuObj.markers <- FindAllMarkers(seuObj, only.pos = TRUE, min.pct = 0, logfc.threshold = 0)
# write.table(seuObj.markers, "SLEEPSEQ_SEURAT_RES_1.4_CLUSTER_MARKERS.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")


```

### Seurat analysis for SD
```
# module purge && module load shared slurm python/3.7.x-anaconda
# module load gcc/8.3.0
# module load hdf5_18/1.8.17
# module load R/4.0.2-gccmkl

##------------------------------------
## Libraries
rm(list = ls())
library(dplyr)
library(patchwork)
library(ggplot2)
library(reshape2)
library(clustree)
# library(Matrix)
library(Seurat)
# library(SeuratObject)
set.seed(10)

##------------------------------------
## reference gene symbol list
ref <- read.table("/work/psychiatry/rgree2/SNRNASEQ/04_SEURAT_DOUBLETFINDER/gencode.vM17.protein_coding_gene_id_symbol.txt.gz", header = TRUE, sep = "\t")

##------------------------------------
## load DoubletFinder retained count table
load("/work/psychiatry/rgree2/SNRNASEQ/04_SEURAT_DOUBLETFINDER/02_SD/SLEEPSEQ_SEURAT_SD_DOUBLETFINDER.RData")
# seuObjDF

seuObjDFremoved <- subset(seuObjDF, subset = DF.classifications_0.25_0.005_6209 == "Singlet")

sd.data <- GetAssayData(seuObjDFremoved, slot = "counts") # 21967 80353

save(sd.data, file = "SLEEPSEQ_COUNTS_SD_FILTERED.RData")


##------------------------------------
rm(list = ls())

## load combined count table
load("SLEEPSEQ_COUNTS_SD_FILTERED.RData")
# sd.data

## Initialize the Seurat object with the raw (non-normalized data).
seuObj <- CreateSeuratObject(counts = sd.data, project = "SLEEPSEQ_SD")
print(dim(seuObj)) ## 21967 80353


##------------------------------------
## Add and update meta data
metaTemp <- as.data.frame(seuObj@meta.data)
metaTemp2 <- as.data.frame(matrix(unlist(strsplit(row.names(metaTemp), "_")), ncol = 3, byrow = TRUE))
row.names(metaTemp2) <- row.names(metaTemp)
metaData <- merge(metaTemp, metaTemp2, by = "row.names")
row.names(metaData) <- metaData$Row.names
metaData$Row.names <- NULL
colnames(metaData) <- c("Ident", "nUMI", "nGenes", "Condition", "Sample", "CellBarcode")

metaData$Batch <- metaData$Sample
metaData$Batch <- gsub("A3", "1", metaData$Batch)
metaData$Batch <- gsub("A4", "1", metaData$Batch)
metaData$Batch <- gsub("A7", "1", metaData$Batch)
metaData$Batch <- gsub("A8", "1", metaData$Batch)
metaData$Batch <- gsub("X5", "2", metaData$Batch)
metaData$Batch <- gsub("X7", "2", metaData$Batch)


metaSample <- metaData$Sample
names(metaSample) <- row.names(metaData)

metaCond <- metaData$Condition
names(metaCond) <- row.names(metaData)

metaBatch <- metaData$Batch
names(metaBatch) <- row.names(metaData)


seuObj$Condition <- seuObj$orig.ident
seuObj$Sample <- metaSample
seuObj$ConditionSample <- paste(seuObj$Condition, seuObj$Sample, sep = "_")
seuObj$Batch <- metaBatch

##------------------------------------
## Calculate Percent Mito
seuObj[["pMito_RNA"]] <- PercentageFeatureSet(seuObj, pattern = "^mt-")

## Set default identities to GenoAgeSample
Idents(seuObj) <- "ConditionSample"

table(seuObj@active.ident)
# SD_A3 SD_A4 SD_A7 SD_A8 SD_X5 SD_X7 
#  9868 10926 12412 10066 17982 19099

seuObj <- subset(seuObj, idents = c("SD_X5", "SD_X7"), invert = TRUE)

table(seuObj@active.ident)
# SD_A3 SD_A4 SD_A7 SD_A8
#  9868 10926 12412 10066

## Visualize Data QC
p1 <- VlnPlot(seuObj, features = c("nFeature_RNA", "nCount_RNA", "pMito_RNA"), ncol = 3, pt.size = 0)
p2 <- FeatureScatter(seuObj, feature1 = "nCount_RNA", feature2 = "pMito_RNA")
p3 <- FeatureScatter(seuObj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

ggsave(filename = "SEURAT_SLEEPSEQ_QC_1.pdf", plot = p1, width = 8, height = 4, units = "in", dpi = 150)
ggsave(filename = "SEURAT_SLEEPSEQ_QC_2.pdf", plot = p2, width = 5, height = 4, units = "in", dpi = 150)
ggsave(filename = "SEURAT_SLEEPSEQ_QC_3.pdf", plot = p3, width = 5, height = 4, units = "in", dpi = 150)


##------------------------------------
## Data Filtering
seuObj <- subset(seuObj, subset = nCount_RNA < 20000 & pMito_RNA < 0.5)
print(dim(seuObj)) ## 21967  80112

##------------------------------------
## Data Normalization
seuObj <- NormalizeData(seuObj)

##------------------------------------
## Identify top 2000 highly variable genes
seuObj <- FindVariableFeatures(seuObj, selection.method = "vst", nfeatures = 2000)

## Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seuObj), 10)

## plot variable features with and without labels
p4 <- VariableFeaturePlot(seuObj)
p5 <- LabelPoints(plot = p4, points = top10, repel = TRUE)

ggsave(filename = "SEURAT_SLEEPSEQ_VARGENES_1.pdf", plot = p4, width = 8, height = 4, units = "in", dpi = 150)
ggsave(filename = "SEURAT_SLEEPSEQ_VARGENES_2.pdf", plot = p5, width = 8, height = 4, units = "in", dpi = 150)


##------------------------------------
## Data Scaling
all.genes <- rownames(seuObj)
seuObj <- ScaleData(seuObj, features = all.genes, vars.to.regress = c("nCount_RNA", "pMito_RNA"))

##------------------------------------
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


# save(seuObj, file = "SLEEPSEQ_SEURAT_CS_NORM_SCALE_PCA.RData")


##------------------------------------
# rm(list = ls())

# load("SLEEPSEQ_SEURAT_CS_NORM_SCALE_PCA.RData")
# # seuObj

## Data Clustering
## IDENTIFY PCS
pcScores <- seuObj@reductions$pca@jackstraw$overall.p.values
pcScoresNoSign <- pcScores[pcScores[,2] > 0.05,]
selpcs <- min(pcScoresNoSign[,1]) - 1
print(selpcs)


seuObj <- FindNeighbors(seuObj, dims = 1:selpcs)
seuObj <- FindClusters(seuObj, resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 0.8, 1.4, 1.6, 1.8, 2.0), n.iter = 1000)
seuObj <- RunUMAP(seuObj, dims = 1:selpcs, n.epochs = 200, min.dist = 0.001, spread = 5)


##------------------------------------
## clustree
seutree <- clustree(seuObj, prefix = "RNA_snn_res.", node_colour = "sc3_stability") # + scale_color_brewer(palette = "Set1") + scale_edge_color_continuous(low = "blue", high = "red")
ggsave(filename = "SLEEPSEQ_SEURAT_SD_CLUSTREE.pdf", plot = seutree, width = 24, height = 24, units = "in", dpi = 150)


##------------------------------------
## UMAP plots and nUMI violin plots
## for all resolutions
ResolutionList <- grep("RNA_snn_res", colnames(seuObj@meta.data), value = TRUE)

for (Resolution in ResolutionList)
    {
    print(paste("====> RESOLUTION ", Resolution, sep = ""))

    pdf(paste0("SLEEPSEQ_SEURAT_SD_UMAP_RES_", Resolution, ".pdf"), width = 7, height = 6)
    g <- DimPlot(object = seuObj, label = TRUE, reduction = "umap", group.by = Resolution)
    print(g)
    dev.off()

    pdf(paste0("SLEEPSEQ_SEURAT_SD_UMAP_VIOLIN_nUMI_RES_", Resolution, ".pdf"), width = 7, height = 3)
    v <- VlnPlot(object = seuObj, features = "nCount_RNA", ncol = 1, pt.size = 0, group.by = Resolution)
    print(v)
    dev.off()
    }



##------------------------------------
## UMAP/Feature plots for marker genes
# mygenes <- c("Mki67", "Sox4", "Sox11", "Foxp1", "Foxp2", "Drd", "Drd2", "Tac1", "Penk", "Casz1", "Aqp4", "Olig1", "Cx3cr1", "Flt1", "Slc17a7", "Dlx2", "Npy", "Ascl1", "Sp9", "Ppp1r1b")

mygenes <- c("Aqp4", "Aldoc", "Tac1", "Clock", "Thy1", "Slc17a6", "Slc17a7", "Gad1", "Foxp1", "Foxp2", "Olig1", "Olig2", "Mobp", "Mbp", "Mag", "Cldn5", "Cx3cr1", "Vtn", "Mrc1")

fpl1 <- FeaturePlot(object = seuObj, features = mygenes, cols = c("gray90", "blue"), reduction = "umap", ncol = 5, pt.size = 0.1, min.cutoff = 0, max.cutoff = 3, order = TRUE, raster = TRUE)
ggsave(filename = "SLEEPSEQ_SEURAT_PLOT_FeaturePlot_orderT.pdf", plot = fpl1, width = 36, height = 24, units = "in", dpi = 300)

fpl2 <- FeaturePlot(object = seuObj, features = mygenes, cols = c("gray90", "blue"), reduction = "umap", ncol = 5, pt.size = 0.1, min.cutoff = 0, max.cutoff = 3, order = FALSE, raster = TRUE)
ggsave(filename = "SLEEPSEQ_SEURAT_PLOT_FeaturePlot_orderF.pdf", plot = fpl2, width = 36, height = 24, units = "in", dpi = 300)

umap1 <- DimPlot(object = seuObj, group.by = "RNA_snn_res.1.4", pt.size = 0.1, label = TRUE, reduction = "umap", raster = TRUE) + NoLegend()
ggsave(filename = "SLEEPSEQ_SEURAT_PLOT_DimPlot_Res_1.4.pdf", plot = umap1, width = 6, height = 6, units = "in", dpi = 300)


##------------------------------------
## Save RData
save(seuObj, file = "SLEEPSEQ_SEURAT_SD.RData")


```

### Integration
```

rm(list = ls())
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(ggrepel)
library(reticulate)
library(clustree)


##-------------------------------------------------------
## SEURAT ANALYSIS | LOAD INDIVIDUAL DATA
##-------------------------------------------------------

print("SLEEPSEQ_CS")
load("/work/psychiatry/rgree2/SNRNASEQ/06_SEURAT_woAmbient_woDoublets/01_CS/SLEEPSEQ_SEURAT_CS.RData")
sleepseq.cs <- seuObj
rm(seuObj)
DefaultAssay(sleepseq.cs) <- "RNA"
print(dim(sleepseq.cs))


print("SLEEPSEQ_SD")
load("/work/psychiatry/rgree2/SNRNASEQ/06_SEURAT_woAmbient_woDoublets/02_SD/SLEEPSEQ_SEURAT_SD.RData")
sleepseq.sd <- seuObj
rm(seuObj)
DefaultAssay(sleepseq.sd) <- "RNA"
print(dim(sleepseq.sd))



##-------------------------------------------------------
## SEURAT ANALYSIS | INTEGRATE DATA
##-------------------------------------------------------
## create a list of seurat objects
sleepseq.list <- c(sleepseq.cs, sleepseq.sd)

## standard preprocessing, identify variable features individually
# normalize and identify variable features for each dataset independently
sleepseq.list <- lapply(X = sleepseq.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
	})


# select features that are repeatedly variable across datasets for integration
sleepseq.features <- SelectIntegrationFeatures(object.list = sleepseq.list)

## identify integration anchors
sleepseq.anchors <- FindIntegrationAnchors(object.list = sleepseq.list, anchor.features = sleepseq.features)

## integrate datasets using anchors
sleepseq.integrated <- IntegrateData(anchorset = sleepseq.anchors)

## scaling, pca and clustering
## switch to integrated assay
DefaultAssay(sleepseq.integrated) <- "integrated"

# save(sleepseq.integrated, sleepseq.anchors, file = "SEURAT_ORGANOIDS_DATA_FILT_NORM_INTEGRATED.RData")


# Run the standard workflow for visualization and clustering
myres <- c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0)
sleepseq.integrated <- ScaleData(sleepseq.integrated, verbose = FALSE)
sleepseq.integrated <- RunPCA(sleepseq.integrated, npcs = 30, verbose = FALSE)
sleepseq.integrated <- RunUMAP(sleepseq.integrated, reduction = "pca", dims = 1:30, n.epochs = 200, min.dist = 0.001, spread = 5)
sleepseq.integrated <- FindNeighbors(object = sleepseq.integrated, dims = 1:30) 
sleepseq.integrated <- FindClusters(object = sleepseq.integrated, resolution = myres, n.iter = 100)


save(sleepseq.integrated, sleepseq.anchors, file = "SLEEPSEQ_DATA_FILT_NORM_PCA_CLU_INTEGRATED.RData")



##-------------------------------------------------------
## SEURAT ANALYSIS | FIGURES
##-------------------------------------------------------

seutree <- clustree(sleepseq.integrated, prefix = "integrated_snn_res.", node_colour = "sc3_stability") # + scale_color_brewer(palette = "Set1") + scale_edge_color_continuous(low = "blue", high = "red")
ggsave(filename = "SLEEPSEQ_DATA_FILT_NORM_PCA_CLU_INTEGRATED_CLUSTREE.pdf", plot = seutree, width = 24, height = 24, units = "in", dpi = 150)


# ResolutionList <- grep("integrated_snn_res", colnames(sleepseq.integrated@meta.data), value = TRUE)

# for (Resolution in ResolutionList){
#     print(Resolution)
#     pdf(paste0("EP_ORG_DATA_FILT_NORM_PCA_CLU_INTEGRATED_UMAP_RES_", Resolution, ".pdf"), width=7, height=7)
#     g <- DimPlot(object = sleepseq.integrated, label = TRUE, reduction = "umap", group.by = Resolution)
#     print(g)
#     dev.off()
#     }

ResolutionList <- grep("integrated_snn_res", colnames(sleepseq.integrated@meta.data), value = TRUE)

for (Resolution in ResolutionList)
    {
    print(Resolution)
    pdf(paste0("SLEEPSEQ_DATA_FILT_NORM_PCA_CLU_INTEGRATED_UMAP_RES_", Resolution, ".pdf"), width=8, height=7)
    g <- DimPlot(object = sleepseq.integrated, label = TRUE, reduction = "umap", group.by = Resolution, raster = TRUE)
    print(g)
    dev.off()

    pdf(paste0("SLEEPSEQ_DATA_INTEGRATE_VIOLIN_nUMI_RES_", Resolution, ".pdf"), width = 9, height = 3)
    v <- VlnPlot(object = sleepseq.integrated, features = "nCount_RNA", ncol = 1, pt.size = 0, group.by = Resolution)
    print(v)
    dev.off()
    }


plotCluUMAP1 <- DimPlot(object = sleepseq.integrated, group.by = "Condition", pt.size = 0.1, reduction = "umap", label = TRUE, label.size = 3, repel = TRUE, raster = TRUE)
ggsave(filename = "SEURAT_SLEEPSEQ_DATA_UMAP_CONDITION.pdf", plot = plotCluUMAP1, width = 9, height = 8, units = "in", dpi = 150)

plotCluUMAP1b <- DimPlot(object = sleepseq.integrated, group.by = "Condition", split.by = "Condition", pt.size = 0.1, reduction = "umap", label = FALSE, label.size = 3, repel = TRUE, raster = TRUE)
ggsave(filename = "SEURAT_SLEEPSEQ_DATA_UMAP_CONDITION_FACET.pdf", plot = plotCluUMAP1b, width = 17, height = 8, units = "in", dpi = 150)


mygenes <- c("Aqp4", "Aldoc", "Tac1", "Clock", "Thy1", "Slc17a6", "Slc17a7", "Gad1", "Foxp1", "Foxp2", "Olig1", "Olig2", "Mobp", "Mbp", "Mag", "Cldn5", "Cx3cr1", "Vtn", "Mrc1")

DefaultAssay(sleepseq.integrated) <- "RNA"

fpl1 <- FeaturePlot(object = sleepseq.integrated, features = mygenes, cols = c("gray90", "blue"), reduction = "umap", ncol = 3, pt.size = 0.05, min.cutoff = 0, max.cutoff = 3, order = TRUE, raster = TRUE)
ggsave(filename = "SEURAT_SLEEPSEQ_DATA_FeaturePlot_orderT.pdf", plot = fpl1, width = 14, height = 28, units = "in", dpi = 150, useDingbats=FALSE, )

fpl2 <- FeaturePlot(object = sleepseq.integrated, features = mygenes, cols = c("gray90", "blue"), reduction = "umap", ncol = 3, pt.size = 0.05, min.cutoff = 0, max.cutoff = 3, order = FALSE, raster = TRUE)
ggsave(filename = "SEURAT_SLEEPSEQ_DATA_FeaturePlot_orderF.pdf", plot = fpl2, width = 14, height = 28, units = "in", dpi = 150, useDingbats=FALSE)




##-------------------------------------------------------
## SEURAT ANALYSIS | CLUSTER MARKERS
##-------------------------------------------------------
rm(list = ls())
load("SLEEPSEQ_DATA_FILT_NORM_PCA_CLU_INTEGRATED.RData")
# sleepseq.integrated

DefaultAssay(object = sleepseq.integrated) <- "integrated"

Idents(sleepseq.integrated) <- "integrated_snn_res.1.2"

seuMarkers <- FindAllMarkers(object = sleepseq.integrated)

write.table(seuMarkers, "SLEEPSEQ_DATA_FILT_NORM_PCA_CLU_INTEGRATED_DEG.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
save(seuMarkers, file = "SLEEPSEQ_DATA_FILT_NORM_PCA_CLU_INTEGRATED_DEG.RData")




##-------------------------------------------------------
## SEURAT ANALYSIS | ANNOTATE CLUSTERS
##-------------------------------------------------------

rm(list = ls())
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(ggrepel)
library(reticulate)
library(WGCNA)
library(dplyr)
library(tidyr)


## Load cluster markers
cluMarkers <- read.table("SLEEPSEQ_DATA_FILT_NORM_PCA_CLU_INTEGRATED_DEG.txt", sep = "\t", header = TRUE)
tab <- cluMarkers
tab <- tab[tab$p_val_adj <= 0.05 & tab$avg_log2FC >= 0.75 & tab$pct.1 >= 0.75,]
# tab <- tab[tab$p_val_adj <= 0.05 & tab$avg_log2FC >= 0.25,]
# tab <- tab[tab$p_val_adj <= 0.05 & tab$pct.1 >= 0.25,]
tab <- tab[c(7,6)]
tab$cluster <- as.factor(paste("Cluster_", sprintf("%02d", as.numeric(as.character(tab$cluster))), sep=""))
colnames(tab)=c("Gene","DEFINITION")
tab$DEFINITION <- as.factor(tab$DEFINITION)
Genes=as.data.frame(table(tab$DEFINITION))


## Load reference data (BICCN)
biccn_master <- read.table("/work/psychiatry/rgree2/RESOURCES/DATASETS/2021_BICCN/2020-03-03645C-s4/Supplementary_Table_6.txt", header = TRUE, sep = ",")
biccn_master_up <- biccn_master[biccn_master$direction == "up",]
biccn_master_up2 <- biccn_master_up[,c("genes", "cl1_cluster_label")]
colnames(biccn_master_up2) <- c("gene", "celltype")
biccn_list <- split(biccn_master_up2, biccn_master_up2$celltype)
biccn_list2 <- lapply(biccn_list, function(x) { x[1:50,] } )

GeneSets <- biccn_list2


for(i in 1:length(GeneSets)){
	colnames(GeneSets[[i]])[1] <- "Gene"
}

ln=length(GeneSets)
cl=length(Genes$Var1)
TEMP=list()
INT=list()
for (i in 1:ln)
{
TEMP[[i]]=tab[tab$Gene %in% GeneSets[[i]]$Gene,]
INT[[i]]=as.data.frame(table(TEMP[[i]]$DEFINITION))
}
names(INT)=names(GeneSets)
names(TEMP)=names(GeneSets)

# Create the matrix for each GeneSet
NROWS <- sapply(GeneSets,nrow)

#
#               GeneSet
#              +       -
#          +-------+-------+
#        + |   a   |   b   |
#  Module  +-------+-------+
#        - |   c   |   d   |
#          +-------+-------+

for (i in 1:length(INT))
{
INT[[i]]$b <- NROWS[[i]]-INT[[i]]$Freq
INT[[i]]$c <- Genes$Freq-INT[[i]]$Freq 
INT[[i]]$d <- 15585-Genes$Freq-nrow(GeneSets[[i]]) #19776 #15585 #13517 #nrow(tab) #8321 # nrow(cluMarkers)(genes in seurat object)
}

# sum(Genes$Freq)
RunFisher <- function(row, alt = 'greater', cnf = 0.85) {
	f <- fisher.test(matrix(row, nrow = 2), alternative = alt, conf.level = cnf)
	return(c(row,
			P_val = f$p.value,
			LogP = -log10(f$p.value), 
			OR = f$estimate[[1]],
			OR_Low = f$conf.int[1],
			OR_Up = f$conf.int[2]))
}

# run
FisherMat=list()
for (i in 1:length(INT))
{
FisherMat[[i]] <- t(apply(INT[[i]][,2:5], 1, RunFisher))
rownames(FisherMat[[i]]) <- INT[[i]]$Var1
FisherMat[[i]] <- FisherMat[[i]][rownames(FisherMat[[i]]) != "grey",]
}
names(FisherMat)<-names(INT)
save(FisherMat, TEMP, file= "SLEEPSEQ_SEURAT_FisherOutput_BICCN_Enrich.RData")


# Create matrices of Pval
tmp<-list()
FisherP<-matrix()
rowNames <- rownames(FisherMat[[1]])
colNames <- names(FisherMat)
for (i in 1:length(INT))
{
tmp[[i]] <- cbind(as.numeric(FisherMat[[i]][,5]))
FisherP <- do.call(cbind,tmp)
}
rownames(FisherP) <- rowNames
colnames(FisherP) <- colNames

# Create matrices of OR
tmp<-list()
FisherOR<-matrix()
rowNames <- rownames(FisherMat[[1]])
colNames <- names(FisherMat)
for (i in 1:length(INT))
{
tmp[[i]] <- cbind(as.numeric(FisherMat[[i]][,6]))
FisherOR <- do.call(cbind,tmp)
}
rownames(FisherOR) <- rowNames
colnames(FisherOR) <- colNames

# Pval Adjusted
library(magrittr)
FisherAdj <- FisherP %>% as.matrix %>% as.vector %>% p.adjust(method='fdr') %>% matrix(ncol=ncol(FisherP))
rownames(FisherAdj) <- rowNames
colnames(FisherAdj) <- colNames

FisherAdj[FisherAdj>0.05]=1
FisherOR[FisherOR < 1]=0


pdf("SLEEPSEQ_SEURAT_FisherOutput_BICCN_Enrich.pdf", width=48, height=24, pointsize=12)
par(mar=c(15, 7, 2, 2))
df=-log10(FisherAdj)
LabelMat = paste(signif(FisherOR, 2), "\n(",signif(FisherAdj, 1), ")", sep = "")
LabelMat[LabelMat == "0\n(1)"] <- NA
labeledHeatmap(Matrix = df, 
xLabels = colnames(df), 
yLabels = rownames(df), 
colorLabels =FALSE,
colors=colorRampPalette(c("white", "red"))(50),
textMatrix=LabelMat, 
setStdMargins = FALSE, 
cex.text = 0.5,
xLabelsAngle = 90)
dev.off()



FisherORt <- as.data.frame(t(FisherOR))
colnames(FisherORt) <- paste(colnames(FisherORt), "OR", sep = "_")

FisherAdjt <- as.data.frame(t(FisherAdj))
colnames(FisherAdjt) <- paste(colnames(FisherAdjt), "Pval", sep = "_")

FisherData <- merge(FisherORt, FisherAdjt, by = "row.names")
row.names(FisherData) <- FisherData$Row.names
FisherData$Row.names <- NULL

FisherData2 <- FisherData[,order(colnames(FisherData))]
write.table(FisherData2, "SLEEPSEQ_SEURAT_FisherOutput_BICCN_Enrich_PlotData.txt", row.names = T, col.names = T, quote = F, sep = "\t")


sink("SLEEPSEQ_SEURAT_BICCN_AnnotationData.txt")
aa <- FisherData2
for(i in 1:ncol(aa))
 {	 
 if(i %% 2 != 0)
  {
  cluor <- i
  clup <- i + 1
  print("------------------")
  bb <- aa[,c(cluor, clup)]
  cc <- bb[order(bb[,1], decreasing = T),]
  dd <- bb[order(bb[,2], decreasing = F),]
  print(gsub("_OR", "", colnames(cc)[1]))
  print(cc[1:5,])
  print(dd[1:5,])
  print("------------------")
  }
 }

sink()



##-------------------------------------------------------
## SEURAT ANALYSIS | PLOTS & DATA
##-------------------------------------------------------
rm(list = ls())
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(ggrepel)
library(reticulate)
library(clustree)
library(scCustomize)

load("SLEEPSEQ_DATA_FILT_NORM_PCA_CLU_INTEGRATED.RData")
# sleepseq.integrated

DefaultAssay(object = sleepseq.integrated) <- "integrated"

Idents(sleepseq.integrated) <- "integrated_snn_res.1.2"

sleepseq.integrated@active.ident <- factor(x = sleepseq.integrated@active.ident, levels = seq(0, max(as.numeric(sleepseq.integrated$integrated_snn_res.1.2) - 1)))


## UMAP PLOTS
p1 <- DimPlot_scCustom(seurat_object = sleepseq.integrated, pt.size = 0.01, group.by = "integrated_snn_res.1.2", raster = TRUE, ggplot_default_colors = TRUE, label = TRUE, label.size = 3) + NoLegend()
ggsave(filename = "SEURAT_SLEEPSEQ_UMAP_CLUSTERS_RES_1.2.pdf", plot = p1, width = 6, height = 6, units = "in", dpi = 300)

p2 <- DimPlot_scCustom(seurat_object = sleepseq.integrated, pt.size = 0.01, group.by = "Condition", ggplot_default_colors = TRUE, raster = TRUE, label = FALSE) #+ NoLegend()
ggsave(filename = "SEURAT_SLEEPSEQ_UMAP_CLUSTERS_CONDITION.pdf", plot = p2, width = 6.5, height = 6, units = "in", dpi = 300)

p3a <- DimPlot_All_Samples(seurat_object = sleepseq.integrated, meta_data_column = "ConditionSample", num_col = 6, pt.size = 0.01, raster = TRUE, label = FALSE, colors_use = DiscretePalette_scCustomize(num_colors = 8, palette = "varibow")) #+ NoLegend()
ggsave(filename = "SEURAT_SLEEPSEQ_UMAP_CLUSTERS_SAMPLE_FACET.pdf", plot = p3a, width = 18, height = 6, units = "in", dpi = 300)

p3b <- DimPlot_scCustom(seurat_object = sleepseq.integrated, pt.size = 0.01, group.by = "ConditionSample", colors_use = DiscretePalette_scCustomize(num_colors = 8, palette = "varibow"), raster = TRUE, label = FALSE) #+ NoLegend()
ggsave(filename = "SEURAT_SLEEPSEQ_UMAP_CLUSTERS_SAMPLE.pdf", plot = p3b, width = 7, height = 6, units = "in", dpi = 300)


## BARPLOT PLOTS
cellsPerCluster <- as.data.frame.matrix(table(sleepseq.integrated$integrated_snn_res.1.2, sleepseq.integrated@meta.data$Condition))
cellsPerCluster$Cluster <- paste("Cluster", sprintf("%02d", as.numeric(row.names(cellsPerCluster))), sep = "_")
write.table(cellsPerCluster, "SEURAT_SLEEPSEQ_CELLS_BY_CLUSTER_BY_CONDITION.tsv", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

cellsPerCluster2 <- melt(cellsPerCluster)
colnames(cellsPerCluster2) <- c("CLUSTER", "CONDITION", "CELLS")
p4 <- ggplot(cellsPerCluster2) +
        geom_bar(aes(x = CLUSTER, y = CELLS, fill = CONDITION), stat = "identity", position = "fill") + 
        scale_y_continuous(labels = scales::percent_format()) + 
        theme_classic() + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
        NULL
ggsave(filename = "SEURAT_SLEEPSEQ_CELLS_BY_CLUSTER_BY_CONDITION.pdf", plot = p4, width = 14, height = 4, units = "in", dpi = 300)


samplePerCluster <- as.data.frame.matrix(table(sleepseq.integrated$integrated_snn_res.1.2, sleepseq.integrated@meta.data$ConditionSample))
samplePerCluster$Cluster <- paste("Cluster", sprintf("%02d", as.numeric(row.names(samplePerCluster))), sep = "_")
write.table(samplePerCluster, "SEURAT_SLEEPSEQ_CELLS_BY_CLUSTER_BY_SAMPLE.tsv", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

samplePerCluster2 <- melt(samplePerCluster)
colnames(samplePerCluster2) <- c("CLUSTER", "SAMPLE", "CELLS")
p5 <- ggplot(samplePerCluster2) +
        geom_bar(aes(x = CLUSTER, y = CELLS, fill = SAMPLE), stat = "identity", position = "fill") + 
        scale_y_continuous(labels = scales::percent_format()) + 
        theme_classic() + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
        NULL
ggsave(filename = "SEURAT_SLEEPSEQ_CELLS_BY_CLUSTER_BY_SAMPLE.pdf", plot = p5, width = 14, height = 4, units = "in", dpi = 300)



## VOILIN PLOTS
p6 <- Stacked_VlnPlot(seurat_object = sleepseq.integrated, features = c("nCount_RNA", "nFeature_RNA"), x_lab_rotate = TRUE, ggplot_default_colors = TRUE)
ggsave(filename = "SEURAT_SLEEPSEQ_VIOLINPLOT_nUMI_nGENES.pdf", plot = p6, width = 14, height = 4, units = "in", dpi = 300)


DefaultAssay(object = sleepseq.integrated) <- "RNA"
mygenes <- c("Aqp4", "Aldoc", "Tac1", "Clock", "Thy1", "Slc17a6", "Slc17a7", "Gad1", "Foxp1", "Foxp2", "Olig1", "Olig2", "Mobp", "Mbp", "Mag", "Cldn5", "Cx3cr1", "Vtn", "Mrc1")
p7 <- Stacked_VlnPlot(seurat_object = sleepseq.integrated, features = mygenes, x_lab_rotate = TRUE, ggplot_default_colors = TRUE)
ggsave(filename = "SEURAT_SLEEPSEQ_VIOLINPLOT_MARKER_GENES.pdf", plot = p7, width = 14, height = 14, units = "in", dpi = 300)



## FEATURE PLOTS
DefaultAssay(object = sleepseq.integrated) <- "RNA"
mygenes <- c("Aqp4", "Aldoc", "Tac1", "Clock", "Thy1", "Slc17a6", "Slc17a7", "Gad1", "Foxp1", "Foxp2", "Olig1", "Olig2", "Mobp", "Mbp", "Mag", "Cldn5", "Cx3cr1", "Vtn", "Mrc1")
p8 <- FeaturePlot_scCustom(seurat_object = sleepseq.integrated, features = mygenes, order = FALSE, pt.size = 0.01, raster = TRUE, num_columns = 5)
ggsave(filename = "SEURAT_SLEEPSEQ_FEATUREPLOT_MARKER_GENES.pdf", plot = p8, width = 21, height = 14, units = "in", dpi = 300)


p9 <- Cluster_Highlight_Plot(seurat_object = sleepseq.integrated, cluster_name = "6", highlight_color = "forestgreen", background_color = "lightgray")
ggsave(filename = "TEMP_UMAP_HIGHLIGHT.pdf", plot = p9, width = 7, height = 6, units = "in", dpi = 300)




##-------------------------------------------------------
## SEURAT ANALYSIS | UPDATE ANNOTATION
##-------------------------------------------------------
rm(list = ls())
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(ggrepel)
library(reticulate)
library(clustree)
library(scCustomize)

load("SLEEPSEQ_DATA_FILT_NORM_PCA_CLU_INTEGRATED.RData")
# sleepseq.integrated

DefaultAssay(object = sleepseq.integrated) <- "integrated"
Idents(sleepseq.integrated) <- "integrated_snn_res.1.2"

sleepseq.integrated@active.ident <- factor(x = sleepseq.integrated@active.ident, levels = seq(0, max(as.numeric(sleepseq.integrated$integrated_snn_res.1.2) - 1)))

sleepseq.meta <- as.data.frame(sleepseq.integrated@meta.data)
sleepseq.meta$CellBarcode <- row.names(sleepseq.meta)
sleepseq.meta$Cluster <- paste("Cluster", sprintf("%02d", as.numeric(sleepseq.integrated@active.ident) - 1), sep = "_")

annotation <- read.table("SEURAT_SLEEPSEQ_ANNOTATION_BY_CLUSTER.tsv", sep = "\t", header = TRUE)

sleepseq.meta.updated <- merge(sleepseq.meta, annotation, by = "Cluster")

annot1 <- sleepseq.meta.updated$Annotation1
names(annot1) <- sleepseq.meta.updated$CellBarcode

annot2 <- sleepseq.meta.updated$Annotation2
names(annot2) <- sleepseq.meta.updated$CellBarcode

sleepseq.integrated$Annotation1 <- annot1
sleepseq.integrated$Annotation2 <- annot2
sleepseq.integrated$Annotation3 <- paste(sleepseq.integrated$Annotation1, sleepseq.integrated$integrated_snn_res.1.2, sep = "_")
sleepseq.integrated$Annotation4 <- paste(sleepseq.integrated$Annotation2, sleepseq.integrated$integrated_snn_res.1.2, sep = "_")

Idents(sleepseq.integrated) <- "Annotation1"

p0 <- Cluster_Highlight_Plot(seurat_object = sleepseq.integrated, cluster_name = "UN", highlight_color = "grey50", background_color = "grey90", raster = TRUE)
ggsave(filename = "SEURAT_SLEEPSEQ_UMAP_CELLTYPE_0.pdf", plot = p0, width = 7, height = 6, units = "in", dpi = 300)

table(sleepseq.integrated$Annotation1)
# Astro  Endo    EX    IN Micro Oligo   OPC  Peri   PVM    UN  VLMC 
#  8788   721 42045 10062  4468  7973  2511   486   228  4481  2312

sleepseq.integrated <- subset(sleepseq.integrated, idents = "UN", invert = TRUE)

table(sleepseq.integrated$Annotation1)
# Astro  Endo    EX    IN Micro Oligo   OPC  Peri   PVM  VLMC 
#  8788   721 42045 10062  4468  7973  2511   486   228  2312


## UMAP PLOTS
p1 <- DimPlot_scCustom(seurat_object = sleepseq.integrated, pt.size = 0.01, group.by = "Annotation1", raster = TRUE, label = TRUE, label.size = 3, colors_use = DiscretePalette_scCustomize(num_colors = length(unique(sort(sleepseq.integrated$Annotation1))), palette = "varibow")) + NoLegend()
ggsave(filename = "SEURAT_SLEEPSEQ_UMAP_CELLTYPE_1.pdf", plot = p1, width = 6, height = 6, units = "in", dpi = 300)

p2 <- DimPlot_scCustom(seurat_object = sleepseq.integrated, pt.size = 0.01, group.by = "Annotation2", raster = TRUE, label = FALSE, colors_use = DiscretePalette_scCustomize(num_colors = length(unique(sort(sleepseq.integrated$Annotation2))), palette = "varibow")) #+ NoLegend()
ggsave(filename = "SEURAT_SLEEPSEQ_UMAP_CELLTYPE_2.pdf", plot = p2, width = 9, height = 6, units = "in", dpi = 300)

p2b <- DimPlot_scCustom(seurat_object = sleepseq.integrated, pt.size = 0.01, group.by = "Annotation2", raster = TRUE, label = TRUE, label.size = 2, colors_use = DiscretePalette_scCustomize(num_colors = length(unique(sort(sleepseq.integrated$Annotation2))), palette = "varibow")) + NoLegend()
ggsave(filename = "SEURAT_SLEEPSEQ_UMAP_CELLTYPE_2b.pdf", plot = p2b, width = 6, height = 6, units = "in", dpi = 300)

p3 <- DimPlot_scCustom(seurat_object = sleepseq.integrated, pt.size = 0.01, group.by = "Annotation4", raster = TRUE, label = TRUE, label.size = 2, colors_use = DiscretePalette_scCustomize(num_colors = length(unique(sort(sleepseq.integrated$Annotation4))), palette = "varibow")) + NoLegend()
ggsave(filename = "SEURAT_SLEEPSEQ_UMAP_CELLTYPE_3.pdf", plot = p3, width = 6, height = 6, units = "in", dpi = 300)



## VOILIN PLOTS
p6 <- Stacked_VlnPlot(seurat_object = sleepseq.integrated, group.by = "Annotation1", features = c("nCount_RNA", "nFeature_RNA"), x_lab_rotate = TRUE, ggplot_default_colors = TRUE)
ggsave(filename = "SEURAT_SLEEPSEQ_CELLTYPE_1_VIOLINPLOT_nUMI_nGENES.pdf", plot = p6, width = 8, height = 4, units = "in", dpi = 300)

p7 <- Stacked_VlnPlot(seurat_object = sleepseq.integrated, group.by = "Annotation2", features = c("nCount_RNA", "nFeature_RNA"), x_lab_rotate = TRUE, ggplot_default_colors = TRUE)
ggsave(filename = "SEURAT_SLEEPSEQ_CELLTYPE_2_VIOLINPLOT_nUMI_nGENES.pdf", plot = p7, width = 12, height = 4, units = "in", dpi = 300)

p8 <- Stacked_VlnPlot(seurat_object = sleepseq.integrated, group.by = "Annotation4", features = c("nCount_RNA", "nFeature_RNA"), x_lab_rotate = TRUE, ggplot_default_colors = TRUE)
ggsave(filename = "SEURAT_SLEEPSEQ_CELLTYPE_3_VIOLINPLOT_nUMI_nGENES.pdf", plot = p8, width = 16, height = 4, units = "in", dpi = 300)


## BARPLOT PLOTS
cellsPerCelltype <- as.data.frame.matrix(table(sleepseq.integrated$Annotation1, sleepseq.integrated@meta.data$Condition))
cellsPerCelltype$CellType <- row.names(cellsPerCelltype)
write.table(cellsPerCelltype, "SEURAT_SLEEPSEQ_CELLS_BY_CELLTYPE_1_BY_CONDITION.tsv", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

cellsPerCelltype2 <- melt(cellsPerCelltype)
colnames(cellsPerCelltype2) <- c("CELLTYPE", "CONDITION", "CELLS")
p4 <- ggplot(cellsPerCelltype2) +
        geom_bar(aes(x = CELLTYPE, y = CELLS, fill = CONDITION), stat = "identity", position = "fill") + 
        scale_y_continuous(labels = scales::percent_format()) + 
        theme_classic() + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
        NULL
ggsave(filename = "SEURAT_SLEEPSEQ_CELLS_BY_CELLTYPE_1_BY_CONDITION.pdf", plot = p4, width = 14, height = 4, units = "in", dpi = 300)


cellsPerCelltype <- as.data.frame.matrix(table(sleepseq.integrated$Annotation2, sleepseq.integrated@meta.data$Condition))
cellsPerCelltype$CellType <- row.names(cellsPerCelltype)
write.table(cellsPerCelltype, "SEURAT_SLEEPSEQ_CELLS_BY_CELLTYPE_2_BY_CONDITION.tsv", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

cellsPerCelltype2 <- melt(cellsPerCelltype)
colnames(cellsPerCelltype2) <- c("CELLTYPE", "CONDITION", "CELLS")
p4 <- ggplot(cellsPerCelltype2) +
        geom_bar(aes(x = CELLTYPE, y = CELLS, fill = CONDITION), stat = "identity", position = "fill") + 
        scale_y_continuous(labels = scales::percent_format()) + 
        theme_classic() + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
        NULL
ggsave(filename = "SEURAT_SLEEPSEQ_CELLS_BY_CELLTYPE_2_BY_CONDITION.pdf", plot = p4, width = 14, height = 4, units = "in", dpi = 300)


cellsPerCelltype <- as.data.frame.matrix(table(sleepseq.integrated$Annotation4, sleepseq.integrated@meta.data$Condition))
cellsPerCelltype$CellType <- row.names(cellsPerCelltype)
write.table(cellsPerCelltype, "SEURAT_SLEEPSEQ_CELLS_BY_CELLTYPE_4_BY_CONDITION.tsv", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

cellsPerCelltype2 <- melt(cellsPerCelltype)
colnames(cellsPerCelltype2) <- c("CELLTYPE", "CONDITION", "CELLS")
p4 <- ggplot(cellsPerCelltype2) +
        geom_bar(aes(x = CELLTYPE, y = CELLS, fill = CONDITION), stat = "identity", position = "fill") + 
        scale_y_continuous(labels = scales::percent_format()) + 
        theme_classic() + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
        NULL
ggsave(filename = "SEURAT_SLEEPSEQ_CELLS_BY_CELLTYPE_4_BY_CONDITION.pdf", plot = p4, width = 14, height = 4, units = "in", dpi = 300)


save(sleepseq.integrated, file = "SLEEPSEQ_DATA_FILT_NORM_PCA_CLU_INTEGRATED_ANNOTATED.RData")



##-------------------------------------------------------
## SEURAT ANALYSIS | UPDATE ANNOTATION
##-------------------------------------------------------
rm(list = ls())
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(ggrepel)
library(reticulate)
library(clustree)
library(scCustomize)

load("SLEEPSEQ_DATA_FILT_NORM_PCA_CLU_INTEGRATED_ANNOTATED.RData")
# sleepseq.integrated

sleepseq.integrated$Annotation5 <- gsub("EX_L2-3_IT", "EX_IT", sleepseq.integrated$Annotation2)
sleepseq.integrated$Annotation5 <- gsub("EX_L4-5_IT", "EX_IT", sleepseq.integrated$Annotation5)
sleepseq.integrated$Annotation5 <- gsub("EX_L5_IT", "EX_IT", sleepseq.integrated$Annotation5)
sleepseq.integrated$Annotation5 <- gsub("EX_L6_IT", "EX_IT", sleepseq.integrated$Annotation5)
sleepseq.integrated$Annotation5 <- gsub("EX_L5_ET", "EX_ET_NP_CT", sleepseq.integrated$Annotation5)
sleepseq.integrated$Annotation5 <- gsub("EX_L5-6_NP", "EX_ET_NP_CT", sleepseq.integrated$Annotation5)
sleepseq.integrated$Annotation5 <- gsub("EX_L6_CT", "EX_ET_NP_CT", sleepseq.integrated$Annotation5)
sleepseq.integrated$Annotation5 <- gsub("EX_L6b", "EX_ET_NP_CT", sleepseq.integrated$Annotation5)
sleepseq.integrated$Annotation5 <- gsub("IN_Lamp5", "IN", sleepseq.integrated$Annotation5)
sleepseq.integrated$Annotation5 <- gsub("IN_Pvalb", "IN", sleepseq.integrated$Annotation5)
sleepseq.integrated$Annotation5 <- gsub("IN_Sncg", "IN", sleepseq.integrated$Annotation5)
sleepseq.integrated$Annotation5 <- gsub("IN_Sst", "IN", sleepseq.integrated$Annotation5)
sleepseq.integrated$Annotation5 <- gsub("IN_Vip", "IN", sleepseq.integrated$Annotation5)


## UMAP PLOTS
p1 <- DimPlot_scCustom(seurat_object = sleepseq.integrated, pt.size = 0.01, group.by = "Annotation5", raster = TRUE, label = TRUE, label.size = 3, colors_use = DiscretePalette_scCustomize(num_colors = length(unique(sort(sleepseq.integrated$Annotation5))), palette = "varibow")) + NoLegend()
ggsave(filename = "SEURAT_SLEEPSEQ_UMAP_CELLTYPE_4.pdf", plot = p1, width = 6, height = 6, units = "in", dpi = 300)


```

### DEG Analysis
```
# module purge && module load shared slurm python/3.7.x-anaconda
# module load gcc/8.3.0
# module load hdf5_18/1.8.17
# module load R/4.0.2-gccmkl


rm(list = ls())
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(ggrepel)
library(reticulate)
library(plyr)
library(tidyverse)
library(tidyr)
library(patchwork)
library(scran)
library(scater)
library(SingleCellExperiment)
library(edgeR)

##-------------------------------------------------------
## DEG | PSEUDOBULK
##-------------------------------------------------------
load("/work/psychiatry/rgree2/SNRNASEQ/06_SEURAT_woAmbient_woDoublets/04_CS_SD_USE_THIS/SLEEPSEQ_DATA_FILT_NORM_PCA_CLU_INTEGRATED_ANNOTATED.RData")
# sleepseq.integrated

sleepseq.integrated$Batch_cDNA <- gsub("CS_A1|CS_A2|SD_A3|SD_A4", 1, sleepseq.integrated$ConditionSample)
sleepseq.integrated$Batch_cDNA <- gsub("CS_A5|CS_A6|SD_A7|SD_A8", 2, sleepseq.integrated$Batch_cDNA)

DefaultAssay(sleepseq.integrated) <- "RNA"
sleepseq.integrated <- NormalizeData(sleepseq.integrated)

table(sleepseq.integrated@meta.data$Condition)
#    CS    SD 
# 38048 41546

table(sleepseq.integrated@meta.data$ConditionSample)
# CS_A1 CS_A2 CS_A5 CS_A6 SD_A3 SD_A4 SD_A7 SD_A8 
#  6812  9496 11916  9824  9541 10447 12030  9528

table(sleepseq.integrated@meta.data$Annotation1)
# Astro  Endo    EX    IN Micro Oligo   OPC  Peri   PVM  VLMC 
#  8788   721 42045 10062  4468  7973  2511   486   228  2312

table(sleepseq.integrated@meta.data$Annotation2)
#       Astrocytes      Endothelial       EX_L2-3_IT       EX_L4-5_IT 
#             8788              721            12448             8660 
#         EX_L5_ET         EX_L5_IT       EX_L5-6_NP         EX_L6_CT 
#             2549             5488             1583             7527 
#         EX_L6_IT           EX_L6b         IN_Lamp5         IN_Pvalb 
#             3154              636              769             3071 
#          IN_Sncg           IN_Sst           IN_Vip        Microglia 
#             2771             2232             1219             4468 
# Oligodendrocytes              OPC        Pericytes              PVM 
#             7973             2511              486              228 
#             VLMC 
#             2312

table(sleepseq.integrated@meta.data$Annotation4)
#        Astrocytes_0       Astrocytes_25       Astrocytes_47      Endothelial_28 
#                7605                1078                 105                 591 
#      Endothelial_43       EX_L2-3_IT_11       EX_L2-3_IT_22        EX_L2-3_IT_4 
#                 130                2800                1386                4188 
#       EX_L2-3_IT_40       EX_L2-3_IT_44        EX_L2-3_IT_5       EX_L4-5_IT_16 
#                 159                 127                3788                2175 
#       EX_L4-5_IT_37       EX_L4-5_IT_59        EX_L4-5_IT_7        EX_L4-5_IT_9 
#                 189                  30                3204                3062 
#         EX_L5_ET_13         EX_L5_IT_15          EX_L5_IT_8       EX_L5-6_NP_19 
#                2549                2350                3138                1583 
#         EX_L6_CT_17          EX_L6_CT_2         EX_L6_CT_48         EX_L6_CT_51 
#                2128                5218                 103                  78 
#         EX_L6_IT_12         EX_L6_IT_29           EX_L6b_27         IN_Lamp5_26 
#                2610                 544                 636                 769 
#         IN_Pvalb_10         IN_Pvalb_36         IN_Pvalb_58          IN_Sncg_21 
#                2828                 206                  37                1464 
#          IN_Sncg_23           IN_Sst_18           IN_Sst_42           IN_Sst_50 
#                1307                1897                 150                  96 
#           IN_Sst_53           IN_Sst_60           IN_Vip_24         Microglia_3 
#                  61                  28                1219                4468 
#  Oligodendrocytes_1 Oligodendrocytes_32 Oligodendrocytes_33              OPC_14 
#                7152                 429                 392                2415 
#              OPC_49        Pericytes_30              PVM_35             VLMC_20 
#                  96                 486                 228                1540 
#             VLMC_31             VLMC_38             VLMC_45 
#                 475                 181                 116


Idents(sleepseq.integrated) <- "Annotation2"

table(sleepseq.integrated@active.ident)
#       Astrocytes      Endothelial       EX_L2-3_IT       EX_L4-5_IT 
#             8788              721            12448             8660 
#         EX_L5_ET         EX_L5_IT       EX_L5-6_NP         EX_L6_CT 
#             2549             5488             1583             7527 
#         EX_L6_IT           EX_L6b         IN_Lamp5         IN_Pvalb 
#             3154              636              769             3071 
#          IN_Sncg           IN_Sst           IN_Vip        Microglia 
#             2771             2232             1219             4468 
# Oligodendrocytes              OPC        Pericytes              PVM 
#             7973             2511              486              228 
#             VLMC 
#             2312

for (myclu in unique(sort(sleepseq.integrated@active.ident)))
	{
    ## select a cell-type    
	cluSelected <- myclu
	print(cluSelected)

	sleepseq.selected <- subset(x = sleepseq.integrated, idents = c(cluSelected))

    ## create SCE object for all genes
    sleepseqSCE <- as.SingleCellExperiment(sleepseq.selected)

    sleepseqGroups <- colData(sleepseqSCE)[, c("ConditionSample", "Condition")]
    sleepseqPseudoSCE <- sumCountsAcrossCells(sleepseqSCE, sleepseqGroups)
    sleepseqAggMat <- sleepseqPseudoSCE@assays@data$sum
    colnames(sleepseqAggMat) <- colData(sleepseqPseudoSCE)[['ConditionSample']]

    ## filter genes
    csExp = which(rowSums(sleepseqAggMat[, grepl('CS', colnames(sleepseqAggMat))] <= 10) == 0) %>% names
    sdExp = which(rowSums(sleepseqAggMat[, grepl('SD', colnames(sleepseqAggMat))] <= 10) == 0) %>% names
    expgns = Reduce(union, list(csExp, sdExp))

    ## create SCE object for filtered genes
    DefaultAssay(sleepseq.selected) <- "RNA"
    selSCE <- as.SingleCellExperiment(sleepseq.selected)

    selGroups <- colData(selSCE)[, c("ConditionSample", "Condition", "Sample", "Batch_cDNA")]
    selPseudoSCE <- sumCountsAcrossCells(selSCE, selGroups)
    selGroups <- colData(selPseudoSCE)[, c("ConditionSample", "Condition", "Sample", "Batch_cDNA")]
    selGroups$ConditionSample <- factor(selGroups$ConditionSample)
    selGroups$Condition <- factor(selGroups$Condition)
    selGroups$Sample <- factor(selGroups$Sample)
    selGroups$Batch_cDNA <- factor(selGroups$Batch_cDNA)

    selAggMat <- selPseudoSCE@assays@data$sum
    colnames(selAggMat) = colData(selPseudoSCE)[['ConditionSample']]
    selAggMat <- selAggMat[expgns,]

    ## perform differential gene expression
    selDGEL <- DGEList(counts = selAggMat)
    selDGEL <- calcNormFactors(selDGEL)
    selDesign <- model.matrix(~Batch_cDNA + Condition, data = selGroups)
    selDGEL <- estimateDisp(selDGEL, selDesign)

    selFit <- glmFit(selDGEL,selDesign)
    selLrt <- glmLRT(selFit)
    selRes <- topTags(selLrt, n = nrow(selAggMat), sort.by = 'PValue') %>% as.data.frame
    degTableSel <- subset(selRes, FDR <= 0.05)

    ## write dge table
    write.table(selRes, paste("SEURAT_SLEEPSEQ_DATA_", cluSelected, "_SD_x_CS_DEG_TABLE_PSEUDOBULK_FULL.txt", sep = ""), row.names = T, col.names = T, quote = F, sep = "\t")

	## average gene expression
	sleepseq.selected$CellType <- paste(Idents(sleepseq.selected), sleepseq.selected$Condition, sep = "_")
	Idents(sleepseq.selected) <- "CellType"

	avg.exp <- as.data.frame(log1p(AverageExpression(sleepseq.selected, verbose = FALSE)$RNA))
	avg.exp$Gene <- row.names(avg.exp)

    # avgExp <- avg.exp[row.names(avg.exp) %in% row.names(selRes),]
    degavg <- merge(avg.exp, selRes, by = "row.names", all.y = TRUE)
    row.names(degavg) <- degavg$Row.names
    degavg$Row.names <- NULL

    ## write average gene expression table
	write.table(avg.exp, paste("SEURAT_SLEEPSEQ_DATA_", cluSelected, "_AVG_EXP.txt", sep = ""), row.names = T, col.names = T, quote = F, sep = "\t")

    ## save dge table and avg gene expression
    save(selRes, avg.exp, degavg, file = paste("SEURAT_SLEEPSEQ_DATA_", cluSelected, "_PSEUDOBULK.RData", sep = ""))
	}


```


### DEG Analysis Plots
```
# module purge && module load shared slurm python/3.7.x-anaconda
# module load gcc/8.3.0
# module load hdf5_18/1.8.17
# module load R/4.0.2-gccmkl

rm(list = ls())
# .libPaths("/home2/akulk1/R/x86_64-redhat-linux-gnu-library/3.5")
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(ggrepel)
library(reticulate)
library(RColorBrewer)
library(UpSetR)

##-------------------------------------------------------
## DEGs
##-------------------------------------------------------

filterDEGs <- function(degtable, pvalcutoff, l2fccutoff) 
	{
	deg.cs.sd.temp <- read.table(degtable, sep = "\t", header = TRUE)
	deg.cs.sd <- deg.cs.sd.temp
	deg.cs.sd.sig.up <- deg.cs.sd[deg.cs.sd$FDR <= pvalcutoff & deg.cs.sd$logFC >= l2fccutoff,] 
	deg.cs.sd.sig.dn <- deg.cs.sd[deg.cs.sd$FDR <= pvalcutoff & deg.cs.sd$logFC <= -l2fccutoff,]
	return(list("UP" = row.names(deg.cs.sd.sig.up), "DOWN" = row.names(deg.cs.sd.sig.dn)))
	}

# degtable <- "SEURAT_SLEEPSEQ_DATA_Astrocytes_SD_x_CS_DEG_TABLE_WILCOX_FULL.txt"
# degtable <- "SEURAT_SLEEPSEQ_DATA_Astrocytes_SD_x_CS_DEG_TABLE_PSEUDOBULK_FULL.txt"
pvalcutoff <- 0.05
l2fccutoff <- 0.1375 # 10% difference

astro_list <- filterDEGs("SEURAT_SLEEPSEQ_DATA_Astrocytes_SD_x_CS_DEG_TABLE_PSEUDOBULK_FULL.txt", pvalcutoff, l2fccutoff)
endo_list <- filterDEGs("SEURAT_SLEEPSEQ_DATA_Endothelial_SD_x_CS_DEG_TABLE_PSEUDOBULK_FULL.txt", pvalcutoff, l2fccutoff)
exl23it_list <- filterDEGs("SEURAT_SLEEPSEQ_DATA_EX_L2-3_IT_SD_x_CS_DEG_TABLE_PSEUDOBULK_FULL.txt", pvalcutoff, l2fccutoff)
exl45it_list <- filterDEGs("SEURAT_SLEEPSEQ_DATA_EX_L4-5_IT_SD_x_CS_DEG_TABLE_PSEUDOBULK_FULL.txt", pvalcutoff, l2fccutoff)
exl5et_list <- filterDEGs("SEURAT_SLEEPSEQ_DATA_EX_L5_ET_SD_x_CS_DEG_TABLE_PSEUDOBULK_FULL.txt", pvalcutoff, l2fccutoff)
exl5it_list <- filterDEGs("SEURAT_SLEEPSEQ_DATA_EX_L5_IT_SD_x_CS_DEG_TABLE_PSEUDOBULK_FULL.txt", pvalcutoff, l2fccutoff)
exl56np_list <- filterDEGs("SEURAT_SLEEPSEQ_DATA_EX_L5-6_NP_SD_x_CS_DEG_TABLE_PSEUDOBULK_FULL.txt", pvalcutoff, l2fccutoff)
exl6b_list <- filterDEGs("SEURAT_SLEEPSEQ_DATA_EX_L6b_SD_x_CS_DEG_TABLE_PSEUDOBULK_FULL.txt", pvalcutoff, l2fccutoff)
exl6ct_list <- filterDEGs("SEURAT_SLEEPSEQ_DATA_EX_L6_CT_SD_x_CS_DEG_TABLE_PSEUDOBULK_FULL.txt", pvalcutoff, l2fccutoff)
exl6it_list <- filterDEGs("SEURAT_SLEEPSEQ_DATA_EX_L6_IT_SD_x_CS_DEG_TABLE_PSEUDOBULK_FULL.txt", pvalcutoff, l2fccutoff)
inlamp5_list <- filterDEGs("SEURAT_SLEEPSEQ_DATA_IN_Lamp5_SD_x_CS_DEG_TABLE_PSEUDOBULK_FULL.txt", pvalcutoff, l2fccutoff)
inpvalb_list <- filterDEGs("SEURAT_SLEEPSEQ_DATA_IN_Pvalb_SD_x_CS_DEG_TABLE_PSEUDOBULK_FULL.txt", pvalcutoff, l2fccutoff)
insncg_list <- filterDEGs("SEURAT_SLEEPSEQ_DATA_IN_Sncg_SD_x_CS_DEG_TABLE_PSEUDOBULK_FULL.txt", pvalcutoff, l2fccutoff)
insst_list <- filterDEGs("SEURAT_SLEEPSEQ_DATA_IN_Sst_SD_x_CS_DEG_TABLE_PSEUDOBULK_FULL.txt", pvalcutoff, l2fccutoff)
invip_list <- filterDEGs("SEURAT_SLEEPSEQ_DATA_IN_Vip_SD_x_CS_DEG_TABLE_PSEUDOBULK_FULL.txt", pvalcutoff, l2fccutoff)
micro_list <- filterDEGs("SEURAT_SLEEPSEQ_DATA_Microglia_SD_x_CS_DEG_TABLE_PSEUDOBULK_FULL.txt", pvalcutoff, l2fccutoff)
oligo_list <- filterDEGs("SEURAT_SLEEPSEQ_DATA_Oligodendrocytes_SD_x_CS_DEG_TABLE_PSEUDOBULK_FULL.txt", pvalcutoff, l2fccutoff)
opc_list <- filterDEGs("SEURAT_SLEEPSEQ_DATA_OPC_SD_x_CS_DEG_TABLE_PSEUDOBULK_FULL.txt", pvalcutoff, l2fccutoff)
peri_list <- filterDEGs("SEURAT_SLEEPSEQ_DATA_Pericytes_SD_x_CS_DEG_TABLE_PSEUDOBULK_FULL.txt", pvalcutoff, l2fccutoff)
pvm_list <- filterDEGs("SEURAT_SLEEPSEQ_DATA_PVM_SD_x_CS_DEG_TABLE_PSEUDOBULK_FULL.txt", pvalcutoff, l2fccutoff)
vlmc_list <- filterDEGs("SEURAT_SLEEPSEQ_DATA_VLMC_SD_x_CS_DEG_TABLE_PSEUDOBULK_FULL.txt", pvalcutoff, l2fccutoff)


degcountlist <- list("astro_up" = length(astro_list$UP), "astro_dn" = length(astro_list$DOWN),
					 "endo_up" = length(endo_list$UP), "endo_dn" = length(endo_list$DOWN),
					 "exl23it_up" = length(exl23it_list$UP), "exl23it_dn" = length(exl23it_list$DOWN),
					 "exl45it_up" = length(exl45it_list$UP), "exl45it_dn" = length(exl45it_list$DOWN),
					 "exl5et_up" = length(exl5et_list$UP), "exl5et_dn" = length(exl5et_list$DOWN),
					 "exl5it_up" = length(exl5it_list$UP), "exl5it_dn" = length(exl5it_list$DOWN),
					 "exl56np_up" = length(exl56np_list$UP), "exl56np_dn" = length(exl56np_list$DOWN),
					 "exl6b_up" = length(exl6b_list$UP), "exl6b_dn" = length(exl6b_list$DOWN),
					 "exl6ct_up" = length(exl6ct_list$UP), "exl6ct_dn" = length(exl6ct_list$DOWN),
					 "exl6it_up" = length(exl6it_list$UP), "exl6it_dn" = length(exl6it_list$DOWN),
					 "inlamp5_up" = length(inlamp5_list$UP), "inlamp5_dn" = length(inlamp5_list$DOWN),
					 "inpvalb_up" = length(inpvalb_list$UP), "inpvalb_dn" = length(inpvalb_list$DOWN),
					 "insncg_up" = length(insncg_list$UP), "insncg_dn" = length(insncg_list$DOWN),
					 "insst_up" = length(insst_list$UP), "insst_dn" = length(insst_list$DOWN),
					 "invip_up" = length(invip_list$UP), "invip_dn" = length(invip_list$DOWN),
					 "micro_up" = length(micro_list$UP), "micro_dn" = length(micro_list$DOWN),
					 "oligo_up" = length(oligo_list$UP), "oligo_dn" = length(oligo_list$DOWN),
					 "opc_up" = length(opc_list$UP), "opc_dn" = length(opc_list$DOWN),
					 "peri_up" = length(peri_list$UP), "peri_dn" = length(peri_list$DOWN),
					 "pvm_up" = length(pvm_list$UP), "pvm_dn" = length(pvm_list$DOWN),
					 "vlmc_up" = length(vlmc_list$UP), "vlmc_dn" = length(vlmc_list$DOWN))

degup <- as.data.frame(t(as.data.frame(degcountlist[grepl("_up", names(degcountlist))])))
degup$celltype <- gsub("_up", "", row.names(degup))
degdn <- as.data.frame(t(as.data.frame(degcountlist[grepl("_dn", names(degcountlist))])))
degdn$celltype <- gsub("_dn", "", row.names(degdn))

degcounts <- merge(degup, degdn, by = "celltype", all.x = TRUE, incomparables = 0)
colnames(degcounts) <- c("CELLTYPE", "UP", "DOWN")
degcounts[is.na(degcounts)] <- 0
write.table(degcounts, "SLEEPSEQ_DEG_COMPARISON_L2FC_0.1375.tsv", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

p1 <- ggplot(degcounts) + 
	  geom_hline(yintercept = 0) + 
	  geom_segment(aes(x = CELLTYPE, xend = CELLTYPE, y = 0, yend = UP)) +
	  geom_point(data = degcounts, aes(x = CELLTYPE, y = UP, color = CELLTYPE)) + 
	  geom_segment(aes(x = CELLTYPE, xend = CELLTYPE, y = 0, yend = -DOWN)) +
	  geom_point(data = degcounts, aes(x = CELLTYPE, y = -DOWN, color = CELLTYPE)) + 	  
	  theme_bw() + 
	  ylim(-500,500) +  #ylim(-150,150) + 
	  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
	  theme(legend.position = "none") + 
	  xlab("CellTypes") + 
	  ylab("# DEGs") + 
	  NULL
ggsave(filename = "SLEEPSEQ_DEG_COMPARISON_L2FC_0.1375.pdf", plot = p1, width = 8, height = 4, units = "in", dpi = 300, useDingbats = FALSE)


```


### DEG Enrichment
```
load("/work/psychiatry/rgree2/SNRNASEQ/06_SEURAT_woAmbient_woDoublets/04_CS_SD_USE_THIS/SLEEPSEQ_DATA_FILT_NORM_PCA_CLU_INTEGRATED_ANNOTATED.RData")
# sleepseq.integrated


##-------------------------------------------------------
## PREPARE REFERENCE GENE LISTS
## ASD genes
load("ASD_SFARI.RData")
# GeneSets
asd_human <- as.vector(GeneSets$ASD$Gene)
asd_mouse <- na.omit(convert_human_to_mouse_symbols(symbols = asd_human))  # 1045 genes

asdscored_human <- as.vector(GeneSets$ASD_Scored$Gene)
asdscored_mouse <- na.omit(convert_human_to_mouse_symbols(symbols = asdscored_human))

rm(GeneSets)


## Chris Cowan Mef2c
load("GeneSets_Mef2c_cKO.RData")
# GeneSets
mef2c_cko_up <- as.vector(GeneSets$Mef2c_cKO_UpReg$Gene)  # 478 genes
mef2c_cko_dn <- as.vector(GeneSets$Mef2c_cKO_DownReg$Gene)  # 598 genes

rm(GeneSets)


## DEGs comparing WT to HDAC4cn
hd4cn_deg <- read.table("DEG_HD4CN_X_MCH_CTX.txt", header = TRUE, sep = "\t")
hd4cn_up <- hd4cn_deg[hd4cn_deg$FDR <= 0.05 & hd4cn_deg$log2FC >= 0.3, "Gene"]  # 323 genes
hd4cn_dn <- hd4cn_deg[hd4cn_deg$FDR <= 0.05 & hd4cn_deg$log2FC <= -0.3, "Gene"]  # 313 genes


## Synaptic plasticity genes
sscgenesTemp <- as.data.frame(scan("Synaptic_Shaping_Components.txt", what = "", sep = "\n")) #279 genes
sscgenes <- as.vector(unique(sscgenesTemp[,1])) ## 865 genes

# ## create a list for all reference gene lists
# myref <- list("01_MEF2C_CKO_DEG_UP" = mef2c_cko_up, 
#               "02_HD4CN_DEG_UP" = hd4cn_up, 
#               "03_MEF2C_CKO_DEG_DOWN" = mef2c_cko_dn, 
#               "04_HD4CN_DEG_DOWN" = hd4cn_dn, 
#               "05_SSC_GENES" = sscgenes,
#               "06_ASD_ALL_GENES" = asd_mouse)

# save(myref, file = "SLEEPSEQ_REFERENCE_GENE_LISTS.RData")

myref <- list("01_SSC_GENES" = sscgenes,
              "02_ASD_GENES" = asd_mouse,
              "03_MEF2C_CKO_DEG" = c(mef2c_cko_up, mef2c_cko_dn),
              "04_HD4CN_DEG" = c(hd4cn_up, hd4cn_dn))


# myref <- list("01_MEF2C_CKO_DEG" = c(mef2c_cko_up, mef2c_cko_dn),
#               "02_HD4CN_DEG" = c(hd4cn_up, hd4cn_dn),
#               "03_SSC_GENES" = sscgenes,
#               "04_ASD_ALL_GENES" = asd_mouse)



##-------------------------------------------------------
## READ & FILTER DEG TABLES
##-------------------------------------------------------
## READ & FILTER DEG TABLES
filterDEGs <- function(degtable, pvalcutoff, l2fccutoff) 
	{
	deg.cs.sd.temp <- read.table(degtable, sep = "\t", header = TRUE)
	deg.cs.sd <- deg.cs.sd.temp
	deg.cs.sd.sig.up <- deg.cs.sd[deg.cs.sd$FDR <= pvalcutoff & deg.cs.sd$logFC >= l2fccutoff,] 
	deg.cs.sd.sig.dn <- deg.cs.sd[deg.cs.sd$FDR <= pvalcutoff & deg.cs.sd$logFC <= -l2fccutoff,]
	return(list("UP" = row.names(deg.cs.sd.sig.up), "DOWN" = row.names(deg.cs.sd.sig.dn)))
	}

pvalcutoff <- 0.05
l2fccutoff <- 0.1375 # 10% difference

astro_list <- filterDEGs("./../SEURAT_SLEEPSEQ_DATA_Astrocytes_SD_x_CS_DEG_TABLE_PSEUDOBULK_FULL.txt", pvalcutoff, l2fccutoff)
endo_list <- filterDEGs("./../SEURAT_SLEEPSEQ_DATA_Endothelial_SD_x_CS_DEG_TABLE_PSEUDOBULK_FULL.txt", pvalcutoff, l2fccutoff)
exl23it_list <- filterDEGs("./../SEURAT_SLEEPSEQ_DATA_EX_L2-3_IT_SD_x_CS_DEG_TABLE_PSEUDOBULK_FULL.txt", pvalcutoff, l2fccutoff)
exl45it_list <- filterDEGs("./../SEURAT_SLEEPSEQ_DATA_EX_L4-5_IT_SD_x_CS_DEG_TABLE_PSEUDOBULK_FULL.txt", pvalcutoff, l2fccutoff)
exl5et_list <- filterDEGs("./../SEURAT_SLEEPSEQ_DATA_EX_L5_ET_SD_x_CS_DEG_TABLE_PSEUDOBULK_FULL.txt", pvalcutoff, l2fccutoff)
exl5it_list <- filterDEGs("./../SEURAT_SLEEPSEQ_DATA_EX_L5_IT_SD_x_CS_DEG_TABLE_PSEUDOBULK_FULL.txt", pvalcutoff, l2fccutoff)
exl56np_list <- filterDEGs("./../SEURAT_SLEEPSEQ_DATA_EX_L5-6_NP_SD_x_CS_DEG_TABLE_PSEUDOBULK_FULL.txt", pvalcutoff, l2fccutoff)
exl6b_list <- filterDEGs("./../SEURAT_SLEEPSEQ_DATA_EX_L6b_SD_x_CS_DEG_TABLE_PSEUDOBULK_FULL.txt", pvalcutoff, l2fccutoff)
exl6ct_list <- filterDEGs("./../SEURAT_SLEEPSEQ_DATA_EX_L6_CT_SD_x_CS_DEG_TABLE_PSEUDOBULK_FULL.txt", pvalcutoff, l2fccutoff)
exl6it_list <- filterDEGs("./../SEURAT_SLEEPSEQ_DATA_EX_L6_IT_SD_x_CS_DEG_TABLE_PSEUDOBULK_FULL.txt", pvalcutoff, l2fccutoff)
inlamp5_list <- filterDEGs("./../SEURAT_SLEEPSEQ_DATA_IN_Lamp5_SD_x_CS_DEG_TABLE_PSEUDOBULK_FULL.txt", pvalcutoff, l2fccutoff)
inpvalb_list <- filterDEGs("./../SEURAT_SLEEPSEQ_DATA_IN_Pvalb_SD_x_CS_DEG_TABLE_PSEUDOBULK_FULL.txt", pvalcutoff, l2fccutoff)
insncg_list <- filterDEGs("./../SEURAT_SLEEPSEQ_DATA_IN_Sncg_SD_x_CS_DEG_TABLE_PSEUDOBULK_FULL.txt", pvalcutoff, l2fccutoff)
insst_list <- filterDEGs("./../SEURAT_SLEEPSEQ_DATA_IN_Sst_SD_x_CS_DEG_TABLE_PSEUDOBULK_FULL.txt", pvalcutoff, l2fccutoff)
invip_list <- filterDEGs("./../SEURAT_SLEEPSEQ_DATA_IN_Vip_SD_x_CS_DEG_TABLE_PSEUDOBULK_FULL.txt", pvalcutoff, l2fccutoff)
micro_list <- filterDEGs("./../SEURAT_SLEEPSEQ_DATA_Microglia_SD_x_CS_DEG_TABLE_PSEUDOBULK_FULL.txt", pvalcutoff, l2fccutoff)
oligo_list <- filterDEGs("./../SEURAT_SLEEPSEQ_DATA_Oligodendrocytes_SD_x_CS_DEG_TABLE_PSEUDOBULK_FULL.txt", pvalcutoff, l2fccutoff)
opc_list <- filterDEGs("./../SEURAT_SLEEPSEQ_DATA_OPC_SD_x_CS_DEG_TABLE_PSEUDOBULK_FULL.txt", pvalcutoff, l2fccutoff)
peri_list <- filterDEGs("./../SEURAT_SLEEPSEQ_DATA_Pericytes_SD_x_CS_DEG_TABLE_PSEUDOBULK_FULL.txt", pvalcutoff, l2fccutoff)
pvm_list <- filterDEGs("./../SEURAT_SLEEPSEQ_DATA_PVM_SD_x_CS_DEG_TABLE_PSEUDOBULK_FULL.txt", pvalcutoff, l2fccutoff)
vlmc_list <- filterDEGs("./../SEURAT_SLEEPSEQ_DATA_VLMC_SD_x_CS_DEG_TABLE_PSEUDOBULK_FULL.txt", pvalcutoff, l2fccutoff)


sleepdeg <-list("01_EX_L23IT_DEG_UP" = exl23it_list$UP, "02_EX_L23IT_DEG_DOWN" = exl23it_list$DOWN,
                "03_EX_L45IT_DEG_UP" = exl45it_list$UP, "04_EX_L45IT_DEG_DOWN" = exl45it_list$DOWN,
                "05_EX_L5IT_DEG_UP" = exl5it_list$UP, "06_EX_L5IT_DEG_DOWN" = exl5it_list$DOWN,
                "07_EX_L6IT_DEG_UP" = exl6it_list$UP, "08_EX_L6IT_DEG_DOWN" = exl6it_list$DOWN,
                "09_EX_L5ET_DEG_UP" = exl5et_list$UP, "10_EX_L5ET_DEG_DOWN" = exl5et_list$DOWN,
                "11_EX_L56NP_DEG_UP" = exl56np_list$UP, "12_EX_L56NP_DEG_DOWN" = exl56np_list$DOWN,
                "13_EX_L6B_DEG_UP" = exl6b_list$UP, "14_EX_L6B_DEG_DOWN" = exl6b_list$DOWN,
                "15_EX_L6CT_DEG_UP" = exl6ct_list$UP, "16_EX_L6CT_DEG_DOWN" = exl6ct_list$DOWN,
                "17_IN_Lamp5_DEG_UP" = inlamp5_list$UP, "18_IN_Lamp5_DEG_DOWN" = inlamp5_list$DOWN,
                "19_IN_Pvalb_DEG_UP" = inpvalb_list$UP, "20_IN_Pvalb_DEG_DOWN" = inpvalb_list$DOWN,
                "21_IN_Sncg_DEG_UP" = insncg_list$UP, "22_IN_Sncg_DEG_DOWN" = insncg_list$DOWN,
                "23_IN_Sst_DEG_UP" = insst_list$UP, "24_IN_Sst_DEG_DOWN" = insst_list$DOWN,
                "25_IN_Vip_DEG_UP" = invip_list$UP, "26_IN_Vip_DEG_DOWN" = invip_list$DOWN)


##-------------------------------------------------------
## Perform Super Exact Test
library(SuperExactTest)
library(ggplot2)
library(ggnewscale)
library(rlist)
library(wesanderson)

setout <- list()
setoutAll <- list()
i <- 0

for(mod in names(sleepdeg))
 {
 i <- i + 1
 newGenes <- list.append(myref, mod = sleepdeg[[mod]])
 names(newGenes)[[length(newGenes)]] <- mod
 setres <- supertest(newGenes, n = 13578, degree = 2)  # 13578 is background number of genes
 setresData <- as.data.frame(summary(setres)$Table)
 setoutAll[[i]] <- setresData
 setresDataSel <- setresData[grep(mod, setresData$Intersections),c(1, 3, 5, 6)]
 setout[[i]] <- setresDataSel
 print(paste(mod, i, sep = " "))
 }

names(setoutAll) <- names(sleepdeg)
names(setout) <- names(sleepdeg)

setoutData <- Reduce("rbind", setout)

save(setout, setoutData, setoutAll, file = "SLEEPSEQ_PSEUDOBULK_DGE_CLASS_SET.RData")


## REORGANIZE DATA FOR PLOTS
setoutDataTemp <- as.data.frame(matrix(unlist(strsplit(setoutData$Intersections, " & ")), byrow = T, ncol = 2))
colnames(setoutDataTemp) <- c("REF_GENES", "CELLTYPE_DEG")
setoutDataTemp$Intersections <- setoutData$Intersections

setoutData2 <- merge(setoutData, setoutDataTemp, by = "Intersections")
setoutData2$NegLog10 <- -log10(setoutData2$P.value + 10^-10)

write.table(setoutData2, "SLEEPSEQ_PSEUDOBULK_DGE_CLASS_SET_FOR_PLOT.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

setoutData3 <- setoutData2
setoutData3$NewPval <- setoutData3$P.value
setoutData3$NewPval[setoutData3$NewPval > 0.05] <- 1
setoutData3$NegLog10 <- -log10(setoutData3$NewPval + 10^-10)

##-------------------------------------------------------
## HEATMAP
pAll <- ggplot(data = setoutData3) + 
        geom_tile(data = setoutData3[setoutData3$REF_GENES == "03_MEF2C_CKO_DEG",], aes(x = CELLTYPE_DEG, y = REF_GENES, fill = NegLog10), color = "white", linewidth = 0.5) + 
        # scale_fill_viridis_c(option = "turbo") + 
        # scale_fill_gradientn(colours = c("red3", "gold", "green4")) +
        scale_fill_gradientn(colours = c("grey95", "green4", "gold", "red3")) +
        theme(panel.background = element_blank(), # bg of the panel
              panel.grid.major = element_blank(), # get rid of major grid
              panel.grid.minor = element_blank(), # get rid of minor grid
              axis.title=element_blank(),
              panel.grid=element_blank(),
              axis.text.x=element_text(size=10, angle=90, vjust=0.5, hjust=1),
              axis.ticks=element_blank(),
              axis.text.y=element_text(size=10)) +
        NULL
pAll <- pAll + 
        geom_tile(data = setoutData3[setoutData3$REF_GENES == "04_HD4CN_DEG",], aes(x = CELLTYPE_DEG, y = REF_GENES, fill = NegLog10), color = "white", linewidth = 0.5) + 
        geom_tile(data = setoutData3[setoutData3$REF_GENES == "02_ASD_GENES",], aes(x = CELLTYPE_DEG, y = REF_GENES, fill = NegLog10), color = "white", linewidth = 0.5) + 
        geom_tile(data = setoutData3[setoutData3$REF_GENES == "01_SSC_GENES",], aes(x = CELLTYPE_DEG, y = REF_GENES, fill = NegLog10), color = "white", linewidth = 0.5) + 
        NULL
ggsave(filename = "SLEEPSEQ_PSEUDOBULK_DGE_CLASS_SET_FOR_PLOT.pdf", plot = pAll, width = 22, height = 5, units = "in", dpi = 300)

```
