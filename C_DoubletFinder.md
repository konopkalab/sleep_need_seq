## Remove doublets using DoubletFinder

### CS
```
## load modules before starting R session
# module purge && module load shared slurm python/3.7.x-anaconda
# module load gcc/8.3.0
# module load hdf5_18/1.8.17
# module load R/4.0.2-gccmkl

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






```

### SD
```

```
