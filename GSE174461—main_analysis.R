setwd("~/desktop/Kharazmi/")

GSM5311807 <-  read.delim("~/desktop/Kharazmi/GSE174461_Proccesed/single_cell_count_mateix_GSM5311807.csv")
GSM5311808 <-  read.delim("~/desktop/Kharazmi/GSE174461_Proccesed/single_cell_count_mateix_GSM5311808.csv")
GSM5311809 <-  read.delim("~/desktop/Kharazmi/GSE174461_Proccesed/single_cell_count_mateix_GSM5311809.csv")
GSM5311810 <-  read.delim("~/desktop/Kharazmi/GSE174461_Proccesed/single_cell_count_mateix_GSM5311810.csv")
GSM5311811 <-  read.delim("~/desktop/Kharazmi/GSE174461_Proccesed/single_cell_count_mateix_GSM5311811.csv")
GSM5311812 <-  read.delim("~/desktop/Kharazmi/GSE174461_Proccesed/single_cell_count_mateix_GSM5311812.csv")

library(GEOquery)
geo_id <- "GSE174461"
geo_data <- getGEO(GEO = geo_id, GSEMatrix = T)

# 1. Constructing object
library(dplyr)
library(Seurat)
library(patchwork)
pbmc.data <- cbind(GSM5311807,GSM5311808,GSM5311809,GSM5311810,GSM5311811,GSM5311812)
dim(pbmc.data)
pbmc <-  CreateSeuratObject(counts =  pbmc.data, project = "Breast Cancer", min.cells = 3, min.features = 200)
pbmc

# 2. QC
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, add.noise = T)

# Now we can trim out data based on violin plot
pbmc <-  subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 3900 & nCount_RNA < 15000 & percent.mt < 5)

library(ggplot2)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")   + geom_smooth(method = "lm")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = "lm")
plot1+plot2


# Normalizing Data
pbmc <- NormalizeData(pbmc)

# 4. Identify highly variable features
pbmc  <-  FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 2000)

top10 <-  head(VariableFeatures(pbmc), 10)
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)
plot1 + plot2

pdf("VariablePlot — GSE174461.pdf", width = 15, height = 15)
plot2
dev.off()

# Scaling the data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc[["RNA"]]@scale.data

# Perform linear dimensional reduction
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")


# Determine the ‘dimensionality’ of the dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)

JackStrawPlot(pbmc, dims = 1:15)

ElbowPlot(pbmc)

pbmc <- FindNeighbors(pbmc, dims = 1:15)
pbmc <- FindClusters(pbmc, resolution = 0.2)
# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)

# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
pbmc <- RunUMAP(pbmc, dims = 1:15)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(pbmc, reduction = "umap")





