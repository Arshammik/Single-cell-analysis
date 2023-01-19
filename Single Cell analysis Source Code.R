# load Library
library(Seurat)
library(tidyverse)


# load the data
nsclc.sparse.m <- Read10X_h5(filename = "~/desktop/Kharazmi/Single Cell/Data/40k_NSCLC_DTC_3p_HT_nextgem_Multiplex_count_raw_feature_bc_matrix.h5")
str(nsclc.sparse.m)
cts <- nsclc.sparse.m$`Gene Expression`

# creating a Seurat object
nsclc.seurat.obj <- CreateSeuratObject(counts = cts, project = "NSCLC", min.cells = 3, min.features = 200) 
# at least we saw a gene a 3 cells and at least 200 gene expressed in each cell 
str(nsclc.seurat.obj)
# Features = Genes

# QC
# 1. Number of Mitocondrial reads
nsclc.seurat.obj[["percent.mt"]] <-  PercentageFeatureSet(nsclc.seurat.obj,pattern = "^MT-") #load it as a `percent` in our meta data
view(nsclc.seurat.obj@meta.data)

VlnPlot(nsclc.seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3 )
FeatureScatter(nsclc.seurat.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = "lm")

# 2. Filtering (Modified)
nsclc.seurat.obj <- subset(nsclc.seurat.obj, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & nCount_RNA < 2500 & percent.mt < 5)
# These filtration usually set by the violin plot results

# 3. Normalizing Data
# nsclc.seurat.obj <-  NormalizeData(nsclc.seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000) #defult value
nsclc.seurat.obj <-  NormalizeData(nsclc.seurat.obj)
view(nsclc.seurat.obj@data)

# 4. Identify highly variable features
nsclc.seurat.obj <- FindVariableFeatures(nsclc.seurat.obj, selection.method = "vst", nfeatures = 2000)
# Kind of a DEGs because it finds the most variable features(genes) in the data set and higher the variation higher the dissimilarity between expression among samples

top20  <- head(VariableFeatures(nsclc.seurat.obj), 20)
Plot20 <- VariableFeaturePlot(nsclc.seurat.obj)
LabelPoints(plot = Plot20, points = top20, repel = T)

# 5. Scaling 
all.genes <-  rownames(nsclc.seurat.obj)
nsclc.seurat.obj <- ScaleData(nsclc.seurat.obj, features = all.genes)
# mean expression = 0 (centering) and standard deviation = 1 (scaling)

# 6. Linear dimensional reduction
nsclc.seurat.obj <- RunPCA(nsclc.seurat.obj, features = VariableFeatures(object = nsclc.seurat.obj))

print(nsclc.seurat.obj[["pca"]], dims = 1:5, nfeatures = 5)
DimHeatmap(nsclc.seurat.obj, dims = 1, cells = 500, balanced = T)

# 7. determine dimensional of the data
ElbowPlot(nsclc.seurat.obj)

# 8. Clustring
nsclc.seurat.obj <- FindNeighbors(nsclc.seurat.obj, dims = 1:15)

# undrestanding reseloutions
nsclc.seurat.obj <-  FindClusters(nsclc.seurat.obj, resolution = c(0.1, 0.3, 0.5, 0.7, 1))
DimPlot(nsclc.seurat.obj, group.by = "RNA_snn_res.0.1", label = T)
DimPlot(nsclc.seurat.obj, group.by = "RNA_snn_res.0.3", label = T)
DimPlot(nsclc.seurat.obj, group.by = "RNA_snn_res.0.5", label = T)
DimPlot(nsclc.seurat.obj, group.by = "RNA_snn_res.0.7", label = T)
DimPlot(nsclc.seurat.obj, group.by = "RNA_snn_res.1", label = T)


Idents(nsclc.seurat.obj)
Idents(nsclc.seurat.obj) <- "RNA_snn_res.0.1"

nsclc.seurat.obj <- RunUMAP(nsclc.seurat.obj, dims = 1:15)
DimPlot(nsclc.seurat.obj, reduction = "umap")






