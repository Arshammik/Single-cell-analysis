setwd("~/desktop/Kharazmi")

barcode_names <- paste0("~/desktop/Kharazmi/GSE174461_RAW/", list.files("~/desktop/Kharazmi/GSE174461_RAW", "_barcodes"))
feature       <- paste0("~/desktop/Kharazmi/GSE174461_RAW/", list.files("~/desktop/Kharazmi/GSE174461_RAW", "_features"))
mtx_names     <- paste0("~/desktop/Kharazmi/GSE174461_RAW/", list.files("~/desktop/Kharazmi/GSE174461_RAW", "_matrix.mtx"))
GSMs          <- sub("(^.*)(GSM\\d*)(.*$)", "\\2", mtx_names)

data_paths <- data.frame(barcode = barcode_names, feature = feature, mtx = mtx_names, GSMs = GSMs)

for (i in 1:nrow(data_paths)){
  
  barcods  <- read.table(data_paths[i,1])
  barcods  <- data.frame(V1 = sub("-1", "", barcods$V1))
  features <- read.table(data_paths[i,2])
  features <- features[,-c(3,4)]
  gene_ENS <- paste0(features$V2, "—", features$V1)

  mtx <- read.table(data_paths[i,3], skip = 3)
  
  count <- matrix(0, nrow(features), nrow(barcods))
  rownames(count) <- gene_ENS
  colnames(count) <- barcods$V1
  
  
  for(j in 1:nrow(mtx)){
    count[mtx[j,1], mtx[j,2]] = mtx[j,3]
  }
  
  count_file_name <- data_paths[i,4]
  print(paste0("====> The file:",count_file_name," Done Processing!"))
  write.table(count, file = paste0("single_cell_count_mateix_",count_file_name, ".csv"), sep = "\t")

}


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
pbmc <-  CreateSeuratObject(counts = pbmc.data, project = "Breast Cancer", min.cells = 3, min.features = 200)
pbmc



# 2. QC
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics as a violin plot
pdf("GSE174461_vilon_plot_V.0.0.pdf")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()



barcods  <- read.table(data_paths[1,1])
barcods  <- data.frame(V1 = sub("-1", "", barcods$V1))
features <- read.table(data_paths[1,2])
features <- features[,-c(3,4)]

dim(features)
features <- distinct(features, V2, .keep_all = T)

mtx <- read.table(data_paths[1,3], skip = 3)

count <- matrix(0, nrow(features), nrow(barcods))
rownames(count) <- features$V2
colnames(count) <- barcods$V1





