setwd("~/desktop/Kharazmi")

barcode_names <- paste0("~/desktop/Kharazmi/GSE174461_RAW/", list.files("~/desktop/Kharazmi/GSE174461_RAW", "_barcodes"))
feature       <- paste0("~/desktop/Kharazmi/GSE174461_RAW/", list.files("~/desktop/Kharazmi/GSE174461_RAW", "_features"))
mtx_names     <- paste0("~/desktop/Kharazmi/GSE174461_RAW/", list.files("~/desktop/Kharazmi/GSE174461_RAW", "_matrix.mtx"))
GSMs <- (c("GSM5311807", "GSM5311808", "GSM5311809", "GSM5311810", "GSM5311811", "GSM5311812"))

data_paths <- data.frame(barcode = barcode_names, feature = feature, mtx = mtx_names, GSMs = GSMs)

for (i in nrow(data_paths)){
  
  barcods  <- read.table(data_paths[i,1])
  barcods  <- data.frame(V1 = sub("-1", "", barcods$V1))
  features <- read.table(data_paths[i,2])
  features <- features[,-c(3,4)]
  
  mtx <- read.table(data_paths[i,3], skip = 3)
  
  count <- matrix(0, nrow(features), nrow(barcods))
  rownames(count) <- features$V1
  colnames(count) <- barcods$V1
  
  dim(features)
  max(mtx$V1)
  dim(barcods)
  max(mtx$V2)
  
  for(j in 1:nrow(mtx)){
    count[mtx[j,1], mtx[j,2]] = mtx[j,3]
  }
  
count_file_name <- data_paths[i,4]
write.table(count, file = paste0(count_file_name, ".csv"), sep = "\t")

}

nrow(data_paths)


