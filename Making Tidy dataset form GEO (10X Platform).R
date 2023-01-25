setwd("~/desktop/Kharazmi")

barcode_names <- paste0("~/desktop/Kharazmi/GSE174461_RAW/", list.files("~/desktop/Kharazmi/GSE174461_RAW", "_barcodes"))
feature       <- paste0("~/desktop/Kharazmi/GSE174461_RAW/", list.files("~/desktop/Kharazmi/GSE174461_RAW", "_features"))
mtx_names     <- paste0("~/desktop/Kharazmi/GSE174461_RAW/", list.files("~/desktop/Kharazmi/GSE174461_RAW", "_matrix.mtx"))
GSMs          <- sub("(^.*)(GSM\\d*)(.*$)", "\\2", mtx_names) # constructing three section and choose the second one

data_paths <- data.frame(barcode = barcode_names, feature = feature, mtx = mtx_names, GSMs = GSMs)

for (i in 1:nrow(data_paths)){
  
  barcods  <- read.table(data_paths[i,1])
  barcods  <- data.frame(V1 = sub("-1", "", barcods$V1))
  barcods  <- data.frame(V1 = sub(".1", "", barcods$V1)) # some of them have `.1` not `-1`
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
