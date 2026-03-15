rm(list = ls())
gc()
setwd("D:/OneDrive - International Campus, Zhejiang University/文档/李牧轩/Zhejiang University/大二春夏学期/挑战杯/数据分析2.0")

library(tidyverse)
library(data.table)
library(DESeq2)

exprSet <- fread("Data/Counts.txt.gz", data.table = FALSE)
rownames(exprSet) <- exprSet$V1
exprSet$V1 <- NULL
rownames(exprSet) <- rownames(exprSet) %>% str_remove("\\.\\d*$")
exprSet <- as.matrix(exprSet)
mode(exprSet) <- "integer"

meta <- read.csv("Data/meta.csv")
rownames(meta) <- meta$sample_id
meta <- meta[colnames(exprSet), ]
meta$cancer_status <- factor(meta$cancer_status)
all(colnames(exprSet) == rownames(meta)) # TRUE
dds <- DESeqDataSetFromMatrix(
  countData = exprSet,
  colData   = meta,
  design    = ~ cancer_status
)

vsd <- vst(dds, blind = TRUE)
vst_mat <- assay(vsd)
fwrite(as.data.frame(vst_mat), "Data/VST.txt", sep = ',',
       row.names = TRUE, col.names = TRUE)

count_to_vst <- function(exprSet, colData = NULL, design = NULL,
                         dir.save = NULL) {
  tryCatch({
    exprSet <- as.matrix(exprSet)
    mode(exprSet) <- "integer"
  })
  if(!is.null(colData)) {
    if(ncol(exprSet) != nrow(colData)) {
      stop("Error: column number of expression
           do not match with row number of colData!")
    }
    colData <- colData[colnames(exprSet), ]
    if(any(colnames(exprSet) != rownames(colData))) {
      stop("Error in matching expression counts with meta data!")
    }
    dds <- DESeqDataSetFromMatrix(
      countData = exprSet,
      colData   = colData,
      design    = design
    )
    vsd <- vst(dds, blind = TRUE)
    vst_mat <- assay(vsd)
  } else {
    dds <- exprSet
    vst_mat <- vst(dds, blind = TRUE)
  }
  if(!is.null(dir.save)) {
    fwrite(as.data.frame(vst_mat), dir.save,
           row.names = TRUE, col.names = TRUE)
  }
}

exprSet <- fread("Data/TCGA-LAML.star_counts.tsv", data.table = F)
rownames(exprSet) <- exprSet$Ensembl_ID
exprSet$Ensembl_ID <- NULL
rownames(exprSet) <- rownames(exprSet) %>% str_remove("\\.\\d*$")

count_to_vst(exprSet, dir.save = "Data/TCGA-LAML-VST.txt")
