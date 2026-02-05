rm(list = ls())
options(future.globals.maxSize=160*1024^3)
library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)
library(copykat)
options(stringsAsFactors=FALSE)
load("/data5/zhangq/ThrombusST/data/scRNAV2/PLAGH_GSE264586_HRA004711_RCCTT_FirstAnnotation.rda")
setwd("/data5/zhangq/ThrombusST/data/scRNAV2/Tumor/copykat")
D301obj <- subset(RCC_SCT, subset = bcluster %in% c("Epi", "Tcell", "Monocytic") & nFeature_RNA >= 500)

listCountMtx <- list(); listNormCells <- list()
for(samplei in unique(as.character(D301obj$orig.ident)))
{
    tmp <- subset(D301obj, subset = orig.ident == samplei)
    listCountMtx[[samplei]] <- tmp@assays$RNA@counts
    listNormCells[[samplei]] <- rownames(tmp@meta.data)[which(tmp@meta.data$bcluster %in% c("Tcell", "Monocytic"))]
    print(samplei)
}

for(samplei in unique(as.character(D301obj$orig.ident)))
{
    #samplei <- "PT3T"
    exp.rawdata <- listCountMtx[[samplei]]
    norm_cell <- listNormCells[[samplei]]
    copykat.test <- copykat(rawmat=exp.rawdata, id.type="S", cell.line="no",  ngene.chr=5, win.size=25, 
                        norm.cell.names=norm_cell, KS.cut=0.1,  sam.name=samplei, distance="euclidean",  n.cores=4)
    save(copykat.test, file=paste0(samplei, "_Filter_copykat.rda"))
    print(samplei)
}

