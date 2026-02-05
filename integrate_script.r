rm(list = ls())
options(future.globals.maxSize=160*1024^3)
library(Seurat)
library(Matrix)
library(dplyr)
library(scDblFinder)
library(ggplot2)
library(SingleCellExperiment)
options(stringsAsFactors=FALSE)
samplesset = list.files("/data5/zhangq/ThrombusST/data/scRNA/scMATRIX/")

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

result.matrix = matrix(0,length(samplesset)[1],5)
i = 0
RCC.list = list()

for(sample_s in samplesset){
  #sample_s="Pre_9_GXL_PT"
  i = i + 1
  print(sample_s)
  datadir = paste("/data5/zhangq/ThrombusST/data/scRNA/scMATRIX/",sample_s,sep="")
  context_ <- list.files(datadir)
  {
   if("filtered_feature_bc_matrix" %in% context_){
        datadir = paste0(datadir, "/filtered_feature_bc_matrix")
        names(datadir) <- sample_s
        RCC.data.i <- Read10X(data.dir = datadir)
  }else{
        names(datadir) <- sample_s
        RCC.data.i <- Read10X(data.dir = datadir)
  }
  }

  result.matrix[i,1] = sample_s
  result.matrix[i,2] = dim(RCC.data.i)[2]
  RCC.seurat <- CreateSeuratObject(counts = RCC.data.i, min.cells = 3, min.features = 200, project = sample_s)
  RCC.seurat$orig.ident <- sample_s
  
  RCC.seurat[["percent.mito"]] <- PercentageFeatureSet(RCC.seurat, pattern = "^MT-")
  print(dim(RCC.seurat))
  HB.genes_total=c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
  HB_m=match(HB.genes_total,rownames(RCC.seurat@assays$RNA))
  HB.genes=rownames(RCC.seurat@assays$RNA)[HB_m]
  HB.genes=HB.genes[!is.na(HB.genes)]
  RCC.seurat[["percent.HB"]]=PercentageFeatureSet(RCC.seurat,features=HB.genes)
  head(RCC.seurat@meta.data)[,c(2,3,4,5)]

  rb.genes <- rownames(RCC.seurat)[grep("^RP[SL]",rownames(RCC.seurat))]
  C<-GetAssayData(object = RCC.seurat, slot = "counts")
  percent.ribo <- Matrix::colSums(C[rb.genes,])/Matrix::colSums(C)*100
  RCC.seurat <- AddMetaData(RCC.seurat, percent.ribo, col.name = "percent.ribo")

  result.matrix[i,3] = dim(RCC.seurat)[2]
  
  #RCC.seurat <- subset(RCC.seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mito < 50)
  RCC.seurat <- subset(RCC.seurat, subset = nFeature_RNA > 100 & nFeature_RNA < 6000 & percent.mito < 50)
  RCC.seurat <- subset(RCC.seurat, subset = nCount_RNA > 300 & nCount_RNA < 30000)

  rb.genes <- rownames(RCC.seurat)[grep("^RP[SL]",rownames(RCC.seurat))]
  mt.genes <- rownames(RCC.seurat)[grep("^MT-",rownames(RCC.seurat))]
  rbmt <- c(rb.genes, mt.genes,"MALAT1","NEAT1")
  filtergene <- setdiff(rownames(RCC.seurat), rbmt)
  RCC.seurat <- RCC.seurat[filtergene,]

  #Doublet Finder
  ##  Detect Doublets 
  set.seed(123) 
  sce <- as.SingleCellExperiment(RCC.seurat) 
  sce <- scDblFinder(sce, dbr=0.1) 
  table(sce$scDblFinder.class)
  scDblObj <- colData(sce)
  scDblObj <- as.data.frame(scDblObj)
  identical(rownames(scDblObj), rownames(RCC.seurat@meta.data)) #TRUE
  RCC.seurat@meta.data$scDblFinder.class <- scDblObj$scDblFinder.class
  RCC.seurat=subset(RCC.seurat, subset = scDblFinder.class == "singlet")

  dim(RCC.seurat)
  RCC.seurat <- NormalizeData(RCC.seurat, normalization.method = "LogNormalize", scale.factor = 10000)
  RCC.seurat <- CellCycleScoring(RCC.seurat, g2m.features=g2m.genes, s.features=s.genes)
  RCC.seurat <- FindVariableFeatures(RCC.seurat, selection.method = "vst", nfeatures = 2000)
  result.matrix[i,4] = dim(RCC.seurat)[2]
  
  RCC.seurat <- SCTransform(RCC.seurat,verbose = FALSE
                            ,vars.to.regress = c("percent.mito","percent.HB","percent.ribo","S.Score","G2M.Score"))
  RCC.seurat <- RunPCA(RCC.seurat, verbose = F, features = VariableFeatures(RCC.seurat),npcs = 50)
  DefaultAssay(RCC.seurat)
  print(i)

  result.matrix[i,5] = dim(RCC.seurat)[2]
  RCC.list[sample_s]=RCC.seurat
  
}

result.dataframe = as.data.frame(result.matrix)
colnames(result.dataframe) = c('sample','before','filter1',"filter2","drop_Doublet")

features <- SelectIntegrationFeatures(object.list = RCC.list, nfeatures = 2000)
RCC.list <- PrepSCTIntegration(object.list = RCC.list, anchor.features = features)
immune.anchors <- FindIntegrationAnchors(object.list = RCC.list, reduction = c("rpca"),dims = 1:50, normalization.method = "SCT",anchor.features = features)
RCC_SCT <- IntegrateData(anchorset = immune.anchors, dims = 1:50, normalization.method = "SCT")

RCC_SCT <- RunPCA(RCC_SCT, verbose = FALSE,npcs = 50)
RCC_SCT <- RunUMAP(RCC_SCT, reduction = "pca", dims = 1:50)
RCC_SCT <- FindNeighbors(RCC_SCT, reduction = "pca", dims = 1:50)
RCC_SCT <- FindClusters(RCC_SCT, resolution = 0.8)