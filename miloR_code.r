options(future.globals.maxSize=160*1024^3)
library(stringr)
library(Seurat)
library(miloR)
library(scater)
library(scran)
library(dplyr)
library(patchwork)
library(qs)
library(BiocParallel)
options(stringsAsFactors=FALSE)

groupi <- "antiTIM3"

tmp <- subset(RCC_SCT, subset = group %in% c("Control", groupi))
tmp$group <- as.character(tmp$group)
tmp$group[which(tmp$group %in% "Control")] <- "con"
tmp$group[which(tmp$group %in% groupi)] <- "treat"

set.seed(123)
subsampled_cells <- tmp@meta.data %>%
rownames_to_column("cell_id") %>%
group_by(orig.ident) %>%
sample_n(5000) %>%
pull(cell_id)
tmp <- subset(tmp, cells = subsampled_cells)

tmp$bcode <- factor(tmp$bcode, levels=c("Normalepi", "Tumorepi","Fib","Endo","Mono","Macro","DC","Neutrophil","CD4T","CD8T","NK","Bcell"))
scRNA_sce <- as.SingleCellExperiment(tmp)

# 获取 Seurat 对象中的所有降维名称
reduction_names <- names(tmp@reductions)

# 将每个降维结果添加到 SCE 对象中
for (name in reduction_names) {
  reducedDim(scRNA_sce, toupper(name)) <- Embeddings(tmp, name)
  print(name)
  print(toupper(name))
}

# 构建miloR对象
scRNA_pre <- Milo(scRNA_sce)


##构建K最邻图谱
# d值与降维时选定的pca值一致
upcs <- 50
# k值与FindNeighbors时选定的值一致
reducedDims(scRNA_pre)
# List of length 4
# names(4): PCA HARMONY UMAP TSNE
scRNA_pre <- buildGraph(scRNA_pre, k = 15, d = upcs,reduced.dim = "PCA")

##在 KNN 图上定义具有代表性的邻域

# 定义数据点的值(Prop)
# 对于小于3万细胞的数据集,Prop设为0.1
# 对于超过3万细胞的数据集,prop设为0.05。
scRNA_pre <- makeNhoods(scRNA_pre, prop = 0.05, k = 15, d = upcs, 
                        refined = TRUE, reduced_dims = "PCA")
plotNhoodSizeHist(scRNA_pre)

##计算neighbourhoods中的细胞
scRNA_pre <- countCells(scRNA_pre, 
                    meta.data = as.data.frame(colData(scRNA_pre)),
                    sample="orig.ident")

# colData(scRNA_pre) 相当于代表了metadata的信息
# sample代表了样本数
head(nhoodCounts(scRNA_pre))

# 实验设计定义
scRNA_design <- data.frame(colData(scRNA_pre))[,c("orig.ident","group")]

## 将批次信息转换为整数
scRNA_design$batch <- as.numeric(as.factor(scRNA_design$orig.ident))
scRNA_design <- distinct(scRNA_design)
rownames(scRNA_design) <- scRNA_design$orig.ident

scRNA_design

##计算邻里连通性
# MiloR使用了由cydar引入的空间FDR校正,这一步可以校正邻里之间重叠的p值
scRNA_pre <- calcNhoodDistance(scRNA_pre, d=upcs, reduced.dim = "PCA")

##检测结果
# 把批次和分组信息加到design中去
da_results <- testNhoods(scRNA_pre, design = ~ batch + group, design.df = scRNA_design)

head(da_results)
da_results %>%
  arrange(SpatialFDR) %>%
  head() 

#可视化
scRNA_pre <- buildNhoodGraph(scRNA_pre)

## Plot single-cell UMAP
umap_pl <- plotReducedDim(scRNA_pre, dimred = "UMAP", 
                          colour_by="group", text_by = "bcode", 
                          text_size = 3, point_size=0.5) + guides(fill="none")+theme(aspect.ratio=1)

## Plot neighbourhood graph
nh_graph_pl <- plotNhoodGraphDA(scRNA_pre, da_results, layout="UMAP",alpha=1) + 
scale_fill_gradient2(low="blue", mid="white", high="red", limits=c(min(da_results$logFC)+0.05,max(da_results$logFC)+0.05))+
#scale_fill_gradient2(low="blue", mid="white", high="red", limits=c(-6,6))+
theme(aspect.ratio=1)  #alpha默认0.1
