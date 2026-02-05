################################GSEA##########################
#################################noICB vs ICBR#########################
rm(list = ls())
library(Seurat)
library(ggplot2)
library(monocle3)

load("D301TNK_EndAnnotation.rda")

expression_matrix <- D301TNK@assays$SCT@data
cell_metadata <- D301TNK@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(D301TNK), row.names = row.names(D301TNK))

cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

cds <- preprocess_cds(cds, num_dim = 8)
cds <- reduce_dimension(cds)

pdf("1.reduce_dimension.pdf", height = 4, width = 4)
plot_cells(cds)
dev.off()

pdf("1.remove.batch.effects.raw.pdf", height = 4, width = 5)
plot_cells(cds, color_cells_by="orig.ident", label_cell_groups=FALSE)
dev.off()

cds <- align_cds(cds, num_dim = 8, alignment_group = "orig.ident")
cds <- reduce_dimension(cds)
pdf("1.remove.batch.effects.pdf", height = 4, width = 5)
plot_cells(cds, color_cells_by="orig.ident", label_cell_groups=FALSE)
dev.off()

pdf("1.remove.batch.effects.CNSgroup.pdf", height = 4, width = 5)
plot_cells(cds, color_cells_by="CNSgroup", label_cell_groups=FALSE)
dev.off()

#Group cells into clusters
cds <- cluster_cells(cds, resolution=1e-5)
pdf("1.clusters.pdf", height = 4, width = 5)
plot_cells(cds)
dev.off()

#plot_cells(cds)
pdf("1.partition.pdf", height = 4, width = 5)
plot_cells(cds, color_cells_by="partition", group_cells_by="partition")
dev.off()

#Find marker genes expressed by each cluster
marker_test_res <- top_markers(cds, group_cells_by="partition", 
                               reference_cells=1000, cores=8)
top_specific_markers <- marker_test_res %>%
                            filter(fraction_expressing >= 0.10) %>%
                            group_by(cell_group) %>%
                            top_n(3, pseudo_R2)

top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))

pdf("1.genes_by_group.pdf", height = 8, width = 4)
plot_genes_by_group(cds,
                    top_specific_marker_ids,
                    group_cells_by="partition",
                    ordering_type="cluster_row_col",
                    max.size=3)
dev.off()

cds <- learn_graph(cds)
pdf("1.learn_graph.pdf", height = 4, width = 4)
plot_cells(cds,
           color_cells_by = "partition",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)
dev.off()

# a helper function to identify the root principal points:
get_earliest_principal_node <- function(cds, time_bin="Naive"){
  cell_ids <- which(colData(cds)[, "CNSgroup"] == time_bin)
  
  closest_vertex <-
  cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
  igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
  (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}

cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))

pdf("1.pseudotime.pdf", height = 4, width = 5)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)
dev.off()