library(data.table)
library(tidyverse)

library(reticulate)
library(cluster)
library(dendextend)
library(ggplot2)
library(pheatmap)
library(grDevices)

static_a375_single_cells_e6 = fread("/dartfs/rc/lab/M/McKennaLab/projects/Maryam/2024_04_07_cancer_single_cell/lineage/2025_10_29_reanalysis/static_E6_L2_S8_R2_combined_collpased_known.txt")
static_a375_single_cells_e5 = fread("/dartfs/rc/lab/M/McKennaLab/projects/Maryam/2024_04_07_cancer_single_cell/lineage/2025_10_29_reanalysis/static_E5_L3_S7_R2_combined_collpased_known.txt")
static_a375_single_cells_e2 = fread("/dartfs/rc/lab/M/McKennaLab/projects/Maryam/2024_04_07_cancer_single_cell/lineage/2025_10_29_reanalysis/static_E2_L2_S6_R2_combined_collpased_known.txt")
static_a375_single_cells_e1 = fread("/dartfs/rc/lab/M/McKennaLab/projects/Maryam/2024_04_07_cancer_single_cell/lineage/2025_10_29_reanalysis/static_E1_L1_S5_R2_combined_collpased_known.txt")

static_a375_single_cells_e6$sample = "static_a375_single_cells_e6"
static_a375_single_cells_e5$sample = "static_a375_single_cells_e5"
static_a375_single_cells_e2$sample = "static_a375_single_cells_e2"
static_a375_single_cells_e1$sample = "static_a375_single_cells_e1"

static_a375_single_cells = rbind(static_a375_single_cells_e6, static_a375_single_cells_e5, static_a375_single_cells_e2, static_a375_single_cells_e1)

static_a375_single_cells$cell_static = paste(static_a375_single_cells$e0,static_a375_single_cells$e2,sep="-")

static_a375_single_cells_dedup = static_a375_single_cells[!duplicated(static_a375_single_cells$cell_static),]

static_all_tags = as.data.frame(table(static_a375_single_cells_dedup$e0))

dups_lane1 = fread("/dartfs/rc/lab/M/McKennaLab/projects/Maryam/2024_04_07_cancer_single_cell/2025_06_24_mouse_and_human/lane1/round1_lane1_calls.txt")
dups_lane2 = fread("/dartfs/rc/lab/M/McKennaLab/projects/Maryam/2024_04_07_cancer_single_cell/2025_06_24_mouse_and_human/lane2/round1_lane2_calls.txt")
full_dups = rbind(dups_lane1,dups_lane2)
full_dups$duplicated = FALSE
full_dups$duplicated[full_dups$x == "Doublet" | full_dups$x == "Negative"] = TRUE
full_dups_undup = full_dups[!full_dups$duplicated,]

translation = fread("/dartfs/rc/lab/M/McKennaLab/resources/cellranger_versions/cellranger-9.0.1/lib/python/cellranger/barcodes/translation/3M-february-2018_NXT.txt.gz",header=F)
full_dups_undup_cell_ids = left_join(full_dups_undup,translation,by=c("V1"="V2"))

revcomp <- function(s) {
  s <- toupper(s)
  s <- strsplit(s, NULL)[[1]]
  s <- rev(s)
  s <- chartr("ACGT", "TGCA", s)
  paste(s, collapse = "")
}
revit <- function(s) {
  s <- toupper(s)
  s <- strsplit(s, NULL)[[1]]
  s <- rev(s)
  # s <- chartr("ACGT", "TGCA", s)
  paste(s, collapse = "")
}



static_a375_single_cells_dedup$e2_rc = sapply(static_a375_single_cells_dedup$e2,revcomp)
static_a375_single_cells_dedup$e2_r = sapply(static_a375_single_cells_dedup$e2,revit)

static_a375_single_cells_dedup = static_a375_single_cells_dedup[is.element(static_a375_single_cells_dedup$e2_rc,full_dups_undup_cell_ids$V1.y),]
top_hits = as.data.frame(table(static_a375_single_cells_dedup$e0)[table(static_a375_single_cells_dedup$e0) > 5])
static_a375_single_cells_dedup = static_a375_single_cells_dedup[is.element(static_a375_single_cells_dedup$e0,top_hits$Var1),]

X <- table(static_a375_single_cells_dedup$e2, static_a375_single_cells_dedup$e0)

# Compute Jaccard distance
jaccard_dist <- dist(X, method = "binary")

# Hierarchical clustering with average linkage
hclust_result <- hclust(jaccard_dist, method = "average")

res = pheatmap(
  X, 
  clustering_distance_rows = "binary",
  clustering_distance_cols = "binary",
  clustering_method = "average",
  show_rownames = FALSE,
  show_colnames = FALSE,
  color = colorRampPalette(c("white", "darkblue"))(100),
  annotation_names_row = TRUE,
  annotation_names_col = TRUE,
  
)
res

cluster_assignments = as.data.frame(cutree(res$tree_col, k=40))
colnames(cluster_assignments) = "cluster"
cluster_assignments$cluster = as.factor(cluster_assignments$cluster)

cluster_assignments_row = as.data.frame(cutree(res$tree_row, k=40))
colnames(cluster_assignments_row) = "cluster"
cluster_assignments_row$cluster = as.factor(cluster_assignments_row$cluster)

# cluster_assignments$row.names = rownames(cluster_assignments)
# https://stackoverflow.com/questions/41628450/r-pheatmap-change-annotation-colors-and-prevent-graphics-window-from-popping-up
newCols <- colorRampPalette(grDevices::rainbow(length(unique(cluster_assignments$cluster))))
mycolors <- newCols(length(unique(cluster_assignments$cluster)))
names(mycolors) <- unique(cluster_assignments$clsuter)
mycolors <- list(category = mycolors)

res2 = pheatmap(
  X, 
  clustering_distance_rows = "binary",
  clustering_distance_cols = "binary",
  clustering_method = "average",
  show_rownames = FALSE,
  show_colnames = TRUE,
  color = colorRampPalette(c("white", "darkblue"))(100),
  annotation_col = cluster_assignments,
  annotation_colors = mycolors
)
res2


# now we want to recluster 'cluster 1' -- just extract those static codes and the cells they map too
just_cluster1 = cluster_assignments
just_cluster1$keep = FALSE
just_cluster1$keep[just_cluster1$cluster == 1] = TRUE

res3 = pheatmap(
  X[,just_cluster1$keep], 
  clustering_distance_rows = "binary",
  clustering_distance_cols = "binary",
  clustering_method = "average",
  show_rownames = FALSE,
  show_colnames = TRUE,
  color = colorRampPalette(c("white", "darkblue"))(100),
  annotation_names_row = TRUE,
  annotation_names_col = TRUE,
)
res3


cluster_assignments2 = as.data.frame(cutree(res3$tree_col, k = 4))
colnames(cluster_assignments2) = "cluster"
cluster_assignments2$cluster = as.factor(cluster_assignments2$cluster)


res4 = pheatmap(
  X[,just_cluster1$keep], 
  clustering_distance_rows = "binary",
  clustering_distance_cols = "binary",
  clustering_method = "average",
  show_rownames = FALSE,
  show_colnames = TRUE,
  color = colorRampPalette(c("white", "darkblue"))(100),
  annotation_names_row = TRUE,
  annotation_names_col = TRUE,
  annotation_col = cluster_assignments2,
)
res4

cluster_assignments$round = 1
cluster_assignments2$round = 2
cluster_assignments$names = rownames(cluster_assignments)
cluster_assignments2$names = rownames(cluster_assignments2)

cluster_assignments_full = rbind(cluster_assignments, cluster_assignments2)
cluster_assignments_full$id = paste(cluster_assignments_full$cluster,cluster_assignments_full$round,sep=".")

write.table(cluster_assignments_full,file="/dartfs/rc/lab/M/McKennaLab/projects/Maryam/2024_04_07_cancer_single_cell/lineage/2025_10_29_reanalysis/cluster_assignments_full_11_05_2025.txt",quote=F,sep="\t",row.names=T)


