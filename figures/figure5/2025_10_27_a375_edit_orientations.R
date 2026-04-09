library(data.table)
library(tidyverse)
#library(seurat)

single_cell = fread("/dartfs/rc/lab/M/McKennaLab/projects/Maryam/2023_02_22_single_cell/combined_collpased_known.txt.gz")


cca = as.data.frame(table(single_cell$orientation_CCAGGAAGTACTCGAGTACTTCCTGG))
cca$sample = "CCAGGAAGTACTCGAGTACTTCCTGG"
cca$prop = cca$Freq / sum(cca$Freq)

cct = as.data.frame(table(single_cell$orientation_CCTGTCATCTTAGCTAAGATGACAGG))
cct$sample = "CCTGTCATCTTAGCTAAGATGACAGG"
cct$prop = cct$Freq / sum(cct$Freq)

proportions = rbind(cca,cct)

gg = ggplot(proportions) + geom_bar(aes(x=sample,y=prop,fill=Var1),stat="identity") + theme_classic() + 
  scale_fill_manual(values = c(
  "LEFT" = "#faa51a",
  "BOTH" = "#278c43",
  "RIGHT" = "#88cdea",
  "NONE" = "#999"
))

ggsave(gg,file="/dartfs/rc/lab/M/McKennaLab/projects/Maryam/2023_02_22_single_cell/proportion_single_cell_2025_10_27.pdf",width=3.5,height=6)


a375_single_cells_e6 = fread("/dartfs/rc/lab/M/McKennaLab/projects/Maryam/2024_04_07_cancer_single_cell/lineage/2025_10_29_reanalysis/E6_L2_S8_R2_combined_collpased_known.txt")
a375_single_cells_e5 = fread("/dartfs/rc/lab/M/McKennaLab/projects/Maryam/2024_04_07_cancer_single_cell/lineage/2025_10_29_reanalysis/E5_L3_S7_R2_combined_collpased_known.txt")
a375_single_cells_e2 = fread("/dartfs/rc/lab/M/McKennaLab/projects/Maryam/2024_04_07_cancer_single_cell/lineage/2025_10_29_reanalysis/E2_L2_S6_R2_combined_collpased_known.txt")
a375_single_cells_e1 = fread("/dartfs/rc/lab/M/McKennaLab/projects/Maryam/2024_04_07_cancer_single_cell/lineage/2025_10_29_reanalysis/E1_L1_S5_R2_combined_collpased_known.txt")

a375_single_cells_e6$sample = "a375_single_cells_e6"
a375_single_cells_e5$sample = "a375_single_cells_e5"
a375_single_cells_e2$sample = "a375_single_cells_e2"
a375_single_cells_e1$sample = "a375_single_cells_e1"

a375_single_cells = rbind(a375_single_cells_e6,a375_single_cells_e5,a375_single_cells_e2,a375_single_cells_e1)



static_a375_single_cells_e6 = fread("/dartfs/rc/lab/M/McKennaLab/projects/Maryam/2024_04_07_cancer_single_cell/lineage/2025_10_29_reanalysis/static_E6_L2_S8_R2_combined_collpased_known.txt")
static_a375_single_cells_e5 = fread("/dartfs/rc/lab/M/McKennaLab/projects/Maryam/2024_04_07_cancer_single_cell/lineage/2025_10_29_reanalysis/static_E5_L3_S7_R2_combined_collpased_known.txt")
static_a375_single_cells_e2 = fread("/dartfs/rc/lab/M/McKennaLab/projects/Maryam/2024_04_07_cancer_single_cell/lineage/2025_10_29_reanalysis/static_E2_L2_S6_R2_combined_collpased_known.txt")
static_a375_single_cells_e1 = fread("/dartfs/rc/lab/M/McKennaLab/projects/Maryam/2024_04_07_cancer_single_cell/lineage/2025_10_29_reanalysis/static_E1_L1_S5_R2_combined_collpased_known.txt")

static_a375_single_cells_e6$sample = "static_a375_single_cells_e6"
static_a375_single_cells_e5$sample = "static_a375_single_cells_e5"
static_a375_single_cells_e2$sample = "static_a375_single_cells_e2"
static_a375_single_cells_e1$sample = "static_a375_single_cells_e1"

static_a375_single_cells = rbind(static_a375_single_cells_e6, static_a375_single_cells_e5, static_a375_single_cells_e2, static_a375_single_cells_e1)


cca = as.data.frame(table(a375_single_cells$orientation_CCAGGAAGTACTCGAGTACTTCCTGG))
cca$sample = "CCAGGAAGTACTCGAGTACTTCCTGG"
cca$prop = cca$Freq / sum(cca$Freq)

cct = as.data.frame(table(a375_single_cells$orientation_CCTGTCATCTTAGCTAAGATGACAGG))
cct$sample = "CCTGTCATCTTAGCTAAGATGACAGG"
cct$prop = cct$Freq / sum(cct$Freq)

a375_single_cells_proportions = rbind(cca,cct)

gg = ggplot(a375_single_cells_proportions) + geom_bar(aes(x=sample,y=prop,fill=Var1),stat="identity") + theme_classic() + 
  scale_fill_manual(values = c(
    "LEFT" = "#faa51a",
    "BOTH" = "#278c43",
    "RIGHT" = "#88cdea",
    "NONE" = "#999"
  ))



a375_single_cells$cell_static = paste(a375_single_cells$e0,a375_single_cells$e2,sep="-")
static_a375_single_cells$cell_static = paste(static_a375_single_cells$e0,static_a375_single_cells$e2,sep="-")

a375_single_cells_dedup = a375_single_cells[!duplicated(a375_single_cells$cell_static),]
static_a375_single_cells_dedup = static_a375_single_cells[!duplicated(static_a375_single_cells$cell_static),]

a375_single_cells[a375_single_cells$cell_static == "TTCCACGTGTGATCTG-GACGCCAGTCGA",]
write.table(a375_single_cells,file="/dartfs/rc/lab/M/McKennaLab/projects/Maryam/2024_04_07_cancer_single_cell/lineage/2025_10_29_reanalysis/all_a375_cells.tsv",sep="\t",row.names=F,quote=F)
write.table(a375_single_cells_dedup,file="/dartfs/rc/lab/M/McKennaLab/projects/Maryam/2024_04_07_cancer_single_cell/lineage/2025_10_29_reanalysis/all_a375_cells_dedup.tsv",sep="\t",row.names=F,quote=F)

write.table(static_a375_single_cells,file="/dartfs/rc/lab/M/McKennaLab/projects/Maryam/2024_04_07_cancer_single_cell/lineage/2025_10_29_reanalysis/static_all_a375_cells.tsv",sep="\t",row.names=F,quote=F)
write.table(static_a375_single_cells_dedup,file="/dartfs/rc/lab/M/McKennaLab/projects/Maryam/2024_04_07_cancer_single_cell/lineage/2025_10_29_reanalysis/static_all_a375_cells_dedup.tsv",sep="\t",row.names=F,quote=F)

all_tags = as.data.frame(table(a375_single_cells_dedup$e2))
static_all_tags = as.data.frame(table(static_a375_single_cells_dedup$e0))
static_tag_merge = full_join(all_tags,static_all_tags,by=c("Var1"="Var1"))
static_tag_merge[is.na(static_tag_merge)] = 0
static_tag_merge_melt = reshape2::melt(static_tag_merge)

top_hits = as.data.frame(table(static_a375_single_cells_dedup$e0)[table(static_a375_single_cells_dedup$e0) > 10])
static_a375_single_cells_dedup_filter = static_a375_single_cells_dedup[is.element(static_a375_single_cells_dedup$e0,top_hits$Var1),]

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


static_a375_single_cells_dedup_filter = static_a375_single_cells_dedup[is.element(static_a375_single_cells_dedup$e0,top_hits$Var1),]
static_a375_single_cells_dedup_filter$e2_rc = sapply(static_a375_single_cells_dedup_filter$e2,revcomp)
static_a375_single_cells_dedup_filter$e2_r = sapply(static_a375_single_cells_dedup_filter$e2,rev)

static_a375_single_cells_dedup_filter_post = static_a375_single_cells_dedup_filter[is.element(static_a375_single_cells_dedup_filter$e2_rc,full_dups$V1.y),]

X <- table(static_a375_single_cells_dedup_filter$e2, static_a375_single_cells_dedup_filter$e0)

# Optional: prune empty/near-empty features
X <- X[, Matrix::colSums(X) > 5]
row_norm <- sqrt(Matrix::rowSums(X^2)); row_norm[row_norm == 0] <- 1
Xn <- X / row_norm

pc <- irlba::prcomp_irlba(Xn, n = 50, center = FALSE, scale. = FALSE)
Z  <- pc$x

# --- Path A: HDBSCAN
hdb <- hdbscan(Z, minPts = 30)
labels_hdb <- hdb$cluster   # 0 = noise

# --- Path B: kNN graph + Louvain
nn <- uwot::nearest_neighbors(Z, k = 30, metric = "cosine", method = "hnsw")
idx <- nn$idx
n <- nrow(Z); edges <- cbind(rep(1:n, each = ncol(idx)), as.vector(idx))
g <- simplify(graph_from_edgelist(as.matrix(edges), FALSE))
lab_louvain <- membership(cluster_louvain(g))




qplot(table(a375_single_cells_dedup$e2)[table(a375_single_cells_dedup$e2) > 100])


library(reticulate)
library(cluster)
library(dendextend)
library(ggplot2)

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


cluster_assignments = as.data.frame(cutree(res$tree_col, h = 0.95))
colnames(cluster_assignments) = "cluster"
cluster_assignments$cluster = as.factor(cluster_assignments$cluster)

cluster_assignments_row = as.data.frame(cutree(res$tree_row, k=40))
colnames(cluster_assignments_row) = "cluster"
cluster_assignments_row$cluster = as.factor(cluster_assignments_row$cluster)

# cluster_assignments$row.names = rownames(cluster_assignments)

res2 = pheatmap(
  X, 
  clustering_distance_rows = "binary",
  clustering_distance_cols = "binary",
  clustering_method = "average",
  show_rownames = FALSE,
  show_colnames = TRUE,
  color = colorRampPalette(c("white", "darkblue"))(100),
  annotation_names_row = TRUE,
  annotation_names_col = TRUE,
  annotation_col = cluster_assignments,
  annotation_row = cluster_assignments_row,
  
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


cluster_assignments2 = as.data.frame(cutree(res3$tree_col, k = 15))
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

cluster_assignments2$cluster_reassigned = cluster_assignments2$cluster
cluster_assignments2$cluster_reassigned[cluster_assignments2$cluster == 3] = 2
cluster_assignments2$cluster_reassigned[cluster_assignments2$cluster == 3] = 2