library(dplyr)
library(Seurat)
library(patchwork)
library(data.table)
library(tidyverse)
library(plotthis)
# Load the PBMC dataset
gex1 <- Read10X(data.dir = "/dartfs/rc/lab/M/McKennaLab/projects/Maryam/2024_04_07_cancer_single_cell/2025_06_24_mouse_and_human/lane1/166383_GEX1_Pool1/outs/filtered_feature_bc_matrix")
gex2 <- Read10X(data.dir = "/dartfs/rc/lab/M/McKennaLab/projects/Maryam/2024_04_07_cancer_single_cell/2025_06_24_mouse_and_human/lane2/166383_GEX1_Pool2/outs/filtered_feature_bc_matrix")
# Initialize the Seurat object with the raw (non-normalized data).
gex1o <- CreateSeuratObject(counts = gex1, project = "gex1", min.cells = 3, min.features = 200)
gex2o <- CreateSeuratObject(counts = gex2, project = "gex2", min.cells = 3, min.features = 200)


doublets_lane_1 = fread("/dartfs/rc/lab/M/McKennaLab/projects/Maryam/2024_04_07_cancer_single_cell/2025_06_24_mouse_and_human/lane1/round1_lane1_calls.txt",sep="\t")
doublets_lane_2 = fread("/dartfs/rc/lab/M/McKennaLab/projects/Maryam/2024_04_07_cancer_single_cell/2025_06_24_mouse_and_human/lane2/round1_lane2_calls.txt",sep="\t")
doublets_lane_1[, col := paste0(V1, "-1")]
doublets_lane_2[, col := paste0(V1, "-1")]
rownames(doublets_lane_1) = doublets_lane_1$col
rownames(doublets_lane_2) = doublets_lane_2$col

gex1o <- AddMetaData(gex1o, metadata = doublets_lane_1)
gex2o <- AddMetaData(gex2o, metadata = doublets_lane_2)

gex.combined <- merge(gex1o, y = gex2o, add.cell.ids = c("gex1o", "gex2o"), project = "gex")
gex.combined

gex.combined.subset <- subset(gex.combined, subset = x != "Doublet")
gex.combined.subset <- subset(gex.combined.subset, subset = x != "Negative")

gex.combined.subset[["percent.mt"]] <- PercentageFeatureSet(gex.combined.subset, pattern = "-MT-")
gex.combined.subset[["percent.mouse"]] <- PercentageFeatureSet(gex.combined.subset, pattern = "GRCm39-")
gex.combined.subset[["percent.human"]] <- PercentageFeatureSet(gex.combined.subset, pattern = "GRCh38-")

VlnPlot(gex.combined.subset, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.mouse","percent.human"), ncol = 3)


plot1 <- FeatureScatter(gex.combined.subset, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(gex.combined.subset, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

gex.combined.subset <- subset(gex.combined.subset, subset = nFeature_RNA > 200 & 
                                percent.mt < 15 &
                                percent.mouse < 5 &
                                percent.human > 95
                              )

gex.combined.subset <- NormalizeData(gex.combined.subset, normalization.method = "LogNormalize", scale.factor = 10000)


gex.combined.subset <- FindVariableFeatures(gex.combined.subset, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(gex.combined.subset), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(gex.combined.subset)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

all.genes <- rownames(gex.combined.subset)
gex.combined.subset <- ScaleData(gex.combined.subset, features = all.genes)


gex.combined.subset <- RunPCA(gex.combined.subset, features = VariableFeatures(object = gex.combined.subset))

# DimPlot(gex.combined.subset, reduction = "pca") + NoLegend()

print(gex.combined.subset[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(gex.combined.subset, dims = 1:2, reduction = "pca")


ElbowPlot(gex.combined.subset)

gex.combined.subset <- FindNeighbors(gex.combined.subset, dims = 1:20)
gex.combined.subset <- FindClusters(gex.combined.subset, resolution = 0.5)

gex.combined.subset <- RunUMAP(gex.combined.subset, dims = 1:10)
# DimPlot(gex.combined.subset, reduction = "umap")

options(future.globals.maxSize = 3e+09)
gex.combined.subset <- SCTransform(gex.combined.subset)
gex.combined.subset <- RunPCA(gex.combined.subset, npcs = 30, verbose = F)
options(future.globals.maxSize = 10e+09) # 1 GB


gex.combined.subset <- IntegrateLayers(
  object = gex.combined.subset,
  method = RPCAIntegration,
  normalization.method = "SCT",
  verbose = F,
  new.reduction = "integrated.cca"
)

# re-join layers after integration
gex.combined.subset[["RNA"]] <- JoinLayers(gex.combined.subset[["RNA"]])

gex.combined.subset <- FindNeighbors(gex.combined.subset, reduction = "integrated.cca", dims = 1:30)
gex.combined.subset <- FindClusters(gex.combined.subset, resolution = 1)
gex.combined.subset <- RunUMAP(gex.combined.subset, dims = 1:30, reduction = "integrated.cca")

#gex.combined.subset <- JoinLayers(gex.combined.subset)
gex.combined.subset <- PrepSCTFindMarkers(object = gex.combined.subset)

gex.combined.subset.markers <- FindAllMarkers(gex.combined.subset, only.pos = TRUE)
write.table(gex.combined.subset.markers,file="/dartfs/rc/lab/M/McKennaLab/projects/Maryam/2024_04_07_cancer_single_cell/2025_06_24_mouse_and_human/2026_02_01_gex.combined.subset.markers.txt",sep="\t")

top_genes = gex.combined.subset.markers %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 3)

###### Plot the UMAP using the dimplots package
gex.combined.subset.dp = as.data.frame(gex.combined.subset[["umap"]]@cell.embeddings)
gex.combined.subset.dp$cell_name = rownames(gex.combined.subset.dp)
clust_temp = as.data.frame(gex.combined.subset$seurat_clusters)


colnames(clust_temp) = "cluster"
clust_temp$cell_name = rownames(clust_temp)
gex.combined.subset.dp = left_join(gex.combined.subset.dp,clust_temp)
gex.combined.subset.dp$cell_name_filtered = sapply(gex.combined.subset.dp$cell_name,function(x){
  xx = strsplit(x,"-")[[1]][1]
  return(strsplit(xx,"_")[[1]][2])
  })


gex_umap = DimPlot(gex.combined.subset.dp, group_by = "cluster", theme = "theme_blank", highlight = TRUE, label = TRUE, pt_size=5, highlight_size=5, highlight_stroke = 1.5, label_size=10)
ggsave(gex_umap,file="/dartfs/rc/lab/M/McKennaLab/projects/Maryam/2024_04_07_cancer_single_cell/2025_06_24_mouse_and_human/2026_01_07_umap.pdf", width=10, height=10)


######################################################################################################
#
#
# Load up the lineage data
#
#
######################################################################################################
# load up the lineage recording side of things
lineage_recordings = fread("/dartfs/rc/lab/M/McKennaLab/projects/Maryam/2024_04_07_cancer_single_cell/lineage/2025_10_29_reanalysis/static_all_a375_cells.tsv")


# also load up our curated list of clones
currated_clone_list = fread("/dartfs/rc/lab/M/McKennaLab/projects/Maryam/2024_04_07_cancer_single_cell/lineage/2025_10_29_reanalysis/cluster_assignments_full_11_05_2025.txt")
colnames(currated_clone_list) = c("static_id","cluster","round","names","id")


rev_comp <- function(seq) {
  comp <- c(
    A = "T", T = "A",
    C = "G", G = "C",
    a = "t", t = "a",
    c = "g", g = "c"
  )
  
  paste0(rev(comp[strsplit(seq, "")[[1]]]), collapse = "")
}

# translate capture cell id sequences to polyA 10x cell ids
translation = fread("/dartfs/rc/lab/M/McKennaLab/resources/cellranger_versions/cellranger-9.0.1/lib/python/cellranger/barcodes/translation/3M-february-2018_NXT.txt.gz",header=F)
colnames(translation) = c("feature_code","cell_id")


lineage_recordings$e2_revcomp = sapply(lineage_recordings$e2,rev_comp)
lineage_recordings_merged = left_join(lineage_recordings, translation, by=c("e2_revcomp"="feature_code"))

lineage_recordings_filtered = lineage_recordings_merged[is.element(lineage_recordings_merged$cell_id,gex.combined.subset.dp$cell_name_filtered)]
lineage_recordings_merged_clone = left_join(lineage_recordings_filtered, currated_clone_list, by=c("e0"="static_id"))

# look at the UMI copies per cell
count_umis = as.data.table(table(lineage_recordings_merged_clone$e2,lineage_recordings_merged_clone$e0))
count_umis = count_umis[count_umis$N > 0,]
umi_counts_per_cell_plot = ggplot(count_umis) + geom_histogram(aes(N),binwidth=10) + theme_classic() + ylim(c(0,30000))
ggsave(umi_counts_per_cell_plot,file="/dartfs/rc/lab/M/McKennaLab/projects/Maryam/2024_04_07_cancer_single_cell/2025_06_24_mouse_and_human/2026_01_11_a375_single_cell_umi_counts_per_cell_plot.pdf",width=4,height=2)

count_umis_table = as.data.frame(table(count_umis$V1)) # count number of integrations across cells
count_umis_table_plot = ggplot(count_umis_table) + geom_histogram(aes(Freq)) + theme_classic() + xlab("Integrations per cell") + ylab("Count")
ggsave(count_umis_table_plot,file="/dartfs/rc/lab/M/McKennaLab/projects/Maryam/2024_04_07_cancer_single_cell/2025_06_24_mouse_and_human/2026_01_11_a375_count_umis_table_plot_plot.pdf",width=4,height=2)


# get the highest read count UMI per cell+UMI
lineage_recordings_merged_clone_dedup <- lineage_recordings_merged_clone %>%
  group_by(cell_static) %>%
  slice_max(order_by = reads_per_umi, n = 1, with_ties = FALSE) %>%
  ungroup()

lineage_recordings_merged_clone_dedup = lineage_recordings_merged_clone_dedup[!is.na(lineage_recordings_merged_clone_dedup$cluster),]

lineage_recordings_merged_clone_dedup_cell <- lineage_recordings_merged_clone_dedup %>%
  group_by(cell_id) %>%
  slice_max(order_by = reads_per_umi, n = 1, with_ties = FALSE) %>%
  ungroup()

lineage_recordings_merged_clone_dedup_cell_merged = left_join(gex.combined.subset.dp,lineage_recordings_merged_clone_dedup_cell,by=c("cell_name_filtered"="cell_id"))
lineage_recordings_merged_clone_dedup_cell_merged = lineage_recordings_merged_clone_dedup_cell_merged[!is.na(lineage_recordings_merged_clone_dedup_cell_merged$cigar_complexity),]
lineage_recordings_merged_clone_dedup_cell_merged$filtered_clones = lineage_recordings_merged_clone_dedup_cell_merged$cluster.y

keep_clusters = as.data.table(table(lineage_recordings_merged_clone_dedup_cell_merged$cluster.y))
keep_clusters_names = keep_clusters[keep_clusters$N > 10]$V1

lineage_recordings_merged_clone_dedup_cell_merged$filtered_clones[!lineage_recordings_merged_clone_dedup_cell_merged$cluster.y %in% keep_clusters_names] <- 0


for (cluster_id in keep_clusters_names) {
  lineage.idx <- which(lineage_recordings_merged_clone_dedup_cell_merged$cluster.y == cluster_id)
  lineage.plot = DimPlot(lineage_recordings_merged_clone_dedup_cell_merged, 
                         group_by = "filtered_clones", 
                         theme = "theme_blank", 
                         highlight = lineage.idx, 
                         pt_size=1, 
                         highlight_size=5, 
                         highlight_stroke = 1.5, 
                         label_size=10)
  
  ggsave(lineage.plot,file=paste("/dartfs/rc/lab/M/McKennaLab/projects/Maryam/2024_04_07_cancer_single_cell/2025_06_24_mouse_and_human/2026_01_07_",cluster_id,"_umap.pdf",sep=""))
}
lineage.idx <- which(lineage_recordings_merged_clone_dedup_cell_merged$cluster.y == 4 | lineage_recordings_merged_clone_dedup_cell_merged$cluster.y == 19)
lineage.plot = DimPlot(lineage_recordings_merged_clone_dedup_cell_merged, 
                       group_by = "filtered_clones", 
                       theme = "theme_blank", 
                       highlight = lineage.idx, 
                       pt_size=1, 
                       highlight_size=5, 
                       highlight_stroke = 1.5, 
                       label_size=10)
ggsave(lineage.plot,file=paste("/dartfs/rc/lab/M/McKennaLab/projects/Maryam/2024_04_07_cancer_single_cell/2025_06_24_mouse_and_human/2026_01_07_4_19_umap.pdf",sep=""))


lineage_recordings_merged_clone_dedup_cell_merged.cl4.idx <- which(lineage_recordings_merged_clone_dedup_cell_merged$cluster.y == 4)
lineage_recordings_merged_clone_dedup_cell_merged.cl19.idx <- which(lineage_recordings_merged_clone_dedup_cell_merged$cluster.y == 19)
lineage_recordings_merged_clone_dedup_cell_merged.cl1.idx <- which(lineage_recordings_merged_clone_dedup_cell_merged$cluster.y == 1)
lineage_recordings_merged_clone_dedup_cell_merged.cl6.idx <- which(lineage_recordings_merged_clone_dedup_cell_merged$cluster.y == 6)
lineage_recordings_merged_clone_dedup_cell_merged.cl16.idx <- which(lineage_recordings_merged_clone_dedup_cell_merged$cluster.y == 16)

lineage_recordings_merged_clone_dedup_cell_merged.cl4p19.idx = c(lineage_recordings_merged_clone_dedup_cell_merged.cl4.idx,lineage_recordings_merged_clone_dedup_cell_merged.cl19.idx)


# now add the lineage annotations back to the cell 
lineage_recordings_merged_clone_dedup_cell_merged_named = lineage_recordings_merged_clone_dedup_cell_merged
rownames(lineage_recordings_merged_clone_dedup_cell_merged_named) = lineage_recordings_merged_clone_dedup_cell_merged_named$cell_name


# now add the annotations to the seurat object
gex.combined.subset.lineage <- AddMetaData(
  object = gex.combined.subset,
  metadata = lineage_recordings_merged_clone_dedup_cell_merged_named
)

# do differential gene expression between clones
cluster.markers <- FindAllMarkers(gex.combined.subset.lineage, only.pos = TRUE)
cluster.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

cluster.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
plot = DoHeatmap(gex.combined.subset.lineage, features = top10$gene) + NoLegend()
ggplot2::ggsave(filename = "/dartfs/rc/lab/M/McKennaLab/projects/Maryam/2024_04_07_cancer_single_cell/2025_06_24_mouse_and_human/feature_by_cluster.pdf", width=15, height=15, plot = plot, device = "pdf")

# do differential gene expression between clones
clone.markers <- FindAllMarkers(gex.combined.subset.lineage, only.pos = TRUE, group.by = "cluster.y")
clone.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

# write the clonal differential expression
write.table(clone.markers,file="/dartfs/rc/lab/M/McKennaLab/projects/Maryam/2024_04_07_cancer_single_cell/2025_06_24_mouse_and_human/2026_02_01_clone_output_tables.txt",sep="\t")

clone.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 25) %>%
  ungroup() -> top10.clones
plot = DoHeatmap(gex.combined.subset.lineage, features = top10.clones$gene) + NoLegend()
ggplot2::ggsave(filename = "/dartfs/rc/lab/M/McKennaLab/projects/Maryam/2024_04_07_cancer_single_cell/2025_06_24_mouse_and_human/feature_by_cluster.pdf", width=15, height=15, plot = plot, device = "pdf")


obs_dr = gex.combined.subset.lineage[[]]

clone_4_names = obs_dr[obs_dr$filtered_clones == 4 & !is.na(obs_dr$filtered_clones),]
clone_19_names = obs_dr[obs_dr$filtered_clones == 19 & !is.na(obs_dr$filtered_clones),]
clone_1_names = obs_dr[obs_dr$filtered_clones == 1 & !is.na(obs_dr$filtered_clones),]

clone_4_names
write.table(clone_4_names,file="/dartfs/rc/lab/M/McKennaLab/projects/Maryam/2024_04_07_cancer_single_cell/2025_06_24_mouse_and_human/clone_4.txt",sep="\t",quote=F)
write.table(clone_19_names,file="/dartfs/rc/lab/M/McKennaLab/projects/Maryam/2024_04_07_cancer_single_cell/2025_06_24_mouse_and_human/clone_19.txt",sep="\t",quote=F)
write.table(clone_1_names,file="/dartfs/rc/lab/M/McKennaLab/projects/Maryam/2024_04_07_cancer_single_cell/2025_06_24_mouse_and_human/clone_1.txt",sep="\t",quote=F)

clone_table = data.table(table(obs_dr$filtered_clones))
colnames(clone_table) = c("og_clone_name","count")

clone_table_plot = ggplot(clone_table) + geom_bar(aes(x=as.factor(clone_name),y=count),stat="identity") + theme_classic()
ggsave(clone_table_plot,file="/dartfs/rc/lab/M/McKennaLab/projects/Maryam/2024_04_07_cancer_single_cell/2025_06_24_mouse_and_human/clone_table_plot.pdf",width=5,height=5 )

# save off the seurat object
saveRDS(object = gex.combined.subset, file = "/dartfs/rc/lab/M/McKennaLab/projects/Maryam/2024_04_07_cancer_single_cell/2025_06_24_mouse_and_human/2026_01_11_a375_single_cell_data.rds")

# run 4 vs 19 diff expression
four_vs_nineteen <- FindMarkers(gex.combined.subset.lineage, ident.1 = "4", ident.2 = "19", only.pos = FALSE, group.by = "cluster.y")

# write the clonal differential expression
write.table(four_vs_nineteen,file="/dartfs/rc/lab/M/McKennaLab/projects/Maryam/2024_04_07_cancer_single_cell/2025_06_24_mouse_and_human/2026_02_01_four_vs_nineteen.txt",sep="\t")


gex.combined.subset.lineage %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(pbmc, features = top10$gene) + NoLegend()