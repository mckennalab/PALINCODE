library(PATH)
library(expm)
library(ggplot2)
library(ggtree)
library(ggtreeExtra)
library("SeuratDisk")
library(tidytree)
library(Matrix)
library(patchwork)
library(parallel)
library(msigdbr)
library(magrittr)
library(tidyverse)
library(fgsea)
library(data.table)
library(ape)
library(ggtreeExtra)
library(ggstar)
library(ggplot2)
library(ggtree)
library(treeio)
library(ggnewscale)
library(Seurat)
library(microViz)

# load global data needed for all trees
gex.combined.subset = readRDS(file = "/dartfs/rc/lab/M/McKennaLab/projects/Maryam/2024_04_07_cancer_single_cell/2025_06_24_mouse_and_human/2026_01_11_a375_single_cell_data.rds")
# translate capture cell id sequences to polyA 10x cell ids
translation = fread("/dartfs/rc/lab/M/McKennaLab/resources/cellranger_versions/cellranger-9.0.1/lib/python/cellranger/barcodes/translation/3M-february-2018_NXT.txt.gz",header=F)
colnames(translation) = c("cell_id","feature_code")

reverse_complement <- function(seq) {
  comp <- c(
    A = "T", C = "G", G = "C", T = "A",
    a = "t", c = "g", g = "c", t = "a"
  )
  paste0(rev(comp[strsplit(seq, "")[[1]]]), collapse = "")
}
reverse_nocomplement <- function(seq) {
  comp <- c(
    A = "T", C = "G", G = "C", T = "A",
    a = "t", c = "g", g = "c", t = "a"
  )
  paste0(rev(strsplit(seq, "")[[1]]), collapse = "")
}

gex.tbl = gex.combined.subset[[]]
gex.tbl$full_cell_id = rownames(gex.tbl)

# Generate 20 distinct colors
n_colors <- 14
my_colors <- scales::hue_pal()(n_colors) # Default ggplot2 hue palette


#######################################################################################################################
#
# Analyze clones 4 and 19 (now 1 and 2 on the plot) together
#
#######################################################################################################################
clones_to_lookup = c("4","19", "5", "6", "7", "8", "14", "16", "17", "20")
for (clone_number in clones_to_lookup) {
  print(clone_number)
  tryCatch({
  clone_number_tree <- read.tree(paste("/dartfs/rc/lab/M/McKennaLab/projects/Maryam/2024_04_07_cancer_single_cell/2025_06_24_mouse_and_human/all_2026_01_29_trees/cluster_",clone_number,"_test_phylip.mix.treefile",sep=""))
  clone_number_mapping <- fread(paste("/dartfs/rc/lab/M/McKennaLab/projects/Maryam/2024_04_07_cancer_single_cell/2025_06_24_mouse_and_human/all_2026_01_29_trees/cluster_",clone_number,"_test_phylip_mapping.txt",sep=""))
  
  # map the tip labels to the cell IDs
  lst <- setNames(as.list(clone_number_mapping$cell), clone_number_mapping$id)
  lookup <- setNames(clone_number_mapping$cell,clone_number_mapping$id)
  clone_number_tree$tip.label <- lookup[clone_number_tree$tip.label]
  
  # tree is a "phylo" from ape
  leaf_order <- as.data.frame(clone_number_tree$tip.label)
  colnames(leaf_order) = c("cell")
  
  leaf_order$rc = sapply(leaf_order$cell,reverse_complement)
  leaf_order$r = sapply(leaf_order$cell,reverse_nocomplement)
  
  #leaf_order_joined = left_join(leaf_order, translation, by=c("cell"="feature_code"))
  #leaf_order_joined2 = left_join(leaf_order_joined, gex.tbl, by=c("cell"="V1"))
  leaf_order_joined = left_join(leaf_order, translation, by=c("rc"="cell_id"))
  leaf_order_joined = left_join(leaf_order_joined, gex.tbl, by=c("feature_code"="V1"))
  
  leaf_order_joined2 = leaf_order_joined2[complete.cases(leaf_order_joined2),]
  
  p <- ggtree(clone_number_tree, layout="circular", branch.length="none")
  p2 <- p + 
    new_scale_fill() + 
    geom_fruit(
      data=leaf_order_joined2,
      geom=geom_point,
      mapping=aes(y=cell, fill=seurat_clusters), size=3, shape = 21, stroke=.25, col="black",
      offset=0.01,   # The distance between external layers, default is 0.03 times of x range of tree.
      pwidth=0.01 # width of the external layer, default is 0.2 times of x range of tree.
    ) + scale_fill_manual(values = my_colors)
  p2
  
  ggsave(p2,file=paste("/dartfs/rc/lab/M/McKennaLab/projects/Maryam/2024_04_07_cancer_single_cell/2025_06_24_mouse_and_human/clone_tree_clades_",clone_number,".pdf",sep=""))
  }, error = function(msg){
    message(paste("Error for list member:", clone_number))
  })
}




  print(paste("clone",clone_number))
  lst <- setNames(as.list(clone_number_mapping$cell), clone_number_mapping$id)
  lookup <- setNames(clone_number_mapping$cell,clone_number_mapping$id)
  clone_number_tree$tip.label <- lookup[clone_number_tree$tip.label]
  print(clone_number_tree)
  # tree distance calculation
  Winv <- inv_tree_dist(clone_number_tree, node = TRUE, norm = FALSE)
  
  # tree is a "phylo" from ape
  leaf_order <- as.data.frame(clone_number_tree$tip.label)
  colnames(leaf_order) = c("cell")
  
  
  leaf_order$rc = sapply(leaf_order$cell,reverse_complement)
  leaf_order$r = sapply(leaf_order$cell,reverse_nocomplement)
  leaf_order_joined = left_join(leaf_order, translation, by=c("cell"="cell_id"))
  gex.tbl = gex.combined.subset[[]]
  gex.tbl$full_cell_id = rownames(gex.tbl)
  leaf_order_joined = left_join(leaf_order_joined, gex.tbl, by=c("cell"="V1"))
  leaf_order_joined = leaf_order_joined[complete.cases(leaf_order_joined),]
  print(leaf_order_joined)
  
  # reorder + subset cells in the Seurat object
  seurat_obj <- gex.combined.subset
  cell_id_list = rownames(seurat_obj[[]][is.element(gex.combined.subset[[]]$V1,leaf_order_joined$cell),])
  
  # keep only cells present in both (and preserve tree order)
  leaf_order <- leaf_order[is.element(leaf_order_joined$cell,gex.combined.subset[[]]$V1),]
  
  seurat_hvg <- subset(seurat_obj, features = VariableFeatures(seurat_obj))
  
  seurat_obj_subset = seurat_hvg[,leaf_order_joined$full_cell_id]
  
  mat <- t(seurat_obj_subset[["RNA"]]$scale.data)
  
  # non_zero_counts = rowSums(GetAssayData(seurat_obj_subset, assay = "RNA", slot = "counts")) > 25
  
  modxcor <- xcor(mat, Winv)
  
  # Process xcor output for plotting and visualization.
  Idf <- reshape2::melt(modxcor$phy_cor, 
                        value.name = "I")
  Zdf <- reshape2::melt(modxcor$Z.score, 
                        value.name = "Z")
  
  df <- full_join(Idf, Zdf, by=c("Var1", "Var2"))
  
  df <- df %>% mutate(Var1 = as.factor(Var1), 
                      Var2 = as.factor(Var2))
  
  # Phylogenetic auto-correlation bar plot.
  df = df[complete.cases(df),]
  df_sorted = df[order(df$Z,decreasing = T),]
  # write.table(df_sorted,file="/dartfs/rc/lab/M/McKennaLab/projects/Maryam/2024_04_07_cancer_single_cell/2025_06_24_mouse_and_human/2026_01_29_clone_4_PATH.txt",sep="\t")
  
  # remap the tip labels again -- this time to the cell data (not the capture sequence)
  
  lookup_translation <- setNames(leaf_order_joined$cell_id,leaf_order_joined$cell)
  clone_number_tree$tip.label <- lookup_translation[clone_number_tree$tip.label]
  

  df_sorted_matched = df_sorted %>% filter(Var1 == Var2)
  dat2 = FetchData(seurat_obj_subset, vars = head(df_sorted_matched,n=400)$Var1)
  #dat2 = FetchData(seurat_obj_subset, vars = rownames(de)[1:2])
  dat2 = dat2[,apply(dat2, 2, median) > 0]
  colnames_dat2 = colnames(dat2)
  
  dat2$ID_full = rownames(dat2)
  dat2$ID = sapply(dat2$ID_full,function(x){
    tt = strsplit(x,"_")[[1]][2]
    return(strsplit(tt,"-")[[1]][1])
  })
  
  dat2_melt = reshape2::melt(dat2,id.vars=c("ID"),measure.vars=colnames_dat2)

  p <- ggtree(clone_number_tree, layout="circular", branch.length="none")
  p2 <- p + 
    new_scale_fill() + 
    geom_fruit(
      data=dat2_melt,
      geom=geom_tile,
      mapping=aes(y=ID, x=variable, fill=value),
      offset=0.08,   # The distance between external layers, default is 0.03 times of x range of tree.
      pwidth=0.1 # width of the external layer, default is 0.2 times of x range of tree.
    ) + scale_fill_gradientn(
      colours = c("blue", "dodgerblue", "white", "orange", "red"),
      rescaler = ~ scales::rescale_mid(.x, mid = 0)
    )
  
  # p2
  ggsave(p2,file=paste("/dartfs/rc/lab/M/McKennaLab/projects/Maryam/2024_04_07_cancer_single_cell/2025_06_24_mouse_and_human/clone_association_",clone_number,".pdf",sep=""))
}