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
library(caret)

# load global data needed for all trees
gex.combined.subset = readRDS(file = "/dartfs/rc/lab/M/McKennaLab/projects/Maryam/2024_04_07_cancer_single_cell/2025_06_24_mouse_and_human/2026_01_11_a375_single_cell_data.rds")
# translate capture cell id sequences to polyA 10x cell ids
translation = fread("/dartfs/rc/lab/M/McKennaLab/resources/cellranger_versions/cellranger-9.0.1/lib/python/cellranger/barcodes/translation/3M-february-2018_NXT.txt.gz",header=F)
colnames(translation) = c("feature_code","cell_id")

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

#######################################################################################################################
#
# Analyze clones 4 and 19 (now 1 and 2 on the plot) together
#
#######################################################################################################################
clone_4_tree <- read.tree("/dartfs/rc/lab/M/McKennaLab/projects/Maryam/2024_04_07_cancer_single_cell/2025_06_24_mouse_and_human/all_2026_01_29_trees/cluster_4_test_phylip.mix.treefile")
clone_4_mapping <- fread("/dartfs/rc/lab/M/McKennaLab/projects/Maryam/2024_04_07_cancer_single_cell/2025_06_24_mouse_and_human/all_2026_01_29_trees/cluster_4_test_phylip_mapping.txt")


lst <- setNames(as.list(clone_4_mapping$cell), clone_4_mapping$id)
# clone_4_tree_renamed <- makeNodeLabel(clone_4_tree, "u", nodeList = lst)

lookup <- setNames(clone_4_mapping$cell,clone_4_mapping$id)
clone_4_tree$tip.label <- lookup[clone_4_tree$tip.label]

##############################################################################
# uncomment this next section to remove 19 from the tree
##############################################################################
#clone_19_tree <- read.tree("/dartfs/rc/lab/M/McKennaLab/projects/Maryam/2024_04_07_cancer_single_cell/2025_06_24_mouse_and_human/all_2026_01_29_trees/cluster_19_test_phylip.mix.treefile")
#clone_19_mapping <- fread("/dartfs/rc/lab/M/McKennaLab/projects/Maryam/2024_04_07_cancer_single_cell/2025_06_24_mouse_and_human/all_2026_01_29_trees/cluster_19_test_phylip_mapping.txt")


#lst <- setNames(as.list(clone_19_mapping$cell), clone_19_mapping$id)
# clone_4_tree_renamed <- makeNodeLabel(clone_4_tree, "u", nodeList = lst)

#lookup <- setNames(clone_19_mapping$cell,clone_19_mapping$id)
#clone_19_tree$tip.label <- lookup[clone_19_tree$tip.label]

#clone_4_tree <- bind.tree(clone_4_tree, clone_19_tree)

Winv <- inv_tree_dist(clone_4_tree, node = TRUE, norm = FALSE)

# load the single cell data if needed - gex.combined.subset
# "/dartfs/rc/lab/M/McKennaLab/projects/Maryam/2024_04_07_cancer_single_cell/2025_06_24_mouse_and_human/2026_01_11_a375_single_cell_data.rds"

# tree is a "phylo" from ape
leaf_order <- as.data.frame(clone_4_tree$tip.label)
colnames(leaf_order) = c("cell")

leaf_order$rc = sapply(leaf_order$cell,reverse_complement)
leaf_order$r = sapply(leaf_order$cell,reverse_nocomplement)
leaf_order_joined = left_join(leaf_order, translation, by=c("rc"="feature_code"))
gex.tbl = gex.combined.subset[[]]
gex.tbl$full_cell_id = rownames(gex.tbl)
leaf_order_joined = left_join(leaf_order_joined, gex.tbl, by=c("cell_id"="V1"))

# reorder + subset cells in the Seurat object
seurat_obj <- gex.combined.subset
cell_id_list = rownames(seurat_obj[[]][is.element(gex.combined.subset[[]]$V1,leaf_order_joined$cell_id),])

# keep only cells present in both (and preserve tree order)
leaf_order <- leaf_order[is.element(leaf_order_joined$cell_id,gex.combined.subset[[]]$V1),]

seurat_hvg <- subset(seurat_obj, features = VariableFeatures(seurat_obj))

seurat_obj_subset = seurat_hvg[,leaf_order_joined$full_cell_id]

seurat_obj_subset

mat <- t(seurat_obj_subset[["RNA"]]$scale.data)

ttt = as.data.table(leaf_order_joined[,"seurat_clusters"])
dummies_model <- dummyVars(~ ., data=ttt)
trainData_mat <- predict(dummies_model, newdata = ttt)

mat2 = cbind(mat, trainData_mat)

# non_zero_counts = rowSums(GetAssayData(seurat_obj_subset, assay = "RNA", slot = "counts")) > 25

modxcor <- xcor(mat2, Winv)

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
df_sorted_matched = df_sorted %>% filter(Var1 == Var2)

# write.table(df_sorted,file="/dartfs/rc/lab/M/McKennaLab/projects/Maryam/2024_04_07_cancer_single_cell/2025_06_24_mouse_and_human/2026_01_29_clone_19_PATH.txt",sep="\t")





# remap the tip labels again -- this time to the cell data (not the capture sequence)

lookup_translation <- setNames(leaf_order_joined$cell_id,leaf_order_joined$cell)
clone_4_tree$tip.label <- lookup_translation[clone_4_tree$tip.label]


# do differential gene expression between clones
#de <- FindMarkers(
#  gex.combined.subset.lineage,
#  group.by = "filtered_clones",
#  ident.1 = "4",
#  ident.2 = "19",
#  test.use = "wilcox"   # default is Wilcoxon; often fine
#)

df_sorted_matched = df_sorted %>% filter(Var1 == Var2)



# dat2 = FetchData(seurat_obj_subset, vars = head(df_sorted_matched,n=400)$Var1)
dat2 = FetchData(seurat_obj_subset, vars = "GRCh38-MKI67")
# dat2 = dat2[,apply(dat2, 2, median) > 0]
colnames_dat2 = colnames(dat2)

dat2$ID_full = rownames(dat2)
dat2$ID = sapply(dat2$ID_full,function(x){
  tt = strsplit(x,"_")[[1]][2]
  return(strsplit(tt,"-")[[1]][1])
  })

dat2_melt = reshape2::melt(dat2,id.vars=c("ID"),measure.vars=colnames_dat2)
# dat2_melt2 = reshape2::melt(dat2,id.vars=c("ID"),measure.vars=c("clone"))

p <- ggtree(clone_4_tree, layout="circular", branch.length="none")

p2 <- p + 
  new_scale_fill() + 
  geom_fruit(
    data=dat2_melt,
    geom=geom_tile,
    mapping=aes(y=ID, x=variable, fill=value),
    offset=0.08,   # The distance between external layers, default is 0.03 times of x range of tree.
    pwidth=0.1 # width of the external layer, default is 0.2 times of x range of tree.
  )# + scale_fill_gradientn(
  #  colours = c("blue", "dodgerblue", "white", "orange", "red"),
  #  rescaler = ~ scales::rescale_mid(.x, mid = 1)
  #)

p2
df_sorted_matched$index = seq(1,nrow(df_sorted_matched))
ggsave(ggplot(df_sorted_matched[df_sorted_matched$index <= 10,]) + geom_bar(aes(x=index,y=Z),stat='identity',fill="#2179b4",col="black") + theme_classic() +  geom_text(aes(label = Var1, x=index,y=Z), angle = 90, vjust = -1, hjust = 1, color = "black"),
       file="/dartfs/rc/lab/M/McKennaLab/projects/Maryam/2024_04_07_cancer_single_cell/2025_06_24_mouse_and_human/2026_01_29_clone_19_PATH_clone_3_PATH_top_genes.pdf",width=8, height=8)


write.table(df_sorted,file="/dartfs/rc/lab/M/McKennaLab/projects/Maryam/2024_04_07_cancer_single_cell/2025_06_24_mouse_and_human/2026_01_29_clone_19_PATH.txt",sep="\t")
df_sorted_clone_4 = fread("/dartfs/rc/lab/M/McKennaLab/projects/Maryam/2024_04_07_cancer_single_cell/2025_06_24_mouse_and_human/2026_01_29_clone_4_PATH.txt")