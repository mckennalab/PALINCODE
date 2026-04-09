
library(ape)

tree <- read.tree("~/Desktop/2025_10_15_mix_input_one_hot.mix.contree")

Ntip <- length(tree$tip.label)
Nnode <- tree$Nnode

# build parent lookup table
parent <- integer(Ntip + Nnode)
for(i in 1:nrow(tree$edge)) {
  parent[tree$edge[i,2]] <- tree$edge[i,1]
}

branch_counts <- integer(Ntip)

for(tip in 1:Ntip) {
  count <- 0
  node <- tip
  while(node != (Ntip + 1)) {  # root is node Ntip+1
    node <- parent[node]
    if(node > Ntip) count <- count + 1  # internal node
  }
  branch_counts[tip] <- count
}

depth_table_293T = data.frame(tip = tree$tip.label, branch_points = branch_counts)
depth_table_293T$sample = "293T"

average_tree_depth = ggplot(depth_table_293T,aes(x=sample, y=branch_points)) + geom_violin() + theme_classic() + ylab("Leaf Depth") + xlab("")

ggsave(average_tree_depth,file="~/Desktop/code/palindrome_base_editing/figures/figure4/depth_per_leaf.pdf",width=2,height=6)