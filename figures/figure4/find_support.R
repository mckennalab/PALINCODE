library(ape)

tr <- read.tree("/Users/aaronmck/Desktop/code/2025_10_24_maryam_trees/2025_10_15_mix_input.mix.contree")
support_table = as.data.table(as.numeric(tr$node.label))         # vector of internal node labels
support_table$name = "twin_encoding"
# compute distance from root to every tip
d <- node.depth.edgelength(tr)  # returns distance for all nodes


tr <- read.tree("/Users/aaronmck/Desktop/code/2025_10_24_maryam_trees/2025_10_15_mix_input_one_hot.mix.contree")
support_table2 = as.data.table(as.numeric(tr$node.label))         # vector of internal node labels
support_table2$name = "one_hot"
full_table = rbind(support_table,support_table2)
full_table = full_table[complete.cases(full_table),]
colnames(full_table) = c("support","sample")

ggg = ggplot(full_table[full_table$sample == "one_hot"]) + geom_boxplot(aes(x=sample, y=support)) + 
  geom_quasirandom(aes(x=sample, y=support),size=2) + theme_classic()


qplot(as.numeric(tr$node.label))


avg_leaf_depth <- average_leaf_depth(tr)
