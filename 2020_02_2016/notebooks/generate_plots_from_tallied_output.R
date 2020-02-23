library(ggplot2)
library(data.table)
library(UpSetR)
library(dplyr)

editing_tbl = fread("~/Desktop/code/palindrome_base_editing/2020_02_2016/data/../notebooks/tallied_outcomes.tsv")
known_editing_tbl = editing_tbl[editing_tbl$isKnown,]

occurances = as.data.frame.matrix(table(known_editing_tbl$guide,known_editing_tbl$lib))

png(file="~/Desktop/code/palindrome_base_editing/2020_02_2016/plots/2020_02_23_library_occurance_plot.png",width = 1200,height=600) # or other device
upset(occurances,sets.bar.color = "#56B4E9",
      order.by = "freq", empty.intersections = "on", text.scale=3, point.size=5)
dev.off()

occurances$sum = rowSums(occurances)
all_captures = row.names(occurances)[occurances$sum == 3]
all_captures_editing_table = known_editing_tbl[is.element(known_editing_tbl$guide,all_captures),]
melted_all_captures_editing_table = melt(all_captures_editing_table,
                                         id.vars = c("guide","lib"), 
                                         measure.vars = c("none","left","right","both","bad"))


wt = ggplot(melted_all_captures_editing_table[melted_all_captures_editing_table$lib == "WT",], aes(variable,guide, fill =value)) + 
  geom_tile() + 
  scale_fill_gradient(low = "white", high = "black") + 
  theme(axis.text.y = element_text(size = 5)) + 
  xlab("Editing outcome") +
  ylab("Targets sequence")
ggsave(wt,file="~/Desktop/code/palindrome_base_editing/2020_02_2016/plots/WT_base_editing.png",width=6,height=14)

abe = ggplot(melted_all_captures_editing_table[melted_all_captures_editing_table$lib == "ABE",], aes(variable,guide, fill =value)) + 
  geom_tile() + 
  scale_fill_gradient(low = "white", high = "black") + 
  theme(axis.text.y = element_text(size = 5))+ 
  xlab("Editing outcome") +
  ylab("Targets sequence")
ggsave(abe,file="~/Desktop/code/palindrome_base_editing/2020_02_2016/plots/ABE_base_editing.png",width=6,height=14)

cbe = ggplot(melted_all_captures_editing_table[melted_all_captures_editing_table$lib == "CBE",], aes(variable,guide, fill =value)) + 
  geom_tile() + 
  scale_fill_gradient(low = "white", high = "black") + 
  theme(axis.text.y = element_text(size = 5))+ 
  xlab("Editing outcome") +
  ylab("Targets sequence")
ggsave(cbe,file="~/Desktop/code/palindrome_base_editing/2020_02_2016/plots/CBE_base_editing.png",width=6,height=14)