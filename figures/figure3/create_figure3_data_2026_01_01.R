library(grid)
library(gtable)

table_all = fread("/Users/aaronmck/Desktop/code/palindrome_base_editing/2025_11_30_twist_analysis/Data/full_twist2_data.csv")

filtered_table = table_all[table_all$closestMatch,]

filtered_table_counts = t(table(filtered_table$editing,filtered_table$sample))

filtered_table_prop = prop.table(filtered_table_counts, margin = 1)

filtered_table_prop_melt = reshape2::melt(filtered_table_prop)

# #############################################
# Load up the original data
# #############################################
table_original = fread("/Users/aaronmck/Desktop/code/palindrome_base_editing/2025_11_30_twist_analysis/Data/all_dat.csv")

filtered_table_original = table_original[table_original$found_known & table_original$guide_target_match,]

# full outcome table - fill in from      8e,  abemax, or control 
full_table_filtered_table_original_all = list()
for (tag in c("abemax","control","8e" )) {
  filtered_table_original_control = filtered_table_original[filtered_table_original$BE == tag]
  filtered_table_original_per_target = as.data.table(prop.table(table(filtered_table_original_control$known_target, filtered_table_original_control$editing), margin = 1))
  filtered_table_original_per_target_editing = filtered_table_original_per_target[filtered_table_original_per_target$V2 == "NEITHER",]
  filtered_table_original_per_target_editing$editing = 1.0 - filtered_table_original_per_target_editing$N
  filtered_table_original_per_target_editing$V1 <- reorder(filtered_table_original_per_target_editing$V1, filtered_table_original_per_target_editing$editing)
  filtered_table_original_per_target_editing$tag = tag
  print(tag)
  g3 = ggplot(filtered_table_original_per_target_editing) + geom_point(aes(V1,editing),size=5,fill="white",color="black", pch=21) + theme_classic()
  g <- ggplotGrob(g3)
  grid::gpar()  # ensure grobs registered
  ggsave(g,file=paste("~/Desktop/code/palindrome_base_editing/figures/figure3/",tag,"_point_plot.pdf",sep=""),width=8,height=4,dpi=300)
  
  full_table_filtered_table_original_all = c(full_table_filtered_table_original_all,list(filtered_table_original_per_target_editing))
}

# just look into the 8e data set
tag = "8e"
filtered_table_original_control = filtered_table_original[filtered_table_original$BE == tag]
filtered_table_original_per_target = as.data.table(prop.table(table(filtered_table_original_control$known_target, filtered_table_original_control$editing), margin = 1))
filtered_table_original_per_target_editing = filtered_table_original_per_target[filtered_table_original_per_target$V2 == "NEITHER",]
filtered_table_original_per_target_editing$editing = 1.0 - filtered_table_original_per_target_editing$N
filtered_table_original_per_target_editing$V1 <- reorder(filtered_table_original_per_target_editing$V1, filtered_table_original_per_target_editing$editing)
filtered_table_original_per_target_editing$tag = tag

# candidate 1:
# CCTAGATCGCGCGATCTAGGGGG, or GTAGACCGGCGCCGGTCTACAGG
table(filtered_table_original[filtered_table_original$known_target == "CCTAGATCGCGCGATCTAGGGGG" & filtered_table_original$BE == "8e",]$editing)

# candidate 2:
# CCGCGTAGCGCGCTACGCGGAGG

eee <- rbindlist(full_table_filtered_table_original_all)
wide <- eee[,c("V1","tag","editing")] %>%
  pivot_wider(names_from = tag, values_from = editing)
wide[is.na(wide)] = 0

fit <- lm(abemax ~ `8e`, wide)
slope <- round(coef(fit)[2], 3)
intercept <- round(coef(fit)[1], 3)
r <- round(cor(wide$abemax, wide$`8e`), 3)

ggg = ggplot(wide, aes(abemax, `8e`)) +
  geom_point() +
  geom_smooth(method="lm", se=FALSE, color="blue") +
  annotate("text",
           x = Inf, y = -Inf,        # bottom-right corner
           hjust = 1.1, vjust = -0.5,
           size = 5,
           label = paste0(
             "Slope = ", slope, "\n",
             "Intercept = ", intercept, "\n",
             "r = ", r
           )) + theme_classic()
ggsave(ggg,file="~/Desktop/code/palindrome_base_editing/figures/sup_figure4/abemax_8e_corr.pdf",width=8,height=7)



filtered_table_original_counts_prop = as.data.table(prop.table(filtered_table_original_counts, margin = 1))

break_up_editor = function(x) {
  splitted = strsplit(x,"_")[1][[1]]
  return(list("editor" = splitted[1], "scaffold" = splitted[2], "length" = splitted[3]))
}

filtered_table_original_counts_prop[, c("editor", "scaffold", "length") := break_up_editor(V1), by = 1:nrow(filtered_table_original_counts_prop)]

editing_rates_plot = ggplot(filtered_table_original_counts_prop[filtered_table_original_counts_prop$V2 != "NEITHER" & filtered_table_original_counts_prop$scaffold != "new"]) + 
  geom_bar(aes(x=V2, y=N, fill=V2), col="black",stat="identity") + 
  facet_wrap(editor ~ length) + 
  theme_classic()


editing_rates_plot = ggplot(filtered_table_original_counts_prop[filtered_table_original_counts_prop$V2 != "NEITHER"]) + 
  geom_bar(aes(x=V2, y=N, fill=V2),stat="identity") + 
  facet_wrap(editor ~ length) + 
  theme_classic()

