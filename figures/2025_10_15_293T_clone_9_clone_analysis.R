library(tidyverse)
library(data.table)


# orientation_CCTGTCATCTTAGCTAAGATGACAGG
process_counts_table = function(clone_name, file_to_load, col_name) {
  
  clfile = fread(file_to_load)
  
  clfile$sample = clone_name
  event_counts = clfile %>%
    count(e0, !!as.name(col_name)) %>%
    pivot_wider(names_from = col_name, values_from = n, values_fill = 0)
  
  event_counts$sum = rowSums(event_counts[,c(2,3,4,5)])
  event_counts = event_counts[event_counts$sum > 100,]
  event_counts$BOTH = event_counts$BOTH / event_counts$sum
  event_counts$LEFT = event_counts$LEFT / event_counts$sum
  event_counts$RIGHT = event_counts$RIGHT / event_counts$sum
  event_counts$NONE = event_counts$NONE / event_counts$sum
  event_counts$entropy = (-1.0 * event_counts$BOTH * log2(event_counts$BOTH + 0.0001)) +
    (-1.0 * event_counts$RIGHT * log2(event_counts$RIGHT + 0.0001)) + 
    (-1.0 * event_counts$LEFT * log2(event_counts$LEFT + 0.0001)) + 
    (-1.0 * event_counts$NONE * log2(event_counts$NONE + 0.0001))
  
  event_counts$target = col_name
  event_counts$sample = clone_name
  event_counts = event_counts[,c("e0","BOTH","LEFT","RIGHT","NONE","sum","entropy","target","sample")]
  #if (is.na(all_event_counts_total)) {
    #all_event_counts_total = event_counts
  #} else {
    #all_event_counts_total = rbind(all_event_counts_total,event_counts)
  #}
  
  #event_counts_melted = reshape2::melt(event_counts[,c(1,2,3,4,5)],id.var="e0")
  #event_counts_melted$sample = clone_name
  return(event_counts)
}


clone_9211_event_counts_melted_CCAGGAAGTACTCGAGTACTTCCTGG = process_counts_table("clone_9211","/dartfs/rc/lab/M/McKennaLab/projects/Maryam/2025_08_28_single_cell_clones/9_2_1_1_S1_R2_001_collapsed.txt","orientation_CCAGGAAGTACTCGAGTACTTCCTGG")
clone_9212_event_counts_melted_CCAGGAAGTACTCGAGTACTTCCTGG = process_counts_table("clone_9212","/dartfs/rc/lab/M/McKennaLab/projects/Maryam/2025_08_28_single_cell_clones/9_2_1_2_S2_R2_001_collapsed.txt")
clone_9213_event_counts_melted_CCAGGAAGTACTCGAGTACTTCCTGG = process_counts_table("clone_9213","/dartfs/rc/lab/M/McKennaLab/projects/Maryam/2025_08_28_single_cell_clones/9_2_1_3_S3_R2_001_collapsed.txt")
clone_9214_event_counts_melted_CCAGGAAGTACTCGAGTACTTCCTGG = process_counts_table("clone_9214","/dartfs/rc/lab/M/McKennaLab/projects/Maryam/2025_08_28_single_cell_clones/9_2_1_4_S4_R2_001_collapsed.txt")
clone_9215_event_counts_melted_CCAGGAAGTACTCGAGTACTTCCTGG = process_counts_table("clone_9215","/dartfs/rc/lab/M/McKennaLab/projects/Maryam/2025_08_28_single_cell_clones/9_2_1_5_S5_R2_001_collapsed.txt")
clone_9216_event_counts_melted_CCAGGAAGTACTCGAGTACTTCCTGG = process_counts_table("clone_9216","/dartfs/rc/lab/M/McKennaLab/projects/Maryam/2025_08_28_single_cell_clones/9_2_1_6_S6_R2_001_collapsed.txt")
clone_921b1_event_counts_melted_CCAGGAAGTACTCGAGTACTTCCTGG = process_counts_table("clone_921b1","/dartfs/rc/lab/M/McKennaLab/projects/Maryam/2025_08_28_single_cell_clones/9_2_72v1_S7_R2_001_collapsed.txt")
clone_921b3_event_counts_melted_CCAGGAAGTACTCGAGTACTTCCTGG = process_counts_table("clone_921b3","/dartfs/rc/lab/M/McKennaLab/projects/Maryam/2025_08_28_single_cell_clones/9_2_72v3_S8_R2_001_collapsed.txt")


all_calls_melted_CCAGGAAGTACTCGAGTACTTCCTGG = rbindlist(list(clone_9211_event_counts_melted_CCAGGAAGTACTCGAGTACTTCCTGG,
                                  clone_9212_event_counts_melted_CCAGGAAGTACTCGAGTACTTCCTGG,
                                  clone_9213_event_counts_melted_CCAGGAAGTACTCGAGTACTTCCTGG,
                                  clone_9214_event_counts_melted_CCAGGAAGTACTCGAGTACTTCCTGG,
                                  clone_9215_event_counts_melted_CCAGGAAGTACTCGAGTACTTCCTGG,
                                  clone_9216_event_counts_melted_CCAGGAAGTACTCGAGTACTTCCTGG,
                                  clone_921b1_event_counts_melted_CCAGGAAGTACTCGAGTACTTCCTGG,
                                  clone_921b3_event_counts_melted_CCAGGAAGTACTCGAGTACTTCCTGG))

clone_9211_event_counts__CCTGTCATCTTAGCTAAGATGACAGG = process_counts_table("1.1","/dartfs/rc/lab/M/McKennaLab/projects/Maryam/2025_08_28_single_cell_clones/9_2_1_1_S1_R2_001_collapsed.txt","orientation_CCTGTCATCTTAGCTAAGATGACAGG")
clone_9212_event_counts__CCTGTCATCTTAGCTAAGATGACAGG = process_counts_table("1.2","/dartfs/rc/lab/M/McKennaLab/projects/Maryam/2025_08_28_single_cell_clones/9_2_1_2_S2_R2_001_collapsed.txt","orientation_CCTGTCATCTTAGCTAAGATGACAGG")
clone_9213_event_counts__CCTGTCATCTTAGCTAAGATGACAGG = process_counts_table("1.3","/dartfs/rc/lab/M/McKennaLab/projects/Maryam/2025_08_28_single_cell_clones/9_2_1_3_S3_R2_001_collapsed.txt","orientation_CCTGTCATCTTAGCTAAGATGACAGG")
clone_9214_event_counts__CCTGTCATCTTAGCTAAGATGACAGG = process_counts_table("1.4","/dartfs/rc/lab/M/McKennaLab/projects/Maryam/2025_08_28_single_cell_clones/9_2_1_4_S4_R2_001_collapsed.txt","orientation_CCTGTCATCTTAGCTAAGATGACAGG")
clone_9215_event_counts__CCTGTCATCTTAGCTAAGATGACAGG = process_counts_table("1.5","/dartfs/rc/lab/M/McKennaLab/projects/Maryam/2025_08_28_single_cell_clones/9_2_1_5_S5_R2_001_collapsed.txt","orientation_CCTGTCATCTTAGCTAAGATGACAGG")
clone_9216_event_counts__CCTGTCATCTTAGCTAAGATGACAGG = process_counts_table("1.6","/dartfs/rc/lab/M/McKennaLab/projects/Maryam/2025_08_28_single_cell_clones/9_2_1_6_S6_R2_001_collapsed.txt","orientation_CCTGTCATCTTAGCTAAGATGACAGG")
clone_921b1_event_counts__CCTGTCATCTTAGCTAAGATGACAGG = process_counts_table("1B","/dartfs/rc/lab/M/McKennaLab/projects/Maryam/2025_08_28_single_cell_clones/9_2_72v1_S7_R2_001_collapsed.txt","orientation_CCTGTCATCTTAGCTAAGATGACAGG")
clone_921b3_event_counts__CCTGTCATCTTAGCTAAGATGACAGG = process_counts_table("3B","/dartfs/rc/lab/M/McKennaLab/projects/Maryam/2025_08_28_single_cell_clones/9_2_72v3_S8_R2_001_collapsed.txt","orientation_CCTGTCATCTTAGCTAAGATGACAGG")



all_calls__CCTGTCATCTTAGCTAAGATGACAGG = rbindlist(list(clone_9211_event_counts__CCTGTCATCTTAGCTAAGATGACAGG,
                                  clone_9212_event_counts__CCTGTCATCTTAGCTAAGATGACAGG,
                                  clone_9213_event_counts__CCTGTCATCTTAGCTAAGATGACAGG,
                                  clone_9214_event_counts__CCTGTCATCTTAGCTAAGATGACAGG,
                                  clone_9215_event_counts__CCTGTCATCTTAGCTAAGATGACAGG,
                                  clone_9216_event_counts__CCTGTCATCTTAGCTAAGATGACAGG,
                                  clone_921b1_event_counts__CCTGTCATCTTAGCTAAGATGACAGG,
                                  clone_921b3_event_counts__CCTGTCATCTTAGCTAAGATGACAGG))

#melted_cor_CCTGTCATCTTAGCTAAGATGACAGG = reshape2::melt(cor_CCTGTCATCTTAGCTAAGATGACAGG)

# all_calls_melted_CCTGTCATCTTAGCTAAGATGACAGG_left = all_calls_melted_CCTGTCATCTTAGCTAAGATGACAGG[all_calls_melted_CCTGTCATCTTAGCTAAGATGACAGG$variable == "LEFT",]
none_values_CCTGTCATCTTAGCTAAGATGACAGG = dcast(all_calls__CCTGTCATCTTAGCTAAGATGACAGG, formula = e0 ~ sample, value.var = "BOTH")
none_values_CCTGTCATCTTAGCTAAGATGACAGG[is.na(none_values_CCTGTCATCTTAGCTAAGATGACAGG)] = 0
cor_CCTGTCATCTTAGCTAAGATGACAGG = cor(none_values_CCTGTCATCTTAGCTAAGATGACAGG[,c(2,3,4,5,6,7,8,9)])
melted_cor_CCTGTCATCTTAGCTAAGATGACAGG = reshape2::melt(cor_CCTGTCATCTTAGCTAAGATGACAGG)



ggplot(melted_cor_CCTGTCATCTTAGCTAAGATGACAGG) + 
  geom_tile(aes(x=Var1,y=Var2,fill=value)) +
  scale_fill_gradientn(colours = c("blue", "white", "red")) #+ 
  #geom_text(aes(x=Var1,y=Var2,label=round(value, 2)))



# put them together into a master table
all_calls = rbindlist(list(clone_9211,clone_9212,clone_9213,clone_9214,clone_9215,clone_9216,clone_921b1,clone_921b3))



e0_table = as.data.table(table(all_calls$e0,all_calls$sample))

all_common = e0_table[e0_table$N > 7000,]
uncommon = table(all_common$V1)[table(all_common$V1) < 8]


# Install if needed
if (!requireNamespace("stringdist", quietly = TRUE)) {
  install.packages("stringdist")
}

library(stringdist)

# All-to-all Levenshtein distance matrix
pairwise_levenstein <- function(sequences) {
  n <- length(sequences)
  dist_matrix <- matrix(0, nrow = n, ncol = n)
  rownames(dist_matrix) <- colnames(dist_matrix) <- sequences
  
  for (i in 1:n) {
    for (j in i:n) {
      d <- stringdist::stringdist(sequences[i], sequences[j], method = "lv")
      dist_matrix[i, j] <- d
      dist_matrix[j, i] <- d  # symmetric
    }
  }
  return(dist_matrix)
}

all_common_distances = pairwise_levenstein(unique(all_common$V1))
all_common_distances_melt = reshape2::melt(all_common_distances)
all_common_tagged_true = all_common[all_common$tagged,]

too_close_barcodes = all_common_distances_melt[as.character(all_common_distances_melt$Var1) > as.character(all_common_distances_melt$Var2) & 
                                                 all_common_distances_melt$Var1 != all_common_distances_melt$Var2 & 
                                                 all_common_distances_melt$value < 4,]

too_close_barcodes_list = c(too_close_barcodes$Var1,too_close_barcodes$Var2)

all_common$tagged = FALSE
all_common$tagged[is.element(all_common$V1,too_close_barcodes_list)] = TRUE


ggplot(df, aes(x = Sample, y = Gene, fill = Expression)) +
  geom_tile() +
  scale_fill_viridis_c() +  # Or use scale_fill_gradient(), scale_fill_gradient2(), etc.
  theme_minimal() +
  labs(title = "Heatmap", x = "Sample", y = "Gene")

counts_AAGGCGATCACA = as.data.table(table(all_calls[all_calls$e0 == "AAGGCGATCACA",]$CCAGGAAGTACTCGAGTACTTCCTGG,all_calls[all_calls$e0 == "AAGGCGATCACA",]$sample))
counts_TCGCCAAGTGTT = as.data.table(table(all_calls[all_calls$e0 == "TCGCCAAGTGTT",]$CCAGGAAGTACTCGAGTACTTCCTGG,all_calls[all_calls$e0 == "TCGCCAAGTGTT",]$sample))

props = all_calls[is.element(all_calls$e0,all_common[!all_common$tagged,]$V1),] %>%
  group_by(sample, e0, CCAGGAAGTACTCGAGTACTTCCTGG) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(sample, e0) %>%
  mutate(prop = n / sum(n))

top_events = props[props$prop > 0.05,]


props_RNF2 = all_calls[is.element(all_calls$e0,all_common[!all_common$tagged,]$V1),] %>%
  group_by(sample, e0, CCTGTCATCTTAGCTAAGATGACAGG) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(sample, e0) %>%
  mutate(prop = n / sum(n))

top_events_RNF2  = props_RNF2[props_RNF2$prop > 0.05,]

ggplot(top_events_RNF2) + geom_tile(aes(x=sample,y=CCTGTCATCTTAGCTAAGATGACAGG,fill=prop)) + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + facet_wrap(. ~ top_events$e0)