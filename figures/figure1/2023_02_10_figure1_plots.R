library(data.table)
library(tidyverse)
library(tidyinftheo)
library(ggpattern)

#guide_length_colors = c("#001845", "#00a6fb", "#F6AE2D", "#F26419", "#ad2831")
#guide_length_colors = c("#c7f9cc", "#80ed99", "#57cc99", "#38a3a5", "#22577a")
#guide_length_colors = c("#c7f9cc", "#80ed99", "#57cc99", "#38a3a5", "#22577a")
guide_length_colors = c("#eaac8b", "#e56b6f", "#b56576", "#6d597a", "#355070")


# read in the master table of editing
cc = fread("~/Desktop/code/palindrome_base_editing/figures/figure1/2023_02_10_full_early_editing.txt.gz")


cc$full_name = paste(cc$sample,cc$experiment,cc$editing_type,cc$tracr_type,sep="_")


# just look at experiment 7, ABEMAX
editing_summary_MF_lib_06 = cc[cc$experiment == "MF_Lib_06"] %>% 
  group_by(full_name,left_right_both,tracr_type,sample,experiment,editing_type,guide_length) %>% 
  summarise(editing_sum=sum(proportion),.groups = 'drop')

editing_summary_MF_lib_06$left_right_both = factor(editing_summary_MF_lib_06$left_right_both, levels = c("LEFT", "BOTH", "RIGHT", "NEITHER"))

# editing by L/R/B/N
g = ggplot(editing_summary_MF_lib_06[editing_summary_MF_lib_06$left_right_both != "NEITHER",]) + 
  geom_bar(aes(x=left_right_both,y=editing_sum,fill=as.factor(guide_length)),col="black",stat='identity',position='dodge') + 
  scale_fill_manual(values=guide_length_colors) + 
  theme_classic() + facet_wrap(. ~ sample)

# just look at experiment 7, ABEMAX
editing_summary_MF_lib_07 = cc[cc$experiment == "MF_lib_07" & cc$editing_type == "ABEMAX"] %>% 
  group_by(full_name,left_right_both,tracr_type,sample,experiment,editing_type,guide_length) %>% 
  summarise(editing_sum=sum(proportion),.groups = 'drop')

editing_summary$left_right_both = factor(editing_summary$left_right_both, levels = c("LEFT", "BOTH", "RIGHT", "NEITHER"))

# editing by L/R/B/N
g = ggplot(editing_summary[editing_summary$left_right_both != "NEITHER",]) + 
  geom_bar(aes(x=left_right_both,y=editing_sum,fill=as.factor(guide_length)),col="black",stat='identity',position='dodge') + 
  scale_fill_manual(values=guide_length_colors) + 
  theme_classic()+ facet_wrap(. ~ sample,nrow = 1)


# plot the left, right, or both editing for the PalT7 
ggplot(editing_summary[editing_summary$left_right_both != "NEITHER",],
       aes(y=as.factor(guide_length),x=left_right_both,fill=sum_salary)) + 
  geom_tile() + 
  scale_fill_gradient(low = "white", high = "red") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_text(aes(label = round(sum_salary,2)), color = "black", size = 4) +
  theme_classic()

# create a function to capture the entropy of the editing outcomes -- something we're trying to maximize
entropy <- function(p) sum(-(p * log2(p)), na.rm = TRUE)

# now summerise the editing entropy and plot it
entropy_values = editing_summary %>% group_by(full_name) %>% summarise(entropy_value=entropy(sum_salary))
editing_summary = left_join(editing_summary,entropy_values)

g = ggplot(editing_summary[editing_summary$left_right_both == "NEITHER",]) + 
  geom_bar(aes(x=as.factor(guide_length),y=entropy_value,fill=as.factor(guide_length)),
           col="black",stat='identity') + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  theme_classic() + 
  scale_fill_manual(values=guide_length_colors) +
  ylab("Entropy value") +
  xlab("CRISPR guide length")

g = ggplot(editing_summary[editing_summary$left_right_both == "NEITHER",]) + 
  geom_bar(aes(x=as.factor(guide_length),y=1.0 - sum_salary,fill=as.factor(guide_length)),
           col="black",stat='identity') + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  theme_classic() + 
  scale_fill_manual(values=guide_length_colors) +
  ylab("Base editing proportion") +
  xlab("CRISPR guide length")


# now work in the melting temperature information to the plots
melting_temp = data.table(
  guide_length = c(17,18,20,24,27),
  tm = c(-3.6,-4.5,-9.6,-9.8,-9.8))

# combine the two
editing_summary_with_melting_temps = left_join(editing_summary,melting_temp)

# and plot their correlation
just_edit_data = editing_summary_with_melting_temps[editing_summary_with_melting_temps$left_right_both == "NEITHER",]
ggplot(just_edit_data,aes(y=1- sum_salary,x=tm)) + 
  geom_line() + 
  geom_point(aes(col=as.factor(guide_length))) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  theme_classic() + 
  scale_color_manual(values=guide_length_colors)

