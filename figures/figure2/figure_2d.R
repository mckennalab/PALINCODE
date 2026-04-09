library(data.table)
library(tidyverse)
library(tidyinftheo)
library(ggpattern)


# read in the master table of editing
cc = fread("~/Desktop/concat.txt.gz")

guide_length_colors = c("#001845", "#00a6fb", "#F6AE2D", "#F26419", "#ad2831")

# create a master name for the samples, combining a number of the parmaters
cc$full_name = paste(cc$sample,cc$experiment,cc$editing_type,cc$tracr_type,sep="_")


# load our samples for plot in figure 2a -- T7 from the genome and PalT7 introduced into the genome
figure_2a_data = cc[
                      (cc$experiment == "MF_Lib_05" & cc$sample == "T7_Control") | 
                      (cc$experiment == "MF_Lib_11" & cc$sample == "NPT7_g20_Omx") |
                      (cc$experiment == "MF_Lib_08" & cc$sample == "RNF2_control") |
                      (cc$experiment == "MF_Lib_08" & cc$sample == "RNF2_max") |
                      (cc$experiment == "MF_Lib_08" & cc$sample == "PalT7-control") |
                      (cc$experiment == "MF_lib_07" & cc$sample == "PalT7-g20" & cc$editing_type == "ABEMAX") |
                      (cc$experiment == "MF_Lib_08" & cc$sample == "Pal_RNF2-control") |
                      (cc$experiment == "MF_Lib_08" & cc$sample == "Pal_RNF2-g20" & cc$editing_type == "ABEMAX"),]
                      
                      


# summerize the wild-type proportion -- some point changes or small indels make it so that we have a set 
# of < 1% unedited that we should account for. 

# first make a map of what editing we'll consider for T7 and PalT7 -- positions that are edited from A->G or from T->C
# PalT7 - CCAGGAAGTACTCGAGTACTTCCTGG
#            ..**..*.*..*.**.**..*..
# T7 -       GGAAGTACTCGCACACGGGATGG
#            ..**.**.*...*.*....**..
t7_columns = c("2_G","3_G","5_C","6_G","8_C","12_G","14_G","19_G")
pal_t7_columns = c("2_G","3_G","5_C","6_G","8_C","11_G","13_C","14_G","16_C","17_C")

# same for RNF2/PalRNF2
#              .*.**.***..***.**.*.
# RNF2_pal	CCTGTCATCTTAGCTAAGATGACAGG
# RNF2	       GTCATCTTAGTCATTACCTGAGG
#              .*.**.***.*.****..*.
pal_rnf2_columns = c("1_C","3_G","4_C","6_C","7_C","8_G","11_C","12_G","13_G","15_G","16_C","18_G")
rnf2_columns = c("1_C","3_G","4_C","6_C","7_C","8_G","10_C","12_G","13_C","14_C","15_G","18_C")


# now figure out what alleles are truely edited and use that to subtract from the total for our editing proportions
figure_2a_data_filtered_editing = figure_2a_data[(rowSums(figure_2a_data[,..t7_columns]) >= 1 & figure_2a_data$palindrome == F) |
                                                 (rowSums(figure_2a_data[,..pal_t7_columns]) >= 1 & figure_2a_data$palindrome == T) |
                                                 (rowSums(figure_2a_data[,..rnf2_columns]) >= 1 & figure_2a_data$palindrome == F) |
                                                 (rowSums(figure_2a_data[,..pal_rnf2_columns]) >= 1 & figure_2a_data$palindrome == T),]

# add up all the editing proportions for each sample
figure_2a_data_collapsed = aggregate(figure_2a_data_filtered_editing$proportion, by=list(sample=figure_2a_data_filtered_editing$sample,
                                                                                         etype=figure_2a_data_filtered_editing$editing_type,
                                                                                         target=figure_2a_data_filtered_editing$target,
                                                                                         palindrome=figure_2a_data_filtered_editing$palindrome), FUN=sum)
# order the samples in a reasonable way
figure_2a_data_collapsed$sample = factor(figure_2a_data_collapsed$sample,
                                         levels=c("T7_Control","NPT7_g20_Omx",
                                                  "PalT7-control","PalT7-g20",
                                                  "RNF2_control","RNF2_max",
                                                  "Pal_RNF2-control","Pal_RNF2-g20"))

# now do a plot with all the samples together, coloring and shading by target / type
ggplot(data = figure_2a_data_collapsed, aes(y = x, x=sample, fill = target, pattern = palindrome)) +
  geom_bar_pattern(position = position_dodge(preserve = "single"),
                   color = "black", 
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_density = 0.1,
                   pattern_spacing = 0.025,
                   pattern_key_scale_factor = 0.6, stat='identity') + 
  scale_fill_manual(values = colorRampPalette(c("#0066CC","#0066CC","#FF8C00","#FF8C00"))(4)) +
  scale_pattern_manual(values = c("TRUE" = "stripe", "FALSE" = "none")) +
  labs(x = "Target", y = "Overall editing rate", pattern = "palindrome sequence") + 
  guides(pattern = guide_legend(override.aes = list(fill = "white")),
         fill = guide_legend(override.aes = list(pattern = "none"))) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

