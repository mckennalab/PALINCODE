neither = fread("~/Desktop/code/python_clique/clique/test_output.txt_unedited")
right = fread("~/Desktop/code/python_clique/clique/test_output.txt_right")
left = fread("~/Desktop/code/python_clique/clique/test_output.txt_left")

left_tab = as.data.frame(prop.table(table(left[left$alignment_rate > .95,"CCAGGAAGTGCTCGAGTACTTCCTGG_nlrb"])))
left_tab$exp = "left"
right_tab = as.data.frame(prop.table(table(right[right$alignment_rate > .95,"CCAGGAAGTACTCGAGCACTTCCTGG_nlrb"])))
right_tab$exp = "right"
neither_tab = as.data.frame(prop.table(table(neither[neither$alignment_rate > .95,"CCAGGAAGTACTCGAGTACTTCCTGG_nlrb"])))
neither_tab$exp = "wild_type"

all_table = rbind(left_tab,right_tab,neither_tab)
write.table(all_table,file="~/Desktop/code/palindrome_base_editing/figures/figure2/all_table.txt",quote=F, sep="\t", row.names=F)

g = ggplot(all_table[all_table$Var1 != "NONE" & all_table$Var1 != "BOTH",]) + 
  geom_bar(aes(y=Freq,x=Var1,fill=exp),stat='identity',position='dodge') + 
  facet_wrap(. ~ exp) + 
  theme_classic() 

ggsave(g,file="~/Desktop/code/palindrome_base_editing/figures/figure2/lnrb_test.pdf",width=5,height=9)