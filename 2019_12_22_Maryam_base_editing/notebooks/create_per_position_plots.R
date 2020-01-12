library("ggplot2")
library("data.table")

tbl <- fread("/Users/aaronmck/Desktop/code/base_editing/2019_12_22_Maryam_base_editing/data/all_samples.txt")
tbl$total = tbl$A + tbl$C + tbl$G + tbl$T
tbl$A.prop = tbl$A / tbl$total
tbl$C.prop = tbl$C / tbl$total
tbl$G.prop = tbl$G / tbl$total
tbl$T.prop = tbl$T / tbl$total
tbl$AT.prop = tbl$A.prop + tbl$T.prop
tbl$CG.prop = tbl$C.prop + tbl$G.prop

pal1.first20 = tbl[tbl$pal_tag == "pal1",][1:20,]
pal3.first20 = tbl[tbl$pal_tag == "pal3",][1:20,]

gg = ggplot(tbl[tbl$index < 20 & tbl$pal_tag == "pal1",]) + 
  geom_bar(aes(x=index,y=1.0 - control_wt,fill=control),stat='identity') + 
  facet_wrap(. ~ sample,ncol = 3) + 
  theme_classic() + 
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  scale_x_continuous(labels = pal1.first20$base,breaks=pal1.first20$index)
ggsave(gg,
       file="/Users/aaronmck/Desktop/code/base_editing/2019_12_22_Maryam_base_editing/plots/editing_per_base_pal1.png",
       width=12,
       height=8)
gg = ggplot(tbl[tbl$index < 20 & tbl$pal_tag == "pal3",]) + 
  geom_bar(aes(x=index,y=1.0 - control_wt,fill=control),stat='identity') + 
  facet_wrap(. ~ sample,ncol = 3) + 
  theme_classic() + 
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  scale_x_continuous(labels = pal3.first20$base,breaks=pal3.first20$index)
ggsave(gg,
       file="/Users/aaronmck/Desktop/code/base_editing/2019_12_22_Maryam_base_editing/plots/editing_per_base_pal3.png",
       width=12,
       height=8)

gg = ggplot(tbl[tbl$index < 20  & tbl$pal_tag == "pal1",]) + 
  geom_bar(aes(x=index,y=control_wt - wt,fill=control),stat='identity') + 
  facet_wrap(. ~ sample,ncol = 3) + 
  theme_classic() + 
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  scale_x_continuous(labels = pal1.first20$base,breaks=pal1.first20$index)
ggsave(gg,
       file="/Users/aaronmck/Desktop/code/base_editing/2019_12_22_Maryam_base_editing/plots/editing_control_subtracted_pal1.png",
       width=12,
       height=8)
gg = ggplot(tbl[tbl$index < 20 & tbl$pal_tag == "pal3",]) + 
  geom_bar(aes(x=index,y=control_wt - wt,fill=control),stat='identity') + 
  facet_wrap(. ~ sample,ncol = 3) + 
  theme_classic() + 
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  scale_x_continuous(labels = pal3.first20$base,breaks=pal3.first20$index)
ggsave(gg,
       file="/Users/aaronmck/Desktop/code/base_editing/2019_12_22_Maryam_base_editing/plots/editing_control_subtracted_pal3.png",
       width=12,
       height=8)

melt.tbl = melt(tbl[tbl$index < 20,],measure.vars = c("A.prop","C.prop","G.prop","T.prop"),id.vars = c("sample","index", "pal_tag"))
melt.CG.tbl = melt(tbl[tbl$index < 20,],measure.vars = c("AT.prop","CG.prop"),id.vars = c("sample","index", "pal_tag"))

gg = ggplot(melt.tbl[melt.tbl$pal_tag == "pal1",]) + 
  geom_bar(aes(x=index,y=value,fill=variable),stat='identity') + 
  facet_wrap(. ~ sample,ncol = 3) + 
  theme_classic() + 
  scale_fill_brewer(palette="Set1") +
  scale_x_continuous(labels = pal1.first20$base,breaks=pal1.first20$index)
ggsave(gg,
       file="/Users/aaronmck/Desktop/code/base_editing/2019_12_22_Maryam_base_editing/plots/editing_by_base_pal1.png",
       width=12,
       height=8)
gg = ggplot(melt.tbl[melt.tbl$pal_tag == "pal3",]) + 
  geom_bar(aes(x=index,y=value,fill=variable),stat='identity') + 
  facet_wrap(. ~ sample,ncol = 3) + 
  theme_classic() + 
  scale_fill_brewer(palette="Set1") +
  scale_x_continuous(labels =pal3.first20$base,breaks=pal3.first20$index)
ggsave(gg,
       file="/Users/aaronmck/Desktop/code/base_editing/2019_12_22_Maryam_base_editing/plots/editing_by_base_pal3.png",
       width=12,
       height=8)



gg = ggplot(melt.CG.tbl[melt.CG.tbl$pal_tag == "pal1",]) + 
  geom_bar(aes(x=index,y=value,fill=variable),stat='identity') + 
  facet_wrap(. ~ sample,ncol = 3) + 
  theme_classic() + 
  scale_fill_manual(values=c("#999999","#56B4E9")) + 
  scale_x_continuous(labels = pal1.first20$base,breaks=pal1.first20$index)
ggsave(gg,
       file="/Users/aaronmck/Desktop/code/base_editing/2019_12_22_Maryam_base_editing/plots/editing_AT_GC_split_pal1.png",
       width=12,
       height=8)

gg = ggplot(melt.CG.tbl[melt.CG.tbl$pal_tag == "pal3",]) + 
  geom_bar(aes(x=index,y=value,fill=variable),stat='identity') + 
  facet_wrap(. ~ sample,ncol = 3) + 
  theme_classic() + 
  scale_fill_manual(values=c("#999999","#56B4E9")) + 
  scale_x_continuous(labels = pal3.first20$base,breaks=pal3.first20$index)
ggsave(gg,
       file="/Users/aaronmck/Desktop/code/base_editing/2019_12_22_Maryam_base_editing/plots/editing_AT_GC_split_pal3.png",
       width=12,
       height=8)

