library(ggplot2)
library(data.table)
library(wesanderson)

plot_loc = "/Users/aaronmck/Desktop/code/base_editing/2019_12_12_control_data/plots/"

tbl_control = fread("/Users/aaronmck/Desktop/code/base_editing/2019_12_12_control_data/notebooks/palindrome_sequence_controls_UMI.txt")

tbl_control = tbl_control[order(tbl_control$percent_correct_target,decreasing = T),]
tbl_control$seq = seq(1,nrow(tbl_control))
tbl_control$number_of_other_targets = (tbl_control$number_of_other_targets / 314) * 100.0
tbl.control.melt = melt(tbl_control,measure.vars = c("percent_correct_target","percent_other_target","number_of_other_targets"))

pal <- wes_palette("Zissou1", 3, type = "continuous")
gg = ggplot(tbl.control.melt) + 
  geom_bar(aes(x=seq,y=value,fill=variable),stat = 'identity',position='dodge') + 
  theme_classic() +
  ylim(c(0,100)) +scale_fill_manual(values=pal) + facet_wrap(. ~ variable)+
  xlab("Guides, ordered by matched target percentage") + 
  ylab("Percent of targets that match the guide")
ggsave(gg,file=paste(plot_loc,"base_editing_proportion_control_UMI.png",sep=""),width=12, height=4)

# cases

tbl_control = fread("/Users/aaronmck/Desktop/code/base_editing/2019_12_12_control_data/notebooks/palindrome_sequence_cases_UMI.txt")

tbl_control = tbl_control[order(tbl_control$percent_correct_target,decreasing = T),]
tbl_control$seq = seq(1,nrow(tbl_control))
tbl_control$number_of_other_targets = (tbl_control$number_of_other_targets / 314) * 100.0
tbl.control.melt = melt(tbl_control,measure.vars = c("percent_correct_target","percent_other_target","number_of_other_targets"))

pal <- wes_palette("Zissou1", 3, type = "continuous")
gg = ggplot(tbl.control.melt) + 
  geom_bar(aes(x=seq,y=value,fill=variable),stat = 'identity',position='dodge') + 
  theme_classic() +
  ylim(c(0,100)) +scale_fill_manual(values=pal) + facet_wrap(. ~ variable) +
  xlab("Guides, ordered by matched target percentage") + 
  ylab("Percent of targets that match the guide")
ggsave(gg,file=paste(plot_loc,"base_editing_proportion_cases_UMI.png",sep=""),width=12, height=4)

# load up the edited and not per-target
control = "/Users/aaronmck/Desktop/code/base_editing/2019_12_12_control_data/notebooks/palindromes_base_counts_control_UMI.txt"
edited = "/Users/aaronmck/Desktop/code/base_editing/2019_12_12_control_data/notebooks/palindromes_base_counts_edited_UMI.txt"
control.counts = fread(control)
edited.counts = fread(edited)

# melted.control.counts = melt(control.counts,measure.vars = c("base1","base2","base3","base4","base5","base6","base7","base8","base9","base10",
#                                                             "base11","base12","base13","base14","base15","base16","base17","base18","base19","base20"))

# from http://www.sthda.com/english/wiki/ggplot2-error-bars-quick-start-guide-r-software-and-data-visualization
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

# process the controls
melted.control.counts = melt(control.counts[control.counts$target_count > 10,],measure.vars = c("1","2","3","4","5","6","7","8","9","10",
                                                                                                "11","12","13","14","15","16","17","18","19","20"))

melted.control.counts.sumary = data_summary(melted.control.counts,"value","variable")
p<- ggplot(melted.control.counts.sumary, aes(x=variable, y=value)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2,
                position=position_dodge(.9)) + 
  theme_classic() +
  xlab("Target position") +
  ylab("proportion edited (with error bars)")
ggsave(p,file=paste(plot_loc,"control_counts_editing_rate_UMI_0.1.png",sep=""),width=8,height=3)


# and the cases
melted.edited.counts = melt(edited.counts[edited.counts$target_count > 10,],measure.vars = c("1","2","3","4","5","6","7","8","9","10",
                                                                                             "11","12","13","14","15","16","17","18","19","20"))

melted.edited.counts.sumary = data_summary(melted.edited.counts,"value","variable")
p<- ggplot(melted.edited.counts.sumary, aes(x=variable, y=value)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2,
                position=position_dodge(.9)) + 
  theme_classic() +
  xlab("Target position") +
  ylab("proportion edited (with error bars)")
ggsave(p,file=paste(plot_loc,"edited_counts_editing_rate_UMI_0.1.png",sep=""),width=8,height=3)