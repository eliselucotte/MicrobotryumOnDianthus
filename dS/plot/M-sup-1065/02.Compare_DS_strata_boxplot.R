library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(ggpubr)
setwd('/Users/eliselucotte/desktop/server0/Assembly/ortholog_reconstruction/')
d2=read.table('./2.make_DS_plot/Mv-sup-1065/database_for_dSplot_Mv-sup-1065_withrank_new_color.txt', header=T, as.is=T)
d2=d2[which(d2$dS<1),]
d2=d2[-which(d2$strata_name=='in_HD'),]
d2=d2%>%mutate(new_color=ifelse(colors=='black',ifelse(chrom_lag=='Mv-lag-1253-A1_MC03A','black_HD','black'),ifelse(strata_name=='PARs' & chrom_lag=='Mv-lag-1253-A1_MC03B', 'PARs2',ifelse(strata_name=='PARs' & chrom_lag=='Mv-lag-1253-A1_MC12','PARs1',strata_name))))
d2=d2%>%mutate(color2=ifelse(new_color=='black_HD','grey',colors))

### MAKING THE PLOTS  
## for each chromosome, 3 plots : the dS plot, the gene rank for A1 and the gene rank for A2

#definition of the coefficient for the rank
size=20
p<-ggplot(data=d2, aes(y=dS, x=new_color, fill=color2), colour="black")+
  geom_violin(width=0.8, scale="width")+
  geom_boxplot(fill='white',width=0.1)+
  geom_jitter(width = 0.2, pch=21)+
  theme_bw(base_size=size)+
  scale_fill_identity()+
  scale_colour_identity()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  scale_x_discrete(limits=c('blue','red','HD_RR','PARs1','light_ruby','dark_ruby','orange','purple','black','turquoise','PARs2'),labels=c('blue'='blue','red'='HD genes','HD_RR'='PARs HD','PARs1'='PARs1','light_ruby'='light ruby','dark_ruby'='dark ruby','orange'='orange','purple'='purple','black'='black','turquoise'='turquoise','PARs2'='PARs2'))+
  xlab('Strata')+
  ylab("dS")

compare_means(dS ~ new_color, data = d2, paired = FALSE)
my_comparisons <- list(c("PARs1", "light_ruby"), c("light_ruby", "dark_ruby"), c("black", "dark_ruby"),c('black','turquoise'),c('turquoise','PARs2'))
p+stat_compare_means(comparisons = my_comparisons)
ggsave("./2.make_DS_plot/Mv-sup-1065/Boxplot_strata_dS.png", height=6, width=10)
ggsave("./2.make_DS_plot/Mv-sup-1065/Boxplot_strata_dS.pdf", device='pdf', height=6, width=10)

## NO HD GENES
#definition of the coefficient for the rank
size=20
p<-ggplot(data=d2, aes(y=dS, x=new_color, fill=color2), colour="black")+
  geom_violin(width=0.8, scale="width", alpha=0.6)+
  geom_jitter(width = 0.3, pch=21)+
  geom_boxplot(fill='white',width=0.1,alpha=0.5)+
  theme_bw(base_size=size)+
  scale_fill_identity()+
  scale_colour_identity()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  scale_x_discrete(limits=c('blue','HD_RR','PARs1','light_ruby','dark_ruby','orange','purple','black','turquoise','PARs2'),labels=c('blue'='blue','HD_RR'='PARs HD','PARs1'='PARs1','light_ruby'='light ruby','dark_ruby'='dark ruby','orange'='orange','purple'='purple','black'='black','turquoise'='turquoise','PARs2'='PARs2'))+
  xlab('Strata')+
  ylab("dS")
compare_means(dS ~ new_color, data = d2, paired = FALSE)
my_comparisons <- list(c("PARs1", "light_ruby"), c("light_ruby", "dark_ruby"), c("black", "dark_ruby"),c('black','turquoise'),c('turquoise','PARs2'))
p+stat_compare_means(comparisons = my_comparisons)
ggsave("./2.make_DS_plot/Mv-sup-1065/Boxplot_strata_dS_noHD.png", height=6, width=15)
ggsave("./2.make_DS_plot/Mv-sup-1065/Boxplot_strata_dS_noHD.pdf", device='pdf', height=6, width=15)

