library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
setwd('/Users/eliselucotte/desktop/server1/Assembly/ortholog_reconstruction/')
d2=read.table('./2.make_DS_plot/Mv-sup-1065/database_for_dSplot_Mv-sup-1065_withrank_new_color.txt', header=T, as.is=T)

### MAKING THE PLOTS  
## for each chromosome, 3 plots : the dS plot, the gene rank for A1 and the gene rank for A2

#definition of the coefficient for the rank
coeff=1000
size_text=20
PR_chrom=c('Mv-lag-1253-A1_MC12','Mv-lag-1253-A1_MC03B','Mv-lag-1253-A1_MC16')
MC12=ggplot(data=d2[which(d2$chrom_lag %in% PR_chrom),])+
  geom_point(aes(x=rank_lag,y=dS,color=as.factor(colors)))+
  scale_colour_identity()+
  #scale_x_continuous(breaks=seq(0,max(d2[which(d2$chrom_lag=='Mv-lag-1253-A1_MC12'),]$rank_lag), by = 0.1))+
  theme_bw()+
  theme(text = element_text(size = size_text))+
  ylim(0,0.3)+
  ggtitle('PR')+
  xlab('')+
  ylab('')+
  geom_vline(aes(xintercept= 64), col='red')

MC12_rankA1=ggplot(data=d2[which(d2$chrom_lag %in% PR_chrom),])+
  geom_point(aes(x=rank_lag,y=rankA1/coeff,color=colors_chromA1))+
  #scale_x_continuous(breaks=seq(0,max(d2[which(d2$chrom_lag=='Mv-lag-1253-A1_MC12'),]$rank_lag), by = 0.1))+
  theme_bw()+
  theme(text = element_text(size = size_text))+
  scale_colour_identity()+
  ylim(min(d2$rankA1/coeff),max(d2$rankA1/coeff))+
  theme(legend.position = "none",axis.text.y=element_blank(),axis.ticks.y=element_blank())+
  xlab('')+
  ylab('')

MC12_rankA2=ggplot(data=d2[which(d2$chrom_lag %in% PR_chrom),])+
  geom_point(aes(x=rank_lag,y=rankA2/coeff,color=colors_chromA2))+
  #scale_x_continuous(breaks=seq(0,max(d2[which(d2$chrom_lag=='Mv-lag-1253-A1_MC12'),]$rank_lag), by = 0.1))+
  theme_bw()+
  theme(text = element_text(size = size_text))+
  scale_colour_identity()+
  ylim(min(d2$rankA2/coeff),max(d2$rankA2/coeff))+
  theme(legend.position = "none",axis.text.y=element_blank(),axis.ticks.y=element_blank())+
  xlab('order along the M. Lagerheimii genome')+
  ylab("")

MC03A=ggplot(data=d2[which(d2$chrom_lag=='Mv-lag-1253-A1_MC03A'),], aes(x=rank_lag,y=dS))+
  geom_point(aes(color=as.factor(colors)))+
  scale_colour_identity()+
  #scale_x_continuous(breaks=seq(0,max(d2[which(d2$chrom_lag=='Mv-lag-1253-A1_MC03A'),]$rank_lag), by = 0.1))+
  theme_bw()+
  theme(text = element_text(size = size_text))+
  ylim(0,0.3)+
  geom_vline(aes(xintercept= 48), col='red')+
  xlab('')+
  ggtitle('HD')

MC03A_rankA1=ggplot(data=d2[which(d2$chrom_lag=='Mv-lag-1253-A1_MC03A'),], aes(x=rank_lag,y=dS))+
  geom_point(aes(x=rank_lag,y=rankA1/coeff,color=colors_chromA1))+
  #scale_x_continuous(breaks=seq(0,max(d2[which(d2$chrom_lag=='Mv-lag-1253-A1_MC03A'),]$rank_lag), by = 0.1))+
  theme_bw()+
  theme(text = element_text(size = size_text))+
  ylim(min(d2$rankA1/coeff),max(d2$rankA1/coeff))+
  scale_colour_identity()+
  theme(legend.position = "none",axis.text.y=element_blank(),axis.ticks.y=element_blank())+
  xlab('')+
  ylab(bquote(a[1]*" gene rank"))

MC03A_rankA2=ggplot(data=d2[which(d2$chrom_lag=='Mv-lag-1253-A1_MC03A'),], aes(x=rank_lag,y=dS))+
  geom_point(aes(x=rank_lag,y=rankA2/coeff,color=colors_chromA2))+
  #scale_x_continuous(breaks=seq(0,max(d2[which(d2$chrom_lag=='Mv-lag-1253-A1_MC03A'),]$rank_lag), by = 0.1))+
  theme_bw()+
  scale_colour_identity()+
  theme(text = element_text(size = size_text))+
  ylim(min(d2$rankA2/coeff),max(d2$rankA2/coeff))+
  theme(legend.position = "none",axis.text.y=element_blank(),axis.ticks.y=element_blank())+
  xlab('')+
  ylab(bquote(a[2]*" gene rank"))


## Making the layout
## the width was calculated proportionally to the length of the chromosome so that the sum is = 1
MC03A+MC12+MC03A_rankA1+MC12_rankA1+
  MC03A_rankA2+MC12_rankA2+
  plot_layout(ncol = 2, nrow=3,width=c(0.21,0.79),heights = c(1,0.3,0.3))


ggsave("./2.make_DS_plot/Mv-sup-1065/dS_plot_Mvsup-1065-nooutlier_strata_with_ranklag.png", height=10, width=15)
ggsave("./2.make_DS_plot/Mv-sup-1065/dS_plot_Mvsup-1065-nooutlier_strata_with_ranklag.pdf", device='pdf', height=10, width=15)
