setwd('/Users/eliselucotte/desktop/server1/Assembly/ortholog_reconstruction/Scorzo/dSplot/')
library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)
library(stringr)
#####################################################################################################

d2=read.table('database_for_dSplot_M-scorzo_withrank.txt', header=T, as.is=T)

### MAKING THE PLOTS  
## for each chromosome, 3 plots : the dS plot, the gene rank for A1 and the gene rank for A2

#definition of the coefficient for the rank
coeff=1000
size=20
MC12=ggplot(data=d2[which(d2$chrom_lag=='Mv-lag-1253-A1_MC12'),])+
  geom_point(aes(x=new_axis,y=dS,color=as.factor(colors)))+
  scale_colour_identity()+
  scale_x_continuous(breaks=seq(0,max(d2[which(d2$chrom_lag=='Mv-lag-1253-A1_MC12'),]$new_axis), by = 0.1))+
  theme_bw(base_size=size)+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ylim(0,0.3)+
  ggtitle('MC12')+
  xlab('')+
  ylab('')+
  geom_vline(aes(xintercept= 343080/10^6), col='red')

MC12_rankA1=ggplot(data=d2[which(d2$chrom_lag=='Mv-lag-1253-A1_MC12'),])+
  geom_point(aes(x=new_axis,y=rankA1/coeff,color=chrom1))+
  scale_x_continuous(breaks=seq(0,max(d2[which(d2$chrom_lag=='Mv-lag-1253-A1_MC12'),]$new_axis), by = 0.1))+
  theme_bw(base_size=size)+
  #scale_colour_identity()+
  ylim(min(d2$rankA1/1000),max(d2$rankA1/1000))+
  theme(legend.position = "none",axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.text.x = element_text(angle = 45, hjust=1))+
  xlab('')+
  ylab('')

MC12_rankA2=ggplot(data=d2[which(d2$chrom_lag=='Mv-lag-1253-A1_MC12'),])+
  geom_point(aes(x=new_axis,y=rankA2/coeff,color=chrom2))+
  scale_x_continuous(breaks=seq(0,max(d2[which(d2$chrom_lag=='Mv-lag-1253-A1_MC12'),]$new_axis), by = 0.1))+
  theme_bw(base_size=size)+
  #scale_colour_identity()+
  ylim(min(d2$rankA2/1000),max(d2$rankA2/1000))+
  theme(legend.position = "none",axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.text.x = element_text(angle = 45, hjust=1))+
  xlab('position (Mb) in M. Lagerheimii ')+
  ylab("")

MC03A=ggplot(data=d2[which(d2$chrom_lag=='Mv-lag-1253-A1_MC03A'),], aes(x=new_axis,y=dS))+
  geom_point(aes(color=as.factor(colors)))+
  scale_colour_identity()+
  scale_x_continuous(breaks=seq(0,max(d2[which(d2$chrom_lag=='Mv-lag-1253-A1_MC03A'),]$new_axis), by = 0.1))+
  theme_bw(base_size=size)+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ylim(0,0.3)+
  #xlim(min(d$start_lag),800038)+
  geom_vline(aes(xintercept= 282495/10^6), col='red')+
  xlab('')+
  ggtitle('MC03A')

MC03A_rankA1=ggplot(data=d2[which(d2$chrom_lag=='Mv-lag-1253-A1_MC03A'),], aes(x=new_axis,y=dS))+
  geom_point(aes(x=new_axis,y=rankA1/coeff,color=chrom1))+
  scale_x_continuous(breaks=seq(0,max(d2[which(d2$chrom_lag=='Mv-lag-1253-A1_MC03A'),]$new_axis), by = 0.1))+
  theme_bw(base_size=size)+
  ylim(min(d2$rankA1/1000),max(d2$rankA1/1000))+
  #scale_colour_identity()+
  theme(legend.position = "none",axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.text.x = element_text(angle = 45, hjust=1))+
  xlab('')+
  ylab(bquote(a[1]*" gene rank"))

MC03A_rankA2=ggplot(data=d2[which(d2$chrom_lag=='Mv-lag-1253-A1_MC03A'),], aes(x=new_axis,y=dS))+
  geom_point(aes(x=new_axis,y=rankA2/coeff,color=chrom2))+
  scale_x_continuous(breaks=seq(0,max(d2[which(d2$chrom_lag=='Mv-lag-1253-A1_MC03A'),]$new_axis), by = 0.1))+
  theme_bw(base_size=size)+
  #scale_colour_identity()+
  ylim(min(d2$rankA2/1000),max(d2$rankA2/1000))+
  theme(legend.position = "none",axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.text.x = element_text(angle = 45, hjust=1))+
  xlab('')+
  ylab(bquote(a[2]*" gene rank"))


MC03=ggplot(data=d2[which(d2$chrom_lag=='Mv-lag-1253-A1_MC03B'),], aes(x=new_axis,y=dS))+
  geom_point(aes(color=as.factor(colors)))+
  scale_colour_identity()+
  scale_x_continuous(breaks=seq(800000/10^6,max(d2[which(d2$chrom_lag=='Mv-lag-1253-A1_MC03B'),]$new_axis), by = 0.1))+
  theme_bw(base_size=size)+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ylim(0,0.3)+
  #xlim(800038,max(d$start_lag))+
  xlab('')+
  ylab('')+
  ggtitle('MC03B')

MC03_rankA1=ggplot(data=d2[which(d2$chrom_lag=='Mv-lag-1253-A1_MC03B'),], aes(x=new_axis,y=dS))+
  geom_point(aes(x=new_axis,y=rankA1/coeff,color=chrom1))+
  scale_x_continuous(breaks=seq(0,max(d2[which(d2$chrom_lag=='Mv-lag-1253-A1_MC03B'),]$new_axis), by = 0.1))+
  theme_bw(base_size=size)+
  ylim(min(d2$rankA1/1000),max(d2$rankA1/1000))+
  #scale_colour_identity()+
  theme(legend.position = "none",axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.text.x = element_text(angle = 45, hjust=1))+
  xlab('')+
  ylab('')

MC03_rankA2=ggplot(data=d2[which(d2$chrom_lag=='Mv-lag-1253-A1_MC03B'),], aes(x=new_axis,y=dS))+
  geom_point(aes(x=new_axis,y=rankA2/coeff,color=chrom2))+
  scale_x_continuous(breaks=seq(0,max(d2[which(d2$chrom_lag=='Mv-lag-1253-A1_MC03B'),]$new_axis), by = 0.1))+
  theme_bw(base_size=size)+
  #scale_colour_identity()+
  ylim(min(d2$rankA2/1000),max(d2$rankA2/1000))+
  theme(legend.position = "none",axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.text.x = element_text(angle = 45, hjust=1))+
  xlab('')+
  ylab('')

MC16=ggplot(data=d2[which(d2$chrom_lag=='Mv-lag-1253-A1_MC16'),], aes(x=new_axis,y=dS))+
  geom_point(aes(color=as.factor(colors)))+
  scale_colour_identity()+
  scale_x_continuous(breaks=seq(0,max(d2[which(d2$chrom_lag=='Mv-lag-1253-A1_MC16'),]$new_axis), by = 0.1))+
  ylim(0,0.3)+
  ggtitle('MC16')+
  xlab('')+
  ylab('')+
  theme_bw(base_size=size)+
  theme(axis.text.x = element_text(angle = 45, hjust=1))


MC16_rankA1=ggplot(data=d2[which(d2$chrom_lag=='Mv-lag-1253-A1_MC16'),], aes(x=new_axis,y=dS))+
  geom_point(aes(x=new_axis,y=rankA1/coeff,color=chrom1))+
  scale_x_continuous(breaks=seq(0,max(d2[which(d2$chrom_lag=='Mv-lag-1253-A1_MC16'),]$new_axis), by = 0.1))+
  #scale_colour_identity()+
  theme_bw(base_size=size)+
  ylim(min(d2$rankA1/1000),max(d2$rankA1/1000))+
  theme(legend.position = "none",axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.text.x = element_text(angle = 45, hjust=1))+
  xlab('')+
  ylab('')

MC16_rankA2=ggplot(data=d2[which(d2$chrom_lag=='Mv-lag-1253-A1_MC16'),], aes(x=new_axis,y=dS))+
  geom_point(aes(x=new_axis,y=rankA2/coeff,color=chrom2))+
  scale_x_continuous(breaks=seq(0,max(d2[which(d2$chrom_lag=='Mv-lag-1253-A1_MC16'),]$new_axis), by = 0.1))+
  #scale_colour_identity()+
  theme_bw(base_size=size)+
  ylim(min(d2$rankA2/1000),max(d2$rankA2/1000))+
  theme(legend.position = "none",axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.text.x = element_text(angle = 45, hjust=1))+
  xlab('')+
  ylab('')

## Making the layout
## the width was calculated proportionally to the length of the chromosome so that the sum is = 1
MC03A+MC03+MC12+MC16+MC03A_rankA1+MC03_rankA1+MC12_rankA1+MC16_rankA1+
  MC03A_rankA2+MC03_rankA2+MC12_rankA2+MC16_rankA2+
  plot_layout(ncol = 4, nrow=3,width=c(0.21,0.365,0.33,0.095),heights = c(1,0.3,0.3))


ggsave("dS_plot_Mscorzo-nooutlier_strata_with_rank2.png", height=10, width=15)
ggsave("dS_plot_Mscorzo-nooutlier_strata_with_rank2.pdf", device='pdf', height=10, width=15)





