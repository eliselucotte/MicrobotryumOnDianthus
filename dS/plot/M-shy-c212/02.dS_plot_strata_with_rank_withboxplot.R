library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
setwd('/Users/eliselucotte/desktop/server1/Assembly/ortholog_reconstruction/')
d2=read.table('./2.make_DS_plot/Mv-cat-C212/database_for_dSplot_Mv-cat-C212_withrank.txt', header=T, as.is=T)

### MAKING THE PLOTS  
## for each chromosome, 3 plots : the dS plot, the gene rank for A1 and the gene rank for A2

#definition of the coefficient for the rank
coeff=1000
size=20

PR_chrom=c('Mv-lag-1253-A1_MC12','Mv-lag-1253-A1_MC03B','Mv-lag-1253-A1_MC16')
MC12=ggplot(data=d2[which(d2$chrom_lag %in% PR_chrom),])+
  geom_point(aes(x=rank_lag,y=dS,color=as.factor(colors)))+
  scale_colour_identity()+
  #scale_x_continuous(breaks=seq(0,max(d2[which(d2$chrom_lag=='Mv-lag-1253-A1_MC12'),]$rank_lag), by = 0.1))+
  theme_bw(base_size=size)+
  ylim(0,0.3)+
  ggtitle('PR')+
  geom_segment(x = 5, y = -0.01, xend = 18,yend = -0.01, color="indianred1")+
  geom_segment(x = 19, y = -0.01, xend = 50,yend = -0.01, color="darkred")+
  geom_segment(x = 205, y = -0.01, xend = 232,yend = -0.01, color="turquoise")+
  xlab('')+
  ylab('')+
  geom_vline(aes(xintercept= 64), col='red')

MC12_rankA1=ggplot(data=d2[which(d2$chrom_lag %in% PR_chrom),])+
  geom_point(aes(x=rank_lag,y=rankA1/coeff,color=colors_chromA1))+
  #scale_x_continuous(breaks=seq(0,max(d2[which(d2$chrom_lag=='Mv-lag-1253-A1_MC12'),]$rank_lag), by = 0.1))+
  theme_bw(base_size=size)+
  scale_colour_identity()+
  ylim(min(d2$rankA1/coeff),max(d2$rankA1/coeff))+
  theme(legend.position = "none",axis.text.y=element_blank(),axis.ticks.y=element_blank())+
  xlab('')+
  ylab('')

MC12_rankA2=ggplot(data=d2[which(d2$chrom_lag %in% PR_chrom),])+
  geom_point(aes(x=rank_lag,y=rankA2/coeff,color=colors_chromA2))+
  #scale_x_continuous(breaks=seq(0,max(d2[which(d2$chrom_lag=='Mv-lag-1253-A1_MC12'),]$rank_lag), by = 0.1))+
  theme_bw(base_size=size)+
  scale_colour_identity()+
  ylim(min(d2$rankA2/coeff),max(d2$rankA2/coeff))+
  theme(legend.position = "none",axis.text.y=element_blank(),axis.ticks.y=element_blank())+
  xlab('order along the M. Lagerheimii genome')+
  ylab("")

MC03A=ggplot(data=d2[which(d2$chrom_lag=='Mv-lag-1253-A1_MC03A'),], aes(x=rank_lag,y=dS))+
  geom_point(aes(color=as.factor(colors)))+
  scale_colour_identity()+
  #scale_x_continuous(breaks=seq(0,max(d2[which(d2$chrom_lag=='Mv-lag-1253-A1_MC03A'),]$rank_lag), by = 0.1))+
  theme_bw(base_size=size)+
  ylim(0,0.3)+
  #xlim(min(d2$rank_lag),)+
  geom_vline(aes(xintercept= 11), col='red')+
  xlab('')+
  ggtitle('HD')

MC03A_rankA1=ggplot(data=d2[which(d2$chrom_lag=='Mv-lag-1253-A1_MC03A'),], aes(x=rank_lag,y=dS))+
  geom_point(aes(x=rank_lag,y=rankA1/coeff,color=colors_chromA1))+
  #scale_x_continuous(breaks=seq(0,max(d2[which(d2$chrom_lag=='Mv-lag-1253-A1_MC03A'),]$rank_lag), by = 0.1))+
  theme_bw(base_size=size)+
  ylim(min(d2$rankA1/coeff),max(d2$rankA1/coeff))+
  scale_colour_identity()+
  theme(legend.position = "none",axis.text.y=element_blank(),axis.ticks.y=element_blank())+
  xlab('')+
  ylab(bquote(a[1]))

MC03A_rankA2=ggplot(data=d2[which(d2$chrom_lag=='Mv-lag-1253-A1_MC03A'),], aes(x=rank_lag,y=dS))+
  geom_point(aes(x=rank_lag,y=rankA2/coeff,color=colors_chromA2))+
  #scale_x_continuous(breaks=seq(0,max(d2[which(d2$chrom_lag=='Mv-lag-1253-A1_MC03A'),]$rank_lag), by = 0.1))+
  theme_bw(base_size=size)+
  scale_colour_identity()+
  ylim(min(d2$rankA2/coeff),max(d2$rankA2/coeff))+
  theme(legend.position = "none",axis.text.y=element_blank(),axis.ticks.y=element_blank())+
  xlab('')+
  ylab(bquote(a[2]))


d2=read.table('./2.make_DS_plot/Mv-cat-C212/database_for_dSplot_Mv-cat-C212_withrank_new_color.txt', header=T, as.is=T)
d2=d2[-which(d2$strata_name=='in_HD'),]

PR_chrom=c('Mv-lag-1253-A1_MC12','Mv-lag-1253-A1_MC03B','Mv-lag-1253-A1_MC16')

library(mcp)
library(magrittr)
library(rjags)
## change the number of 1~1 according to the number of strata you think there are
model = list(dS ~ 1, 1~ 1, 1 ~ 1, 1~1)

df=select(d2[which(d2$chrom_lag %in% PR_chrom),], c(rank_lag,dS))
#df est un dataframe avec deux colonnes: 
# 1 - les positions le long du chromosome (ici: order)
# 2 - la valeur de Ds
df = df %>% filter(df$dS <0.3 )

fit_mcp = mcp(model, data = df, par_x = "rank_lag")
summary(fit_mcp)

size=20
MCP=plot(fit_mcp) + ylab(expression(italic(dS))) + theme_bw(base_size=size) +
  xlab('')+
  ylim(-1,0.3)+
  geom_point(size = 1, shape = 3) 

MC03A+MC12+plot_spacer()+MCP+MC03A_rankA1+MC12_rankA1+
  MC03A_rankA2+MC12_rankA2+
  plot_layout(ncol = 2, nrow=4,width=c(0.21,0.79),heights = c(1,1,0.3,0.3))#+
plot_annotation(tag_levels = 'A')

ggsave("./2.make_DS_plot/Mv-cat-C212/Figure_dS_changepoint.png", height=10, width=13)
ggsave("./2.make_DS_plot/Mv-cat-C212/Figure_dS_changepoint.pdf", device='pdf', height=10, width=13)


