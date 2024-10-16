library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
setwd('/Users/eliselucotte/desktop/server1/Assembly/ortholog_reconstruction/')
#R –arch x86_64 CMD INSTALL rjags_4-4.tar.gz 
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

size=18
MCP=plot(fit_mcp) + ylab(expression(italic(dS))) + theme_bw(base_size=size) +
  xlab('')+
  ylim(-1,0.2)+
  geom_point(size = 1, shape = 3) 


MC12=ggplot(data=d2[which(d2$chrom_lag %in% PR_chrom),])+
  geom_point(aes(x=rank_lag,y=dS,color=as.factor(colors)))+
  scale_colour_identity()+
  #scale_x_continuous(breaks=seq(0,max(d2[which(d2$chrom_lag=='Mv-lag-1253-A1_MC12'),]$rank_lag), by = 0.1))+
  theme_bw(base_size=size)+
  ylim(0,0.2)+
  ggtitle('PR')+
  xlab('')+
  ylab(expression(italic(dS)))
#+
#geom_vline(aes(xintercept= 343080/10^6), col='red')

MC12_rankA1=ggplot(data=d2[which(d2$chrom_lag %in% PR_chrom),])+
  geom_point(aes(x=rank_lag,y=rankA1/coeff,color=colors_chromA1))+
  #scale_x_continuous(breaks=seq(0,max(d2[which(d2$chrom_lag=='Mv-lag-1253-A1_MC12'),]$rank_lag), by = 0.1))+
  theme_bw(base_size=size)+
  scale_colour_identity()+
  ylim(min(d2$rankA1/coeff),max(d2$rankA1/coeff))+
  theme(legend.position = "none",axis.text.y=element_blank(),axis.ticks.y=element_blank())+
  xlab('')+
  ylab(bquote(a[1]*" gene rank"))

MC12_rankA2=ggplot(data=d2[which(d2$chrom_lag %in% PR_chrom),])+
  geom_point(aes(x=rank_lag,y=rankA2/coeff,color=colors_chromA2))+
  #scale_x_continuous(breaks=seq(0,max(d2[which(d2$chrom_lag=='Mv-lag-1253-A1_MC12'),]$rank_lag), by = 0.1))+
  theme_bw(base_size=size)+
  scale_colour_identity()+
  ylim(min(d2$rankA2/coeff),max(d2$rankA2/coeff))+
  theme(legend.position = "none",axis.text.y=element_blank(),axis.ticks.y=element_blank())+
  xlab('order along the M. Lagerheimii genome')+
  ylab(bquote(a[2]*"\n gene rank"))

MC12+MCP+MC12_rankA1+MC12_rankA2+plot_layout(ncol=1,heights = c(1,1,0.4,0.4))

ggsave("./2.make_DS_plot/Mv-cat-C212/dS_plot_Mvsup-1065-PR_with_model.png", height=10, width=10)
ggsave("./2.make_DS_plot/Mv-cat-C212/dS_plot_Mvsup-1065-PR_with_model.pdf", device='pdf', height=10, width=10)

pdf(file = './2.make_DS_plot/Mv-cat-C212/fit_mcp.pdf', width=10, height=25)
plot_pars(fit_mcp)
dev.off()
write.table(fitted(fit_mcp), file='./2.make_DS_plot/Mv-cat-C212/fitted_fit_mcp.txt')
predict(fit_mcp)
write.table(predict(fit_mcp), file='./2.make_DS_plot/Mv-cat-C212/predict_fit_mcp.txt')
