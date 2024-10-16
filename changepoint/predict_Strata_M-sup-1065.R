library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
setwd('/Users/eliselucotte/desktop/server0/Assembly/ortholog_reconstruction/')
#R –arch x86_64 CMD INSTALL rjags_4-4.tar.gz 
d2=read.table('./2.make_DS_plot/Mv-sup-1065/database_for_dSplot_Mv-sup-1065_withrank_new_color.txt', header=T, as.is=T)

#remove the genes that are now located on the HD chromosomes
d2=d2[-which(d2$strata_name=='in_HD'),]
PR_chrom=c('Mv-lag-1253-A1_MC12','Mv-lag-1253-A1_MC03B','Mv-lag-1253-A1_MC16')

library(mcp)
library(magrittr)
library(rjags)
## change the number of 1~1 according to the number of strata you think there are
model = list(dS ~ 1, 1~ 1, 1 ~ 1, 1~1)

df=select(d2[which(d2$chrom_lag %in% PR_chrom),], c(rank_lag,dS))
df=df[order(df$rank_lag),]
head(df)
df=df[-which(df$rank_lag==1),]
head(df)
#df est un dataframe avec deux colonnes: 
# 1 - les positions le long du chromosome (ici: order)
# 2 - la valeur de Ds
df = df %>% filter(df$dS <0.2 )
colnames(df)=c('ordre','Ds')
write.table(df, "./2.make_DS_plot/find_strata/Mv-sup-1065/Table_for_changepoint.txt", quote=F, row.names = F)
colnames(df)=c('rank_lag','dS')

fit_mcp = mcp(model, data = df, par_x = "rank_lag")
summary(fit_mcp)

coeff=1000
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
  ylim(0,0.3)+
  ggtitle('PR')+
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


MC12+MCP+MC12_rankA1+MC12_rankA2+plot_layout(ncol=1,heights = c(1,1,0.4,0.4))

ggsave("./2.make_DS_plot/Mv-sup-1065/with_extra_contig/dS_plot_Mvsup-1065-PR_with_model.png", height=10, width=10)
ggsave("./2.make_DS_plot/Mv-sup-1065/with_extra_contig/dS_plot_Mvsup-1065-PR_with_model.pdf", device='pdf', height=10, width=10)

pdf(file = './2.make_DS_plot/Mv-sup-1065/with_extra_contig/fit_mcp.pdf', width=10, height=25)
plot_pars(fit_mcp)
dev.off()
write.table(fitted(fit_mcp), file='./2.make_DS_plot/Mv-sup-1065/with_extra_contig/fitted_fit_mcp.txt')
predict(fit_mcp)
write.table(predict(fit_mcp), file='./2.make_DS_plot/Mv-sup-1065/with_extra_contig/predict_fit_mcp.txt')
