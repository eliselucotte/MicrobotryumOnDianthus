library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
MAT='A1'
window=1000
size=18
path="/Users/eliselucotte/desktop/server1/Assembly/sexchrom/detect_telomeres/"
species='Mv-cat-C212'
dchrom <- c(
  "h1tg000040l"='HD',
  "h1tg000006l"='PR',
  "h2tg000018l"='HD',
  "h2tg000005l"="PR"
)
# path="/Users/eliselucotte/desktop/server1/Assembly/sexchrom/detect_telomeres/Mv-sup/"
# species="MvDp-1065"
# dchrom <- c(
#   "tig00000040"='HD',
#   "tig00000080"='PR1',
#   "tig00000084"='PR',
#   "tig00000031"='HD',
#   "tig00000068"="PR"
# )



t=read.table(paste(path,species,'-',MAT,'-count_telomeres_motifs_',window,'bp.txt',sep=''), as.is=T, header=T )
t$interval=t$interval/window
t=mutate(t,contig2=dchrom[t$contig])
#positions of the telomeres
View(t[which(t$contig=='h2tg000005l'),])
#13275*1000

A1=ggplot(data=t, aes(y=count, x=interval))+
  geom_line()+
  theme_bw(base_size = size)+
  ylim(c(0,17))+
  labs(title='A1 TTAGGG',x='position along the contig', y='motif count')+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  facet_wrap(~contig2, scales='free_x')

t=read.table(paste(path,species,'-',MAT,'-count_telomeres_motifs_inv_compl_',window,'bp.txt',sep=''), as.is=T, header=T )
t$interval=t$interval/window
t=mutate(t,contig2=dchrom[t$contig])
#positions of the telomeres
View(t[which(t$contig=='h2tg000005l'),])
# 37*10000

A1_inv_compl=ggplot(data=t, aes(y=count, x=interval))+
  geom_line()+
  theme_bw(base_size = size)+
  ylim(c(0,17))+
  labs(title='A1 CCCTAA',x='position along the contig', y='motif count')+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  facet_wrap(~contig2, scales='free_x')

MAT='A2'
t=read.table(paste(path,species,'-',MAT,'-count_telomeres_motifs_',window,'bp.txt',sep=''), as.is=T, header=T )
t$interval=t$interval/window
t=mutate(t,contig2=dchrom[t$contig])

A2=ggplot(data=t, aes(y=count, x=interval))+
  geom_line()+
  ylim(c(0,17))+
  theme_bw(base_size = size)+
  labs(title='A2 TTAGGG',x='position along the contig', y='motif count')+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  facet_wrap(~contig2, scales='free_x')

t=read.table(paste(path,species,'-',MAT,'-count_telomeres_motifs_inv_compl_',window,'bp.txt',sep=''), as.is=T, header=T )
t$interval=t$interval/window
t=mutate(t,contig2=dchrom[t$contig])

A2_inv_compl=ggplot(data=t, aes(y=count, x=interval))+
  geom_line()+
  theme_bw(base_size = size)+
  ylim(c(0,17))+
  labs(title='A2 CCCTAA',x='position along the contig', y='motif count')+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  facet_wrap(~contig2, scales='free_x')


p=A1+A2+A1_inv_compl+A2_inv_compl+plot_layout(ncol = 2, nrow=2)#,width=c(0.6,0.4))
p
ggsave(paste(path,species,'-count_telomeres_motifs_all_',window,'bp.png',sep=''),plot=p,width=10, height=8)
