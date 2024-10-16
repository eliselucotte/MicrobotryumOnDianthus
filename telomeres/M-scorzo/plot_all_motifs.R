library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
MAT='A1'

size=18
path="/Users/eliselucotte/desktop/server1/Assembly/ortholog_reconstruction/Scorzo/telomeres/"
species='M-scorzo'
dchromA2 <- c(
    "tig00000185"="M-scor-A2-HD",
    "tig00000008"="M-scor-A2-PR2",
    "tig00000533"="M-scor-A2-PR1b",
    "tig00000531"="M-scor-A2-PR3",
    "tig00005667"="M-scor-A2-PR4",
    "tig00005663"="M-scor-A2-PR1a"
  )
dchromA1 <- c(
    "tig00000123"="M-scor-A1-HD",
    "tig00000006"="M-scor-A1-PR1",
    "tig00000008"="M-scor-A1-PR2",
    "tig00000064"="M-scor-A1-PR3",
    "tig00000163"="M-scor-A1-PR4"
  )
window=1000


t=read.table(paste(path,species,'-',MAT,'-count_telomeres_motifs_',window,'bp.txt',sep=''), as.is=T, header=T )
t$interval=t$interval/window
t=mutate(t,contig2=dchromA1[t$contig])
#positions of the telomeres
#View(t[which(t$contig=='h2tg000005l'),])
# 1327*window

A1=ggplot(data=t, aes(y=count, x=interval))+
  geom_line()+
  theme_bw(base_size = size)+
  ylim(c(0,17))+
  labs(title='A1 TTAGGG',x='position along the contig', y='motif count')+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  facet_wrap(~contig2, scales='free_x')

t=read.table(paste(path,species,'-',MAT,'-count_telomeres_motifs_inv_compl_',window,'bp.txt',sep=''), as.is=T, header=T )
t$interval=t$interval/window
t=mutate(t,contig2=dchromA1[t$contig])
#positions of the telomeres
#View(t[which(t$contig=='h2tg000005l'),])
# 37*window

A1_inv_compl=ggplot(data=t, aes(y=count, x=interval))+
  geom_line()+
  theme_bw(base_size = size)+
  ylim(c(0,17))+
  labs(title='A1 CCCTAA',x='position along the contig', y='motif count')+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  facet_wrap(~contig2, scales='free_x')

t=read.table(paste(path,species,'-',MAT,'-count_telomeres_motifs_inv_',window,'bp.txt',sep=''), as.is=T, header=T )
t$interval=t$interval/window
t=mutate(t,contig2=dchromA1[t$contig])
A1_inv=ggplot(data=t, aes(y=count, x=interval))+
  geom_line()+
  theme_bw(base_size = size)+
  ylim(c(0,17))+
  labs(title='A1 GGGATT',x='position along the contig', y='motif count')+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  facet_wrap(~contig2, scales='free_x')

t=read.table(paste(path,species,'-',MAT,'-count_telomeres_motifs_compl_',window,'bp.txt',sep=''), as.is=T, header=T )
t$interval=t$interval/window
t=mutate(t,contig2=dchromA1[t$contig])
A1_compl=ggplot(data=t, aes(y=count, x=interval))+
  geom_line()+
  theme_bw(base_size = size)+
  ylim(c(0,17))+
  labs(title='A1 AATCCC',x='position along the contig', y='motif count')+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  facet_wrap(~contig2, scales='free_x')


p=A1+A1_inv_compl+A1_compl+A1_inv+plot_layout(ncol = 2, nrow=2)

ggsave(paste(path,species,'-count_telomeres_motifs_all_A1.png',sep=''),plot=p,width=20, height=16)



MAT='A2'
t=read.table(paste(path,species,'-',MAT,'-count_telomeres_motifs_',window,'bp.txt',sep=''), as.is=T, header=T )
t$interval=t$interval/window
t=mutate(t,contig2=dchromA2[t$contig])
#positions of the telomeres
#View(t[which(t$contig=='h2tg000005l'),])
# 1327*window

A2=ggplot(data=t, aes(y=count, x=interval))+
  geom_line()+
  theme_bw(base_size = size)+
  ylim(c(0,17))+
  labs(title='A2 TTAGGG',x='position along the contig', y='motif count')+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  facet_wrap(~contig2, scales='free_x')

t=read.table(paste(path,species,'-',MAT,'-count_telomeres_motifs_inv_compl_',window,'bp.txt',sep=''), as.is=T, header=T )
t$interval=t$interval/window
t=mutate(t,contig2=dchromA2[t$contig])
#positions of the telomeres
#View(t[which(t$contig=='h2tg000005l'),])
# 37*window

A2_inv_compl=ggplot(data=t, aes(y=count, x=interval))+
  geom_line()+
  theme_bw(base_size = size)+
  ylim(c(0,17))+
  labs(title='A2 CCCTAA',x='position along the contig', y='motif count')+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  facet_wrap(~contig2, scales='free_x')

t=read.table(paste(path,species,'-',MAT,'-count_telomeres_motifs_inv_',window,'bp.txt',sep=''), as.is=T, header=T )
t$interval=t$interval/window
t=mutate(t,contig2=dchromA2[t$contig])
A2_inv=ggplot(data=t, aes(y=count, x=interval))+
  geom_line()+
  theme_bw(base_size = size)+
  ylim(c(0,17))+
  labs(title='A2 GGGATT',x='position along the contig', y='motif count')+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  facet_wrap(~contig2, scales='free_x')

t=read.table(paste(path,species,'-',MAT,'-count_telomeres_motifs_compl_',window,'bp.txt',sep=''), as.is=T, header=T )
t$interval=t$interval/window
t=mutate(t,contig2=dchromA2[t$contig])
A2_compl=ggplot(data=t, aes(y=count, x=interval))+
  geom_line()+
  theme_bw(base_size = size)+
  ylim(c(0,17))+
  labs(title='A2 AATCCC',x='position along the contig', y='motif count')+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  facet_wrap(~contig2, scales='free_x')


p=A2+A2_inv_compl+A2_compl+A2_inv+plot_layout(ncol = 2, nrow=2)

ggsave(paste(path,species,'-count_telomeres_motifs_all_A2.png',sep=''),plot=p,width=20, height=16)



# 
# p=A1+A2+A1_inv_compl+A2_inv_compl+plot_layout(ncol = 2, nrow=2,width=c(0.6,0.4))
# p
# ggsave(paste(path,species,'-count_telomeres_motifs_all.png',sep=''),plot=p,width=20, height=16)
