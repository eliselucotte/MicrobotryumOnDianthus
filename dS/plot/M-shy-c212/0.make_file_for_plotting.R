setwd('/Users/eliselucotte/desktop/server1/Assembly/ortholog_reconstruction/')
library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)
library(stringr)
d=read.table('./2.make_DS_plot/Mv-cat-C212/table_all_DS.txt', header=T, as.is=T, sep='\t')
sexchrom2=c('h2tg000018l', 'h2tg000005l')
sexchrom1=c('h1tg000040l', 'h1tg000006l')

#d<-mutate(d, line=ifelse(chrom_lag=='Mv-lag-1253-A1_MC03', 282495,ifelse(chrom_lag=='Mv-lag-1253-A1_MC12',343080,0)))
d=d[which(d$chrom1 %in% sexchrom1 & d$chrom2 %in% sexchrom2),]
d=d%>%select(c('ortho_gp','name_lag','chrom_lag','start_lag','end_lag','chrom1','real_id1','chrom2','real_id2','dS','dS_se'))

strat=read.table('./2.make_DS_plot/find_strata/Old_strata/Strata_correspondance_new_annot.txt', header=T, as.is=T, sep='\t')
colnames(strat)=c('gene_old','chrom_old','start_old','end_old','name_lag','chrom_new','start_new','end_new','type_strata')
strat=strat%>% select(c('name_lag','type_strata'))
strat=unique(strat)

d2=left_join(d,strat,c("name_lag"))

## RECODE THE GENE NAMES
d2$real_id1=d2$real_id1 %>% str_replace(".t1", "")
d2$real_id1=d2$real_id1 %>% str_replace(".t2", "")
d2$real_id2=d2$real_id2 %>% str_replace(".t1", "")
d2$real_id2=d2$real_id2 %>% str_replace(".t2", "")

### ADD THE ORDER OF THE GENES
A1=read.table('/Users/eliselucotte/desktop/server1/Assembly/ortholog_reconstruction/0.data/GTF_final/Mv-cat-C212-A1_final_fixed.bed', header=F, as.is=T, sep='\t')
chrom_A1=c('Mv-cat-C212-A1_h2tg000018l','Mv-cat-C212-A1_h2tg000005l')
A1=A1[which(A1$V1 %in% chrom_A1),]
A1=A1[which(A1$V8=='gene'),]
colnames(A1)=c('contig2','start2','end2','real_id2','V5','V6','V7','V8','V9','V10')
A1=select(A1,c('contig2','start2','end2','real_id2'))
## remove the amber stratum
A1=A1[which(A1$start2>372152),]
A1[which(A1$contig2=="Mv-cat-C212-A1_h2tg000005l"),]$start2=A1[which(A1$contig2=="Mv-cat-C212-A1_h2tg000005l"),]$start2-372152
A1[which(A1$contig2=="Mv-cat-C212-A1_h2tg000005l"),]$end2=A1[which(A1$contig2=="Mv-cat-C212-A1_h2tg000005l"),]$end2-372152

new_A1=NULL
for(chrom in chrom_A1)
{
  temp=A1[which(A1$contig2==chrom),]
  temp=temp[order(temp$start2),]
  row.names(temp)<-NULL
  temp=mutate(temp,rankA1=row.names(temp))
  new_A1=rbind.data.frame(new_A1,temp)
}
d2=merge(d2,new_A1,by='real_id2')
d2$rankA1=as.numeric(as.character(d2$rankA1))

A2=read.table('/Users/eliselucotte/desktop/server1/Assembly/ortholog_reconstruction/0.data/GTF_final/Mv-cat-C212-A2_final_fixed.bed', header=F, as.is=T, sep='\t')
chrom_A2=c('Mv-cat-C212-A2_h1tg000040l', 'Mv-cat-C212-A2_h1tg000006l')
A2=A2[which(A2$V1 %in% chrom_A2),]
A2=A2[which(A2$V8=='gene'),]
colnames(A2)=c('contig1','start1','end1','real_id1','V5','V6','V7','V8','V9','V10')
A2=select(A2,c('contig1','start1','end1','real_id1'))
new_A2=NULL
for(chrom in chrom_A2)
{
  temp=A2[which(A2$contig1==chrom),]
  temp=temp[order(temp$contig1),]
  row.names(temp)<-NULL
  temp=mutate(temp,rankA2=row.names(temp))
  new_A2=rbind.data.frame(new_A2,temp)
}
d2=merge(d2,new_A2,by='real_id1')
d2$rankA2=as.numeric(as.character(d2$rankA2))

#### DEFINE THE COLOR OF THE PLOTS

dcolors <- c(
  "purpleToCEN"='black',
  "HD_RR"="black",
  "PR_RR"="black",
  "blueToCEN"="black",
  "FCinMAT"="black",
  "purple"="purple",
  "blue"="blue",
  "orange"="orange",
  "HD_blue"="red",
  "FCnoAnnot"="black"
)

## coloring the old strata
d2=mutate(d2,colors=dcolors[d2$type_strata])
d2=mutate(d2,chrom_lag=ifelse(chrom_lag=='Mv-lag-1253-A1_MC03' & start_lag<800038,'Mv-lag-1253-A1_MC03A',ifelse(chrom_lag=='Mv-lag-1253-A1_MC03' & start_lag>800038,'Mv-lag-1253-A1_MC03B',chrom_lag)))
d2=mutate(d2,new_axis=d2$start_lag/10^6)

## coloring the chromosomes for the gene order
d2=mutate(d2,colors_chromA2=ifelse(d2$chrom1=='h1tg000006l','purple',ifelse(d2$chrom1=='h1tg000040l','navyblue','black')))
d2=mutate(d2,colors_chromA1=ifelse(d2$chrom2=='h2tg000005l','purple',ifelse(d2$chrom2=='h2tg000018l','navyblue','black')))


## coloring the genes that are not in the "right" chromosome
d2[which(d2$chrom_lag=="Mv-lag-1253-A1_MC03A" & d2$chrom2=='h2tg000005l'),]$colors='green'
#d2[which(d2$chrom_lag=="Mv-lag-1253-A1_MC03B" & d2$chrom2=='h2tg000018l'),]$colors='green'
#d2[which(d2$chrom_lag=="Mv-lag-1253-A1_MC16" & d2$chrom2=='h2tg0000018l'),]
#d2[which(d2$name_lag=='Mv-lag-1253-A1_MC03_g2112'),]
#d2[which(d2$name_lag=='Mv-lag-1253-A1_MC03_g2111'),]


## coloring the PARs
PARs=read.table('./2.make_DS_plot/find_strata/Mv-cat-C212/PARS_Mv-Cat-C212.txt', header=T, as.is=T, sep='\t')
d2[which(d2$real_id1 %in% PARs$name1),]$colors='grey'



## adding the dark ruby strata 
d2[which(d2$start_lag>200000 & d2$start_lag<286224 & d2$chrom_lag=="Mv-lag-1253-A1_MC12"),]$colors='firebrick'

write.table(d2, file='./2.make_DS_plot/Mv-cat-C212/database_for_dSplot_Mv-cat-C212.txt', quote = F, row.names = F, sep='\t')


###ADD THE RANK OF THE GENES FOR LAGERHEIMII
new_d=NULL
for(chrom in unique(d2$chrom_lag))
{
  temp=d2[which(d2$chrom_lag==chrom),]
  temp=temp[order(temp$start_lag),]
  row.names(temp)<-NULL
  temp=mutate(temp,rank_lag=row.names(temp))
  new_d=rbind.data.frame(new_d,temp)
}

new_d$rank_lag=as.numeric(as.character(new_d$rank_lag))
d2=new_d

rank_MC12=max(d2[which(d2$chrom_lag=="Mv-lag-1253-A1_MC12"),]$rank_lag)
d2[which(d2$chrom_lag=="Mv-lag-1253-A1_MC16"),]$rank_lag=d2[which(d2$chrom_lag=="Mv-lag-1253-A1_MC16"),]$rank_lag+rank_MC12
rank_MC16=max(d2[which(d2$chrom_lag=="Mv-lag-1253-A1_MC16"),]$rank_lag)
d2[which(d2$chrom_lag=="Mv-lag-1253-A1_MC03B"),]$rank_lag=d2[which(d2$chrom_lag=="Mv-lag-1253-A1_MC03B"),]$rank_lag+rank_MC16

write.table(d2, file='./2.make_DS_plot/Mv-cat-C212/database_for_dSplot_Mv-cat-C212_withrank.txt', quote = F, row.names = F, sep='\t')

####GIVE A NAME TO EACH STRATA
dstrata <- c(
  "black"='black',
  "purple"="purple",
  "blue"="blue",
  "orange"="orange",
  "red"="red",
  'firebrick'='dark_ruby',
  "green"="in_HD",
  "grey"="PARs"
)

d2=mutate(d2,strata_name=dstrata[d2$colors])
write.table(d2, file='./2.make_DS_plot/Mv-cat-C212/database_for_dSplot_Mv-cat-C212_withrank_new_color.txt', quote = F, row.names = F, sep='\t')

###RENAME THE COLUMNS
d3=d2%>%select(c('ortho_gp','name_lag','chrom_lag','start_lag','end_lag','rank_lag','real_id2','contig2','start2','end2','rankA1','real_id1','contig1','start1','end1','rankA2',"dS",'type_strata',"colors",'strata_name'))
colnames(d3)=c('ortho_gp','gene_lag','chrom_lag','start_lag','end_lag','rank_lag','gene_A1','contig_A1','start_A1','end_A1','rank_A1','gene_A2','contig_A2','start_A2','end_A2','rank_A2','dS','strata_lag','colors','strata_name')
write.table(d3, file='./2.make_DS_plot/Mv-cat-C212/database_for_strata_Mv-cat-C212.txt', quote = F, row.names = F, sep='\t')

