setwd('/Users/eliselucotte/desktop/server1/Assembly/ortholog_reconstruction/Scorzo/dSplot/')
library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)
library(stringr)

d=read.table('table_all_DS.txt', header=T, as.is=T, sep='\t')

d=d%>%select(c('name_lag','chrom_lag','start_lag','end_lag','chrom1','real_id1','chrom2','real_id2','dS','dS_se'))

strat=read.table('../../2.make_DS_plot/find_strata/Old_strata/Strata_correspondance_new_annot.txt', header=T, as.is=T, sep='\t')
colnames(strat)=c('gene_old','chrom_old','start_old','end_old','name_lag','chrom_new','start_new','end_new','type_strata')
strat=strat%>% select(c('name_lag','type_strata'))
## CAREFUL HERE IF THE GENE HAS SEVERAL OLD GENE THEN THE LINE WILL BE DUPLICATED !
strat=unique(strat)
#strat[which(strat$name_lag=='Mv-lag-1253-A1_MC03_g2165'),]

d2=left_join(d,strat,c("name_lag"))

### ADD THE ORDER OF THE GENES
A1=read.table('../0.data/bed/M-scorzo-A1.bed', header=F, as.is=T, sep='\t')
chrom_A1=c("MscorzoA1_tig00000022","MscorzoA1_tig00000006","MscorzoA1_tig00000008","MscorzoA1_tig00000036","MscorzoA1_tig00000123","MscorzoA1_tig00000064","MscorzoA1_tig00000140","MscorzoA1_tig00000163")
A1=A1[which(A1$V1 %in% chrom_A1),]
colnames(A1)=c('contig1','start1','end1','real_id1')
A1$real_id1=A1$real_id1 %>% str_replace(".t1", "")
A1$real_id1=A1$real_id1 %>% str_replace(".t2", "")
new_A1=NULL
for(chrom in chrom_A1)
{
  temp=A1[which(A1$contig1==chrom),]
  temp=temp[order(temp$start1),]
  row.names(temp)<-NULL
  temp=mutate(temp,rankA1=row.names(temp))
  new_A1=rbind.data.frame(new_A1,temp)
}
d2=merge(d2,new_A1,by='real_id1')
d2$rankA1=as.numeric(as.character(d2$rankA1))

##flips the ranks for tig00000080
# max_rankPR1=max(d2[which(d2$chrom2=='tig00000084'),]$rankA1)
# max_rankPR2=max(d2[which(d2$chrom2=='tig00000080'),]$rankA1)
# d2=mutate(d2, rankA1=ifelse(d2$chrom2=='tig00000080',max_rankPR2-rankA1+max_rankPR1,rankA1))

#inversion of the start and end for tig00000084, like in the circos
# end=max(d2[d2$contig2=='Mv-sup-1065-A1_tig00000084',]$end2)
# d2[d2$contig2=='Mv-sup-1065-A1_tig00000084',]$start2=end-d2[d2$contig2=='Mv-sup-1065-A1_tig00000084',]$start2
# d2[d2$contig2=='Mv-sup-1065-A1_tig00000084',]$end2=end-d2[d2$contig2=='Mv-sup-1065-A1_tig00000084',]$end2
# 
# max_rankPR1=max(d2[which(d2$chrom2=='tig00000084'),]$rankA1)
# d2=mutate(d2, rankA1=ifelse(d2$chrom2=='tig00000084',max_rankPR1-rankA1,rankA1))


A2=read.table('../0.data/bed/M-scorzo-A2.bed', header=F, as.is=T, sep='\t')
chrom_A2=c("MscorzoA2_tig00000533","MscorzoA2_tig00000008", "MscorzoA2_tig00000007", "MscorzoA2_tig00000010", "MscorzoA2_tig00000531", "MscorzoA2_tig00000185", "MscorzoA2_tig00000557", "MscorzoA2_tig00005667")
A2=A2[which(A2$V1 %in% chrom_A2),]
colnames(A2)=c('contig2','start2','end2','real_id2')
A2$real_id2=A2$real_id2 %>% str_replace(".t1", "")
A2$real_id2=A2$real_id2 %>% str_replace(".t2", "")
new_A2=NULL
for(chrom in chrom_A2)
{
  temp=A2[which(A2$contig2==chrom),]
  temp=temp[order(temp$start2),]
  row.names(temp)<-NULL
  temp=mutate(temp,rankA2=row.names(temp))
  new_A2=rbind.data.frame(new_A2,temp)
}
d2=merge(d2,new_A2,by='real_id2')
d2$rankA2=as.numeric(as.character(d2$rankA2))

#### DEFINE THE COLOR OF THE PLOTS

dcolors <- c(
  "purpleToCEN"='black',
  "HD_RR"="grey",
  "PR_RR"="grey",
  "blueToCEN"="grey",
  "FCinMAT"="grey",
  "purple"="purple",
  "blue"="blue",
  "orange"="orange",
  "HD_blue"="red",
  "FCnoAnnot"="grey"
)

## coloring the old strata
d2=mutate(d2,colors=dcolors[d2$type_strata])
d2=mutate(d2,chrom_lag=ifelse(chrom_lag=='Mv-lag-1253-A1_MC03' & start_lag<594367,'Mv-lag-1253-A1_MC03A',ifelse(chrom_lag=='Mv-lag-1253-A1_MC03' & start_lag>594367,'Mv-lag-1253-A1_MC03B',chrom_lag)))
d2=mutate(d2,new_axis=d2$start_lag/10^6)

## coloring the chromosomes for the gene order
# d2=mutate(d2,colors_chromA2=ifelse(d2$chrom1=='tig00000068','purple',ifelse(d2$chrom1=='tig00000031','navyblue','black')))
# d2=mutate(d2,colors_chromA1=ifelse(d2$chrom2=='tig00000084','purple',ifelse(d2$chrom2=='tig00000080','plum2',ifelse(d2$chrom2=='tig00000040','navyblue','black'))))


## coloring the genes that are not in the "right" chromosome
#in HD of Mv-lag but PR of Mv-sup
# HD scorzo = MscorzoA1_tig00000123
d2[which(d2$chrom_lag=="Mv-lag-1253-A1_MC03A" & d2$chrom1!='tig00000123'),]$colors='green'
#in PR of Mv-lag but HD of Mv-sup, 
#d2[which(d2$chrom_lag=="Mv-lag-1253-A1_MC03B" & d2$chrom1=='tig00000123'),]$colors='green'
#d2[which(d2$chrom_lag=="Mv-lag-1253-A1_MC12" & d2$chrom1=='tig00000123'),]$colors='green'


# coloring the emerald strata
d2[which(d2$chrom_lag=="Mv-lag-1253-A1_MC03B" & d2$start_lag>1.35*10^6 & d2$start_lag<1.45*10^6 & d2$dS>0),]$colors='springgreen4'

#coloring the quartz rose strata
d2[which(d2$chrom_lag=="Mv-lag-1253-A1_MC12" & d2$start_lag>0.57*10^6 & d2$start_lag<0.75*10^6 & d2$dS>0 & d2$colors=='black'),]$colors='springgreen4'

#coloring the amethyste strata
#d2[which(d2$chrom_lag=="Mv-lag-1253-A1_MC12" & d2$start_lag>0.2*10^6 & d2$start_lag<0.3*10^6 & d2$dS>0 & d2$colors=='grey'),]$colors='plum'


#finding the fragment in PR2
d2[which(d2$chrom1=="tig00000008" & d2$chrom2=="tig00000533" &d2$colors=='grey'),]$colors='mistyrose'



#changing the color of the black stratum
d2[which(d2$chrom_lag=="Mv-lag-1253-A1_MC12" & d2$colors=="black"),]$colors='grey'
write.table(d2, file='database_for_dSplot_M-scorzo.txt', quote = F, row.names = F, sep='\t')

d2=read.table('database_for_dSplot_M-scorzo.txt', header=T, as.is=T)

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

#add the rank of MC12-MC16 and MC03B to paste them together in the plot instead of having 3 distinct chromosomes
rank_MC12=max(d2[which(d2$chrom_lag=="Mv-lag-1253-A1_MC12"),]$rank_lag)
d2[which(d2$chrom_lag=="Mv-lag-1253-A1_MC16"),]$rank_lag=d2[which(d2$chrom_lag=="Mv-lag-1253-A1_MC16"),]$rank_lag+rank_MC12
rank_MC16=max(d2[which(d2$chrom_lag=="Mv-lag-1253-A1_MC16"),]$rank_lag)
d2[which(d2$chrom_lag=="Mv-lag-1253-A1_MC03B"),]$rank_lag=d2[which(d2$chrom_lag=="Mv-lag-1253-A1_MC03B"),]$rank_lag+rank_MC16

write.table(d2, file='database_for_dSplot_M-scorzo_withrank.txt', quote = F, row.names = F, sep='\t')

####GIVE A NAME TO EACH STRATA
dstrata <- c(
  "black"='black',
  "purple"="purple",
  "blue"="blue",
  "orange"="orange",
  "red"="red",
  "indianred1"='light_ruby',
  'firebrick'='dark_ruby',
  'turquoise'='turquoise',
  "green"="in_HD",
  "grey"="PARs",
  "mistyrose"="quartz",
  "springgreen4"="emerald"
  
)

d2=mutate(d2,strata_name=dstrata[d2$colors])
HD_lag=c("Mv-lag-1253-A1_MC03A","Mv-lag-1253-A1_MC03B")
HD_sup=c("Mv-sup-1065-A1_tig00000040")
PR_sup=c("Mv-sup-1065-A1_tig00000084","Mv-sup-1065-A1_tig00000080")
write.table(d2, file='database_for_dSplot_M-scorzo_withrank_new_color.txt', quote = F, row.names = F, sep='\t')


###RENAME THE COLUMNS
d3=d2%>%select(c('name_lag','chrom_lag','start_lag','end_lag','rank_lag','real_id1','contig1','start1','end1','rankA1','real_id2','contig2','start2','end2','rankA2',"dS",'type_strata',"colors",'strata_name'))
colnames(d3)=c('gene_lag','chrom_lag','start_lag','end_lag','rank_lag','gene_A1','contig_A1','start_A1','end_A1','rank_A1','gene_A2','contig_A2','start_A2','end_A2','rank_A2','dS','strata_lag','colors','strata_name')
write.table(d3, file='database_for_strata_M-scorzo.txt', quote = F, row.names = F, sep='\t')
