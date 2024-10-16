setwd('/Users/eliselucotte/desktop/server0/Assembly/ortholog_reconstruction/')
library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)
library(stringr)

d=read.table('./2.make_DS_plot/Mv-sup-1065/table_all_DS.txt', header=T, as.is=T, sep='\t')
sexchrom1=c('tig00000068', 'tig00000031')
sexchrom2=c('tig00000080', 'tig00000084','tig00000040')

#d<-mutate(d, line=ifelse(chrom_lag=='Mv-lag-1253-A1_MC03', 282495,ifelse(chrom_lag=='Mv-lag-1253-A1_MC12',343080,0)))
d=d[which(d$chrom1 %in% sexchrom1 & d$chrom2 %in% sexchrom2),]
d=d%>%select(c('ortho_gp','name_lag','chrom_lag','start_lag','end_lag','chrom1','real_id1','chrom2','real_id2','dS','dS_se'))

strat=read.table('./2.make_DS_plot/find_strata/Old_strata/Strata_correspondance_new_annot.txt', header=T, as.is=T, sep='\t')
colnames(strat)=c('gene_old','chrom_old','start_old','end_old','name_lag','chrom_new','start_new','end_new','type_strata')
strat=strat%>% select(c('name_lag','type_strata'))
## CAREFUL HERE IF THE GENE HAS SEVERAL OLD GENE THEN THE LINE WILL BE DUPLICATED !
strat=unique(strat)
#strat[which(strat$name_lag=='Mv-lag-1253-A1_MC03_g2165'),]

d2=left_join(d,strat,c("name_lag"))
#d2[which(d2$real_id1=="Mv-sup-1065-A2_tig00000068_g10173.t2"),]

## RECODE THE GENE NAMES
d2$real_id1=d2$real_id1 %>% str_replace(".t1", "")
d2$real_id1=d2$real_id1 %>% str_replace(".t2", "")
d2$real_id2=d2$real_id2 %>% str_replace(".t1", "")
d2$real_id2=d2$real_id2 %>% str_replace(".t2", "")

#d2[which(d2$real_id1=="Mv-sup-1065-A2_tig00000068_g10173"),]

### ADD THE ORDER OF THE GENES
A1=read.table('/Users/eliselucotte/desktop/server0/Assembly/ortholog_reconstruction/0.data/GTF_final/Mv-sup-1065-A1_final_fixed.bed', header=F, as.is=T, sep='\t')
chrom_A1=c('Mv-sup-1065-A1_tig00000084','Mv-sup-1065-A1_tig00000080','Mv-sup-1065-A1_tig00000040')
A1=A1[which(A1$V1 %in% chrom_A1),]
A1=A1[which(A1$V8=='gene'),]
colnames(A1)=c('contig2','start2','end2','real_id2','V5','V6','V7','V8','V9','V10')
A1=select(A1,c('contig2','start2','end2','real_id2'))
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

##flips the ranks for tig00000080
max_rankPR1=max(d2[which(d2$chrom2=='tig00000084'),]$rankA1)
max_rankPR2=max(d2[which(d2$chrom2=='tig00000080'),]$rankA1)
d2=mutate(d2, rankA1=ifelse(d2$chrom2=='tig00000080',max_rankPR2-rankA1+max_rankPR1,rankA1))

#inversion of the start and end for tig00000084, like in the circos
# end=max(d2[d2$contig2=='Mv-sup-1065-A1_tig00000084',]$end2)
# d2[d2$contig2=='Mv-sup-1065-A1_tig00000084',]$start2=end-d2[d2$contig2=='Mv-sup-1065-A1_tig00000084',]$start2
# d2[d2$contig2=='Mv-sup-1065-A1_tig00000084',]$end2=end-d2[d2$contig2=='Mv-sup-1065-A1_tig00000084',]$end2
# 
# max_rankPR1=max(d2[which(d2$chrom2=='tig00000084'),]$rankA1)
# d2=mutate(d2, rankA1=ifelse(d2$chrom2=='tig00000084',max_rankPR1-rankA1,rankA1))


A2=read.table('/Users/eliselucotte/desktop/server0/Assembly/ortholog_reconstruction/0.data/GTF_final/Mv-sup-1065-A2_final_fixed.bed', header=F, as.is=T, sep='\t')
chrom_A2=c('Mv-sup-1065-A2_tig00000068','Mv-sup-1065-A2_tig00000031')
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
d2=mutate(d2,colors_chromA2=ifelse(d2$chrom1=='tig00000068','purple',ifelse(d2$chrom1=='tig00000031','navyblue','black')))
d2=mutate(d2,colors_chromA1=ifelse(d2$chrom2=='tig00000084','purple',ifelse(d2$chrom2=='tig00000080','plum2',ifelse(d2$chrom2=='tig00000040','navyblue','black'))))


## coloring the genes that are not in the "right" chromosome
#in HD of Mv-lag but PR of Mv-sup
d2[which(d2$chrom_lag=="Mv-lag-1253-A1_MC03A" & d2$chrom1=='tig00000068'),]$colors='green'
#in HD of Mv-lag but PR of Mv-sup, 
d2[which(d2$chrom_lag=="Mv-lag-1253-A1_MC03B" & d2$chrom1=='tig00000031'),]$colors='green'
d2[which(d2$chrom_lag=="Mv-lag-1253-A1_MC12" & d2$chrom1=='tig00000031'),]$colors='green'

## coloring the turquoise strata
new_strata=read.table('./2.make_DS_plot/find_strata/Mv-sup-1065/Turquoise_strata_genes.txt', header=T, as.is=T, sep='\t')
d2[which(d2$name_lag %in% new_strata$name_lag),]$colors='turquoise'
# ggplot(data=d2[which(d2$new_axis>=1.45 & d2$new_axis<=1.611763 & d2$dS<0.05),], aes(x=new_axis, y=dS, color=colors))+geom_point()+ scale_colour_identity()+theme_bw()
# d2[which(d2$chrom_lag=='Mv-lag-1253-A1_MC03B' & d2$new_axis>=1.45 & d2$new_axis<=1.611763 & d2$dS<0.05),]$colors='turquoise'
#write.table(d2[which(d2$real_id1 %in% new_strata$name1),], file=('./2.make_DS_plot/find_strata/new_strata_genes_MvLagA1.txt'),quote=F, row.names = F, sep=('\t'))

## coloring the PARs
PARs=read.table('./2.make_DS_plot/find_strata/Mv-sup-1065/PARS_Mv-Sup-1065.txt', header=T, as.is=T, sep='\t')
d2[which(d2$real_id1 %in% PARs$name1),]$colors='grey'

## coloring the RUBY strata
new_strata2=read.table('./2.make_DS_plot/find_strata/Mv-sup-1065/Ruby_strata_genes.txt', header=T, as.is=T, sep='\t')
d2[which(d2$name_lag %in% new_strata2$name_lag),]$colors='firebrick'

# d2[which(d2$colors=='firebrick'),]
# 
# ggplot(data=d2[which(d2$chrom_lag=='Mv-lag-1253-A1_MC12' & d2$new_axis>=0.035 & d2$new_axis<=0.24),], aes(x=new_axis, y=dS, color=colors))+geom_point()+ scale_colour_identity()+theme_bw()
# d2[which(d2$chrom_lag=='Mv-lag-1253-A1_MC12' & d2$new_axis>=0.035 & d2$new_axis<=0.24),]$colors='firebrick'


write.table(d2, file='./2.make_DS_plot/Mv-sup-1065/database_for_dSplot_Mv-sup-1065.txt', quote = F, row.names = F, sep='\t')

d2=read.table('./2.make_DS_plot/Mv-sup-1065/database_for_dSplot_Mv-sup-1065.txt', header=T, as.is=T)

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

write.table(d2, file='./2.make_DS_plot/Mv-sup-1065/database_for_dSplot_Mv-sup-1065_withrank.txt', quote = F, row.names = F, sep='\t')

#####SEPARATE THE RUBY STRATA IN TWO
d2=read.table('./2.make_DS_plot/Mv-sup-1065/database_for_dSplot_Mv-sup-1065_withrank.txt', header=T, as.is=T)

ruby=d2[which(d2$colors=='firebrick'),]

MC12=ggplot(data=ruby, aes(x=rank_lag,y=dS))+
  geom_point()+
  theme_bw()+
  geom_vline(aes(xintercept= 20), col='red')

MC12_rankA1=ggplot(data=ruby)+
  geom_point(aes(x=rank_lag,y=rankA1/1000))+
  theme_bw()+
  scale_colour_identity()+
  geom_vline(aes(xintercept= 20), col='red')

MC12+MC12_rankA1+plot_layout(ncol = 1)

d2[which(d2$colors=='firebrick' & d2$rank_lag<=20),]$colors='indianred1'
ruby=d2[which(d2$colors=='firebrick' | d2$colors=='indianred1'),]

MC12=ggplot(data=ruby, aes(x=rank_lag,y=dS,color=colors))+
  scale_colour_identity()+
  geom_point()+
  theme_bw()+
  geom_vline(aes(xintercept= 20), col='red')

MC12_rankA1=ggplot(data=ruby)+
  geom_point(aes(x=rank_lag,y=rankA1/1000,color=colors))+
  scale_colour_identity()+
  theme_bw()+
  geom_vline(aes(xintercept= 20), col='red')

MC12+MC12_rankA1+plot_layout(ncol = 1)

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
  "grey"="PARs"
)

#d2=mutate(d2,strata_name=dstrata[d2$colors])
HD_lag=c("Mv-lag-1253-A1_MC03A","Mv-lag-1253-A1_MC03B")
HD_sup=c("Mv-sup-1065-A1_tig00000040")
PR_sup=c("Mv-sup-1065-A1_tig00000084","Mv-sup-1065-A1_tig00000080")
##Recombining HD
HD_RR=d2[which(d2$chrom_lag %in% HD_lag & d2$contig2 %in% HD_sup & d2$colors=="black"),]
d2=mutate(d2,strata_name=ifelse(real_id1 %in% HD_RR$real_id1, "HD_RR",dstrata[d2$colors]))
d2=mutate(d2,colors=ifelse(strata_name=='HD_RR',"grey",colors))
nrow(d2[which(d2$strata_name=='HD_RR'),])
write.table(d2, file='./2.make_DS_plot/Mv-sup-1065/database_for_dSplot_Mv-sup-1065_withrank_new_color.txt', quote = F, row.names = F, sep='\t')

###RENAME THE COLUMNS
d3=d2%>%select(c('ortho_gp','name_lag','chrom_lag','start_lag','end_lag','rank_lag','real_id2','contig2','start2','end2','rankA1','real_id1','contig1','start1','end1','rankA2',"dS",'type_strata',"colors",'strata_name'))
colnames(d3)=c('ortho_gp','gene_lag','chrom_lag','start_lag','end_lag','rank_lag','gene_A1','contig_A1','start_A1','end_A1','rank_A1','gene_A2','contig_A2','start_A2','end_A2','rank_A2','dS','strata_lag','colors','strata_name')
write.table(d3, file='./2.make_DS_plot/Mv-sup-1065/database_for_strata_Mv-sup-1065.txt', quote = F, row.names = F, sep='\t')
