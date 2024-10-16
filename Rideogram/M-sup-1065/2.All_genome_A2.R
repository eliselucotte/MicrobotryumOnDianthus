#install.packages("RIdeogram")
require(RIdeogram)
library(dplyr)
library(tidyr)
library(stringr)
setwd("/Users/eliselucotte/desktop/server1/Assembly/ortholog_reconstruction/")
specie1='Mv-lag-1253-A2'
specie2='Mv-sup-1065-A2'

PR_chroms=c('MC12','MC16','tig00000068')
HD_chroms=c('MC03','tig00000031')
# dchrom <- c(
#   "tig00000031"="M-sup-HD",
#   "tig00000068"="M-sup-PR",
#   "MC03"="HD",
#   "MC12"="PR1",
#   "MC16"="PR2"
# )

dspecies <-c(
  "Mv-lag-1253-A2"="M-lag-A2",
  "Mv-sup-1065-A2"="M-sup-A2"
)

#define the HD and PR genes' names
HD_genes= c('Mv-lag-1253-A2_MC03_g2006')
PR_gene='Mv-lag-1253-A2_MC12_PR'
#which contig to change direction on the plot
to_inv=c("tig00000031","MC13","MC14-1","MC14-2","MC14-3","MC13","MC09","MC05-1","MC05-2")
#to_inv=c()

d=read.table('3.circos/synteny_Mv-lag-1253-A2_Mv-sup-1065-A2.txt', header=T, as.is=T, sep='\t')
d_PR_genes=read.table('5.Rideogram/Mv-sup/PR_genes.txt', header=T, as.is=T, sep='\t')
A1=rbind(d,d_PR_genes[2,])


contig1=read.table(paste('5.Rideogram/Mv-sup/',specie1,'_contig.bed',sep=''), header=T, as.is=T, sep='\t')
contig1=contig1[which(as.numeric(contig1$End) > 100000),]

contig2=read.table(paste('5.Rideogram/Mv-sup/',specie2,'_contig.bed',sep=''), header=T, as.is=T, sep='\t')
contig2=contig2[which(as.numeric(contig2$End) > 100000),]

A1=A1[which(A1$chrom_lag%in%contig1$Chr & A1$chrom1%in%contig2$Chr),]



contigs=bind_rows(contig1,contig2)


for (contig in to_inv)
{
  end=contigs[which(contigs$Chr==contig),]$End
  #change the coordinates for the synteny databases
  A1[which(A1$chrom1==contig),]$start1=end-A1[which(A1$chrom1==contig),]$start1
  A1[which(A1$chrom1==contig),]$end1=end-A1[which(A1$chrom1==contig),]$end1
  
}

# ##Change the names of the contigs for plotting, according to the dictionnary at the begining
# A1=A1 %>% mutate(chrom_lag=ifelse(chrom_lag %in% names(dchrom),dchrom[A1$chrom_lag],chrom_lag))
# A1=A1 %>% mutate(chrom1=ifelse(chrom1 %in% names(dchrom),dchrom[A1$chrom1],chrom1))
# 
# 
# contigs=contigs %>% mutate(Chr=ifelse(Chr %in% names(dchrom),dchrom[contigs$Chr],Chr))
# contig1=contig1 %>% mutate(Chr=ifelse(Chr %in% names(dchrom),dchrom[contig1$Chr],Chr))
# contig2=contig2 %>% mutate(Chr=ifelse(Chr %in% names(dchrom),dchrom[contig2$Chr],Chr))

contig1<-contig1 %>%
  mutate(Chr =  factor(Chr, levels = c('MC03','MC12','MC16',"MC01","MC02-1","MC02-2","MC11-1","MC11-2","MC09","MC05-1","MC05-2","MC08","MC04","MC07","MC06","MC14-3","MC14-2","MC14-1","MC13","MC10","MC15-1","MC15-2","C23"))) %>%
  arrange(Chr)  

contig2<-contig2 %>%
  mutate(Chr =  factor(Chr, levels = c("tig00000031","tig00000068","tig00000089","tig00000012","tig00000013","tig00000014","tig00000016","tig00000017","tig00000020","tig00000022","tig00000023","tig00000030","tig00000033","tig00000065","tig00000071","tig00000092","tig00000093"))) %>%
  arrange(Chr)
contigs=bind_rows(contig1,contig2)

#----  RIdeogram plot genomes of the same size
#---- this is a probelm if genome size if different. 
#----  so we create a false chromosome that a the size of the difference between the 2 genomes

gapsize <- contigs %>% group_by(species) %>% summarise(size = sum(End) ) %>% summarise(diff = max(size) - min(size))
sp <- contigs %>% group_by(species) %>% summarise(size = sum(End) ) %>%  filter(size == min(size)) %>% select(species)

small  <- data.frame(Chr = "none",  
                     Start = 0, 
                     End = gapsize$diff, 
                     fill = "#0000FF00",  #"#FF0000CC", # "fffffff", #empty fill
                     species = sp$species, #name of 
                     size = 12, 
                     color = 25252525) 

if(small$species==contig1$species[1]) {
  contigs <- rbind(contig1, small, contig2)
} else { contigs <- rbind(contig1, contig2, small) }


A1=A1 %>% mutate(fill=ifelse(name_lag %in% HD_genes,'ff0000',
                             ifelse(name_lag %in% PR_gene,'ff0000',
                                    ifelse(chrom_lag %in% PR_chroms,'800080',
                                           ifelse(chrom_lag %in% HD_chroms,'003399','CCCCCC')))))

t_A1=A1 %>% select(c(chrom_lag,start_lag,end_lag,chrom1,start1,end1,fill))

contigs=contigs%>% mutate(fill=ifelse(Chr %in% PR_chroms,'800080',
                                      ifelse(Chr %in% HD_chroms,'003399','CCCCCC')))


#get the ranks of the contig on the contig dataframe (I guess it is the ranks)
species1=NULL
species2=NULL
for(i in seq(1,nrow(t_A1)))
{
  species1=c(species1,which(contig1$Chr==t_A1$chrom_lag[i]))
  species2=c(species2,which(contig2$Chr==t_A1$chrom1[i]))
}
t_A1=cbind(t_A1,species1,species2)
contigs$Chr=as.character(contigs$Chr)
#change species name
contigs=contigs %>% mutate(species=dspecies[contigs$species])

t_plot=select(t_A1, c(species1,start_lag,end_lag,species2,start1,end1,fill))
colnames(t_plot)=c("Species_1","Start_1","End_1","Species_2","Start_2","End_2","fill")
t_plot=t_plot[order(t_plot$fill,decreasing=T),]
end=t_plot[which(t_plot$fill=='ff0000'),]
rest=t_plot[-which(t_plot$fill=='ff0000'),]
t_plot=rbind(rest,end)
ideogram(karyotype = contigs, synteny = t_plot, output=paste('5.Rideogram/Mv-sup/All_contigs_',specie1,specie2,'.svg', sep=''))
convertSVG(paste('5.Rideogram/Mv-sup/All_contigs_',specie1,specie2,'.svg', sep=''), device = "png")
