#install.packages("RIdeogram")
require(RIdeogram)
library(dplyr)
library(tidyr)
library(stringr)
setwd("/Users/eliselucotte/desktop/server1/Assembly/ortholog_reconstruction/Scorzo/Rideogram/")

specie1='M-scorzo-A1'
specie2='M-scorzo-A2'
chroms_keep_species1=c("tig00000123","tig00000006","tig00000064","tig00000008","tig00000163")
chroms_keep_species2=c("tig00000185","tig00000008","tig00000533","tig00000531","tig00005667","tig00005663")

PR_chroms=c('M-scor-A2-PR1a','M-scor-A2-PR1b','M-scor-A2-PR2','M-scor-A2-PR3','M-scor-A2-PR4','M-scor-A1-PR1','M-scor-A1-PR2','M-scor-A1-PR3','M-scor-A1-PR4')
HD_chroms=c('M-scor-A1-HD','M-scor-A2-HD')
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
dspecies <-c(
  "M-scorzo-A1"="M-scorzo-A1",
  "M-scorzo-A2"="M-scorzo-A2"
)

#define the HD and PR genes' names
HD_genes= c('MscorzoA1_tig00000123_g5608')
PR_gene='MscorzoA1_tig00000006_g593'
#which contig to change direction on the plot
#to_inv=c("tig00000185","tig00000533","tig00000531","tig00005667")
to_inv_A2=c("M-scor-A2-HD","M-scor-A2-PR2","M-scor-A2-PR1b","M-scor-A2-PR3","M-scor-A2-PR4","M-scor-A2-PR1a")
to_inv_A1=c("M-scor-A1-PR2")

d=read.table('synteny_M-scorzo-A1_M-scorzo-A2.txt', header=T, as.is=T, sep='\t')
# d_PR_genes=read.table('PR_genes.txt', header=T, as.is=T, sep='\t')
# d=rbind(d,d_PR_genes[2,])

A1=d[which(d$chrom1 %in% chroms_keep_species2),]
A1=A1[which(A1$chrom_lag %in% chroms_keep_species1),]

contig1=read.table(paste(specie1,'_contig.bed',sep=''), header=T, as.is=T, sep='\t')
contig1=contig1[which(contig1$Chr %in% unique(A1$chrom_lag)),]

contig2=read.table(paste(specie2,'_contig.bed',sep=''), header=T, as.is=T, sep='\t')
contig2=contig2[which(contig2$Chr %in% unique(A1$chrom1)),]

contig1<-contig1 %>%
  mutate(Chr =  factor(Chr, levels = c("tig00000123","tig00000022","tig00000008","tig00000006","tig00000064","tig00000036","tig00000163"))) %>%
  arrange(Chr)

contig2<-contig2 %>%
  mutate(Chr =  factor(Chr, levels = c("tig00000185","tig00000008","tig00005663","tig00000533","tig00000531","tig00000010","tig00005667"))) %>%
  arrange(Chr)

contigs1=contig1 %>% mutate(Chr=dchromA1[as.character(Chr)])
contigs2=contig2 %>% mutate(Chr=dchromA2[as.character(Chr)])
contigs=bind_rows(contigs1,contigs2)
A1=A1 %>% mutate(chrom1=dchromA2[A1$chrom1])
A1=A1 %>% mutate(chrom_lag=dchromA1[A1$chrom_lag])

for (contig in to_inv_A2)
{
  end=contigs2[which(contigs2$Chr==contig),]$End
  #change the coordinates for the synteny databases
  A1[which(A1$chrom1==contig),]$start1=end-A1[which(A1$chrom1==contig),]$start1
  A1[which(A1$chrom1==contig),]$end1=end-A1[which(A1$chrom1==contig),]$end1
}
for (contig1 in to_inv_A1)
{
  end=contigs1[which(contigs1$Chr==contig1),]$End
  #change the coordinates for the synteny databases
  A1[which(A1$chrom_lag==contig1),]$start_lag=end-A1[which(A1$chrom_lag==contig1),]$start_lag
  A1[which(A1$chrom_lag==contig1),]$end_lag=end-A1[which(A1$chrom_lag==contig1),]$end_lag
}
##Change the names of the contigs for plotting, according to the dictionnary at the begining




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

if(small$species==contigs1$species[1]) {
  contigs <- rbind(contigs1, small, contigs2)
} else { contigs <- rbind(contigs1, contigs2, small) }


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
  species1=c(species1,which(contigs1$Chr==t_A1$chrom_lag[i]))
  species2=c(species2,which(contigs2$Chr==t_A1$chrom1[i]))
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
ideogram(karyotype = contigs, synteny = t_plot, output=paste(specie1,specie2,'.svg', sep=''))
convertSVG(paste(specie1,specie2,'.svg', sep=''), device = "png")
