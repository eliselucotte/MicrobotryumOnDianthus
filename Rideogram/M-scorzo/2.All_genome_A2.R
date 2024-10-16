#install.packages("RIdeogram")
require(RIdeogram)
library(dplyr)
library(tidyr)
library(stringr)
setwd("/Users/eliselucotte/desktop/server1/Assembly/ortholog_reconstruction/Scorzo/Rideogram/")
specie1='Mv-lag-1253-A2'
specie2='M-scorzo-A2'

# PR a2: MscorzoA2_tig00000008_MP006917
# hd1 a2: MscorzoA2_tig00000185_g4967.t1
# hd2 a2: MscorzoA2_tig00000185_g4968.t1
# PR a1: MscorzoA1_tig00000006_g593.t1
# hd1 a1: MscorzoA1_tig00000123_g5608.t1
# hd2 a1: MscorzoA1_tig00000123_g5609.t1

PR_chroms=c('M-scor-PR','M-scor-PR2','M-scor-PR3','M-scor-PR4','M-scor-PR5','M-scor-PR6','PR1','PR2')
HD_chroms=c('M-scor-HD','HD')
dchrom <- c(
  "tig00000185"="M-scor-HD",
  "tig00000008"="M-scor-PR",
  "tig00000533"="M-scor-PR2",
  "tig00000531"="M-scor-PR3",
  "tig00000007"='M-scor-PR4',
  "tig00005667"="M-scor-PR5",
  "tig00000010"="M-scor-PR6",
  "MC03"="HD",
  "MC12"="PR1",
  "MC16"="PR2"
)

dspecies <-c(
  "Mv-lag-1253-A2"="M-lag-A2",
  "M-scorzo-A2"="M-scorzo-A2"
)

#define the HD and PR genes' names
HD_genes= c('Mv-lag-1253-A1_MC03_g1964','Mv-lag-1253-A2_MC03_g2006')
PR_gene='Mv-lag-1253-A2_MC12_PR'
#which contig to change direction on the plot
#to_inv=c("tig00000006")
to_inv=c()

d=read.table('synteny_Mv-lag-1253-A2_M-scorzo-A2.txt', header=T, as.is=T, sep='\t')
d_PR_genes=read.table('PR_genes.txt', header=T, as.is=T, sep='\t')
A1=rbind(d,d_PR_genes[2,])
A1=d

contig1=read.table(paste(specie1,'_contig.bed',sep=''), header=T, as.is=T, sep='\t')
contig1=contig1[which(as.numeric(contig1$End) > 100000),]

contig2=read.table(paste(specie2,'_contig.bed',sep=''), header=T, as.is=T, sep='\t')
contig2=contig2[which(as.numeric(contig2$End) > 100000),]

A1=A1[which(A1$chrom_lag%in%contig1$Chr & A1$chrom1%in%contig2$Chr),]

contig1<-contig1 %>%
  arrange(Chr)  

contig2<-contig2 %>%
  arrange(Chr)

contigs=bind_rows(contig1,contig2)


for (contig in to_inv)
{
  end=contigs[which(contigs$Chr==contig),]$End
  #change the coordinates for the synteny databases
  A1[which(A1$chrom1==contig),]$start1=end-A1[which(A1$chrom1==contig),]$start1
  A1[which(A1$chrom1==contig),]$end1=end-A1[which(A1$chrom1==contig),]$end1
  
}

##Change the names of the contigs for plotting, according to the dictionnary at the begining
A1=A1 %>% mutate(chrom_lag=ifelse(chrom_lag %in% names(dchrom),dchrom[A1$chrom_lag],chrom_lag))
A1=A1 %>% mutate(chrom1=ifelse(chrom1 %in% names(dchrom),dchrom[A1$chrom1],chrom1))


contigs=contigs %>% mutate(Chr=ifelse(Chr %in% names(dchrom),dchrom[contigs$Chr],Chr))
contig1=contig1 %>% mutate(Chr=ifelse(Chr %in% names(dchrom),dchrom[contig1$Chr],Chr))
contig2=contig2 %>% mutate(Chr=ifelse(Chr %in% names(dchrom),dchrom[contig2$Chr],Chr))


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
ideogram(karyotype = contigs, synteny = t_plot, output=paste('All_contigs_',specie1,specie2,'.svg', sep=''))
convertSVG(paste('All_contigs_',specie1,specie2,'.svg', sep=''), device = "png")




table(t_A1[which(t_A1$chrom_lag=='HD'),]$chrom1)
table(t_A1[which(t_A1$chrom_lag=='PR1'),]$chrom1)
table(t_A1[which(t_A1$chrom_lag=='PR2'),]$chrom1)


