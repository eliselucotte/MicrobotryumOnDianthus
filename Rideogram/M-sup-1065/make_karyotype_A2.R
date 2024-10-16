#install.packages("RIdeogram")
require(RIdeogram)
library(dplyr)
library(tidyr)
library(stringr)
setwd("/Users/eliselucotte/desktop/server0/Assembly/ortholog_reconstruction/")
specie1='Mv-lag-1253-A2'
specie2='Mv-sup-1065-A2'
chroms_keep_species1=c('MC12','MC16','MC03')
chroms_keep_species2=c('tig00000031','tig00000068')

PR_chroms=c('M-sup-PR','PR1','PR2')
HD_chroms=c('M-sup-HD','HD')
dchrom <- c(
  "tig00000031"="M-sup-HD",
  "tig00000068"="M-sup-PR",
  "MC03"="HD",
  "MC12"="PR1",
  "MC16"="PR2"
)

dspecies <-c(
  "Mv-lag-1253-A2"="M-lag-A2",
  "Mv-sup-1065-A2"="M-sup-A2"
)

#define the HD and PR genes' names
HD_genes= c('Mv-lag-1253-A2_MC03_g2006')
PR_gene='Mv-lag-1253-A2_MC12_PR'
#which contig to change direction on the plot
to_inv=c("tig00000031")
#to_inv=c()

d_PR_genes=read.table('5.Rideogram/Mv-sup/PR_genes.txt', header=T, as.is=T, sep='\t')
d_centro=read.table('5.Rideogram/Mv-sup/bedfiles/Positions_centromeres_A2.txt', header=T, as.is=T, sep='\t')
d_legend=read.table('5.Rideogram/Mv-sup/bedfiles/legend_marker_A2.txt', header=T, as.is=T, sep='\t')

contig1=read.table(paste('5.Rideogram/Mv-sup/',specie1,'_contig.bed',sep=''), header=T, as.is=T, sep='\t')
contig1=contig1[which(contig1$Chr%in%chroms_keep_species1),]
contig2=read.table(paste('5.Rideogram/Mv-sup/',specie2,'_contig.bed',sep=''), header=T, as.is=T, sep='\t')
contig2=contig2[which(contig2$Chr%in%chroms_keep_species2),]

contig1<-contig1 %>%
  arrange(Chr)  

contig2<-contig2 %>%
  #mutate(Chr =  factor(Chr, levels = c('tig00000040','tig00000084','tig00000080'))) %>%
  arrange(Chr)

contigs=bind_rows(contig1,contig2)


#----  RIdeogram plot genomes of the same size
#---- this is a problem if genome size if different. 
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


contigs=contigs%>% mutate(fill=ifelse(Chr %in% PR_chroms,'800080',
                                      ifelse(Chr %in% HD_chroms,'003399','0000FF00')))





contigs$Chr=as.character(contigs$Chr)
#change species name
#ideogram(karyotype = contigs, label=d_legend, label_type = "marker", output=paste('5.Rideogram/Mv-sup/',specie1,specie2,'karyotype.svg', sep='')) #overlaid=d_centro
ideogram(karyotype = contigs, overlaid=d_centro, label=d_legend, label_type = "marker", output=paste('5.Rideogram/Mv-sup/',specie1,specie2,'karyotype_A2.svg', sep=''))
#ideogram(karyotype = contigs, overlaid=d_centro, output=paste('5.Rideogram/Mv-sup/',specie1,specie2,'karyotype.svg', sep=''))

convertSVG(paste('5.Rideogram/Mv-sup/',specie1,specie2,'karyotype_A2.svg', sep=''), device = "png")
