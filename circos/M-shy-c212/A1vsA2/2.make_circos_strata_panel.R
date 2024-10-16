###documentation links:
##https://www.royfrancis.com/beautiful-circos-plots-in-r/#data
##https://jokergoo.github.io/circlize_book/book/genomic-plotting-region.html#genomic-rectangles
##centromere file : @GEEserver03 /home/marine/MS_newSpecies/allGenomes_centromereCoordinates_newSpecies.tab
setwd("/Users/eliselucotte/desktop/server0/Assembly/ortholog_reconstruction/3.circos/Mv-Cat/A1vsA2/")
library(circlize)
library(dplyr)
library(tidyr)
library(stringr)

species1=paste("Mv-cat-C212-A1")
species2=paste("Mv-cat-C212-A2")

chroms_lag=c('h2tg000005l','h2tg000018l')
chroms2=c('h1tg000040l','h1tg000006l')

#other strata
strata=read.table("/Users/eliselucotte/desktop/server0/Assembly/ortholog_reconstruction/2.make_DS_plot/Mv-cat-C212/database_for_strata_Mv-cat-C212.txt", header=T, as.is=T, sep='\t')
strata1=select(strata, c("contig_A1","start_A1","end_A1","colors"))
colnames(strata1)=c('contig','start','end','colors')
strata1$contig=strata1$contig %>% str_replace("Mv-cat-C212-A1_", "")

strata2=select(strata, c("contig_A2","start_A2","end_A2","colors"))
colnames(strata2)=c('contig','start','end','colors')
strata2$contig=strata2$contig %>% str_replace("Mv-cat-C212-A2_", "")
strata=rbind(strata1, strata2)

##Import the synteny data
d=read.table(paste('synteny_',species1,'_',species2,'-noamber.txt',sep=''), header=T, as.is=T, sep='\t')
d=d[which(d$chrom_lag %in% chroms_lag & d$chrom1 %in% chroms2),]

#Import the contig informations
contig1=read.table(paste(species1,'_contig.bed',sep=''), header=T, as.is=T, sep='\t')
contig1=contig1[which(contig1$Chr %in% unique(d$chrom_lag)),]
contig2=read.table(paste(species2,'_contig.bed',sep=''), header=T, as.is=T, sep='\t')
contig2=contig2[which(contig2$Chr %in% unique(d$chrom1)),]

##Reorder the contig in the order I want them to appear in the plot
contig1=contig1[c(2,1),]
#contig2=contig2[c(2,1),]

#bind both tables
contigs=rbind.data.frame(contig1,contig2)

#Import the gene positions + the centromere positions
data_genes=read.table(paste('Positions_genes.txt',sep=''), header=T, as.is=T, sep='\t')
data_centromere_A1=read.table(paste('../Positions_centromeres_A1.txt',sep=''), header=T, as.is=T, sep='\t')
data_centromere_A2=read.table(paste('../Positions_centromeres_A2.txt',sep=''), header=T, as.is=T, sep='\t')
data_centromeres=rbind(data_centromere_A1[which(data_centromere_A1$contig %in% chroms_lag),],data_centromere_A2[which(data_centromere_A2$contig %in% chroms2),])

data_telomere_A1=read.table(paste('../Positions_telomeres_A1.txt',sep=''), header=T, as.is=T, sep='\t')
data_telomere_A2=read.table(paste('../Positions_telomeres_A2.txt',sep=''), header=T, as.is=T, sep='\t')
data_telomeres=rbind(data_telomere_A1[which(data_telomere_A1$contig %in% chroms_lag),],data_telomere_A2[which(data_telomere_A2$contig %in% chroms2),])

data_pheromones_A1=read.table(paste('../Positions_pheromones_A1.txt',sep=''), header=T, as.is=T, sep='\t')
data_pheromones_A2=read.table(paste('../Positions_pheromones_A2.txt',sep=''), header=T, as.is=T, sep='\t')
data_pheromones=rbind(data_pheromones_A1[which(data_pheromones_A1$contig %in% chroms_lag),],data_pheromones_A2[which(data_pheromones_A2$contig %in% chroms2),])

#Make the database for the genomic links
nuc1=select(d, c('chrom_lag', 'start_lag','end_lag'))
nuc2=select(d, c('chrom1', 'start1','end1'))

##Which contig to you want to invert the orientation for?
to_inv=c("h1tg000006l")
for (contig in to_inv)
{
  end=contigs[which(contigs$Chr==contig),]$End
  #change the coordinates for the genes
  data_genes[which(data_genes$contig==contig),]$start=end-data_genes[which(data_genes$contig==contig),]$start
  data_genes[which(data_genes$contig==contig),]$end=end-data_genes[which(data_genes$contig==contig),]$end
  #change the coordinates for the centromeres
  data_centromeres[which(data_centromeres$contig==contig),]$start=end-data_centromeres[which(data_centromeres$contig==contig),]$start
  data_centromeres[which(data_centromeres$contig==contig),]$end=end-data_centromeres[which(data_centromeres$contig==contig),]$end
  #change the coordinates for the telomeres
  data_telomeres[which(data_telomeres$contig==contig),]$start=end-data_telomeres[which(data_telomeres$contig==contig),]$start
  data_telomeres[which(data_telomeres$contig==contig),]$end=end-data_telomeres[which(data_telomeres$contig==contig),]$end
  #change the coordinates for the pheromones
  data_pheromones[which(data_pheromones$contig==contig),]$start=end-data_pheromones[which(data_pheromones$contig==contig),]$start
  data_pheromones[which(data_pheromones$contig==contig),]$end=end-data_pheromones[which(data_pheromones$contig==contig),]$end
  # #change the coordinates for the synteny databases
  nuc1[which(nuc1$chrom_lag==contig),]$start_lag=end-nuc1[which(nuc1$chrom_lag==contig),]$start_lag
  nuc1[which(nuc1$chrom_lag==contig),]$end_lag=end-nuc1[which(nuc1$chrom_lag==contig),]$end_lag
  nuc2[which(nuc2$chrom1==contig),]$start1=end-nuc2[which(nuc2$chrom1==contig),]$start1
  nuc2[which(nuc2$chrom1==contig),]$end1=end-nuc2[which(nuc2$chrom1==contig),]$end1
  strata[which(strata$contig==contig),]$start=end-strata[which(strata$contig==contig),]$start
  strata[which(strata$contig==contig),]$end=end-strata[which(strata$contig==contig),]$end

}


# ###CHANGE THE CONTIG NAMES###
new_names <- c(
  "h2tg000005l"="a1-PR",
  "h2tg000018l"="a1-HD",
  "h1tg000040l"="a2-HD",
  "h1tg000006l"="a2-PR"
)
d=mutate(d,chrom_lag=new_names[chrom_lag])
d=mutate(d,chrom1=new_names[chrom1])
nuc1=mutate(nuc1,chrom_lag=new_names[chrom_lag])
nuc2=mutate(nuc2,chrom1=new_names[chrom1])
contigs=mutate(contigs,Chr=new_names[Chr])
data_centromeres=mutate(data_centromeres,contig=new_names[contig])
data_telomeres=mutate(data_telomeres,contig=new_names[contig])
data_genes=mutate(data_genes,contig=new_names[contig])
data_pheromones=mutate(data_pheromones,contig=new_names[contig])
strata=mutate(strata,contig=new_names[contig])


#make a matrix with the contigs start and end to initialize the circos
nb_contig=nrow(contigs)
m=matrix(c(rep(0, nb_contig), c(contigs$End)), ncol=2)


pdf(file = paste("/Users/eliselucotte/desktop/server0/Assembly/ortholog_reconstruction/3.circos/Mv-Cat/A1vsA2/circos_",species1,"-",species2,"_stratacolor.pdf",sep=''))#, width=200, height=200)

##START OF CIRCOS CODE##
##----------------------
circos.clear()
#color of the name of the contigs
col_text <- "grey40"
circos.par("track.height"=0.8, "canvas.xlim"=c(-1.1,1.1),"canvas.ylim"=c(-1.1,1.1),gap.degree=5, cell.padding=c(0, 0, 0, 0))
circos.initialize(factors=contigs$Chr, 
                  xlim=m)

##Plot the contigs
##-----------------
#definition of the color of each contig (following the order of the plot)
#boder_color=c("blue","purple","purple","mediumpurple","dodgerblue")
#border_color=c("blue","purple","purple","purple","blue")

#contig_color=c("#B2D0E6","#BFB9DC","#BFB9DC","#E1DEFC","#E1DEFC","#DEEDFF")
contig_color=c("lightgrey","lightgrey")
circos.track(ylim=c(0, 1), panel.fun=function(x, y) {
  chr=CELL_META$sector.index
  xlim=CELL_META$xlim
  ylim=CELL_META$ylim
  circos.text(mean(xlim), max(ylim)+1, chr, cex=1, col='black', 
              facing="bending.inside", niceFacing=TRUE)
  
}, bg.col=contig_color, bg.border="grey40", track.height=0.09)

## Trace genomes x axis
##----------------------
#size of the biggest genome rounded
max_size=ceiling(max(c(contig1$End/10^6,contig2$End/10^6)))
brk <- seq(0,max_size,1)*10^6

circos.track(track.index = get.current.track.index(), panel.fun=function(x, y) {
  circos.axis(h="top", major.at=brk, labels=round(brk/10^6, 1), labels.cex=0.8,
              col=col_text, labels.col=col_text, lwd=0.7, labels.facing="clockwise")
}, bg.border=F)

##Plot the positions of the genes HD and PR
##-----------------------------------------
for(i in data_genes$contig)
{print(i)
  circos.genomicRect(data_genes[which(data_genes$contig==i),], sector.index=i, track.index=1, ytop = 1, ybottom = 0,col='red',border='red')}

##Plot the positions of the centromeres
##-------------------------------------
for(i in data_centromeres$contig)
{print(i)
  circos.genomicRect(data_centromeres[which(data_centromeres$contig==i),], sector.index=i, track.index=1, ytop = 1, ybottom = 0,col='yellow',border='grey40')}

##Plot the positions of the telomeres
##-------------------------------------
for(i in data_telomeres$contig)
{print(i)
  circos.genomicRect(data_telomeres[which(data_telomeres$contig==i),], sector.index=i, track.index=1, ytop = 1, ybottom = 0,col='forestgreen',border='forestgreen')}

##Plot the positions of the pheromones
##-------------------------------------
for(i in data_pheromones$contig)
{print(i)
  circos.genomicRect(data_pheromones[which(data_pheromones$contig==i),], sector.index=i, track.index=1, ytop = 1, ybottom = 0,col='orange',border='orange')}

##Plot the positions of the strata 
circos.track(ylim=c(0, 1), panel.fun=function(x, y) {
  chr=CELL_META$sector.index
  xlim=CELL_META$xlim
  ylim=CELL_META$ylim
}, bg.col="white", bg.border="grey40", track.height=0.06)

for(col in unique(strata$colors)) {
  print(col)
  data_strata=strata[which(strata$colors==col),]
  for(i in unique(data_strata$contig))
  {
    circos.genomicRect(data_strata[which(data_strata$contig==i),], sector.index=i, track.index=2, ytop = 1, ybottom = 0,col=col,border=col)}
}


##Genomic rearrangements
##------------------------
##to color the inversion differently
#rcols <- scales::alpha(ifelse(sign(nuc1$start_lag - nuc1$end_lag) != sign(nuc2$start1 - nuc2$end1), "#f46d43", "#66c2a5"), alpha=1)
rcols=scales::alpha(ifelse(d$chrom_lag=='a1-HD',"blue","purple"),alpha=1)
circos.genomicLink(nuc1, nuc2, col=rcols, border=NA)



##END OF THE CIRCOS CODE#
##----------------------#
dev.off()


