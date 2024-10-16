###documentation links:
##https://www.royfrancis.com/beautiful-circos-plots-in-r/#data
##https://jokergoo.github.io/circlize_book/book/genomic-plotting-region.html#genomic-rectangles
##centromere file : @GEEserver03 /home/marine/MS_newSpecies/allGenomes_centromereCoordinates_newSpecies.tab
setwd("/Users/eliselucotte/desktop/server0/Assembly/ortholog_reconstruction/Scorzo/circos/")
library(circlize)
library(dplyr)
library(tidyr)
species1=paste("M-scorzo-A1")
species2=paste("M-scorzo-A2")

chroms_lag=c("tig00000022","tig00000123","tig00000006","tig00000064","tig00000008","tig00000163","tig00000036")
chroms_1=c("tig00000185","tig00000008","tig00000533","tig00000531","tig00000007","tig00005667","tig00000010")


##Import the synteny data
d=read.table(paste('synteny_',species1,'_',species2,'.txt',sep=''), header=T, as.is=T, sep='\t')
d=d[which(d$chrom_lag %in% chroms_lag & d$chrom1 %in% chroms2),]

#Import the contig informations
contig1=read.table(paste(species1,'_contig.bed',sep=''), header=T, as.is=T, sep='\t')
contig1=contig1[which(contig1$Chr %in% unique(d$chrom_lag)),]
contig2=read.table(paste(species2,'_contig.bed',sep=''), header=T, as.is=T, sep='\t')
contig2=contig2[which(contig2$Chr %in% unique(d$chrom1)),]

##Reorder the contig in the order I want them to appear in the plot
contig1=contig1[c(5,6,1,4,7,3,2),]
contig2=contig2[c(6,2,1,5,4,3,7),]

#bind both tables
contigs=rbind.data.frame(contig1,contig2)

#Import the gene positions + the centromere positions
data_genes=read.table(paste('Positions_genes.txt',sep=''), header=T, as.is=T, sep='\t')
#data_centromere=read.table(paste('Positions_centromeres.txt',sep=''), header=T, as.is=T, sep='\t')
#remove the centromere for the splitted contigs
#to_rm=c('tig00000084')
#data_centromere=data_centromere[-which(data_centromere$contig %in% to_rm ),]
#Make the database for the genomic links
nuc1=select(d, c('chrom_lag', 'start_lag','end_lag'))
nuc2=select(d, c('chrom1', 'start1','end1'))

##Which contig to you want to invert the orientation for?
to_inv=c("tig00000533","tig00000007","MC02-2","MC02-1")
for (contig in to_inv)
{
  end=contigs[which(contigs$Chr==contig),]$End
  #change the coordinates for the genes
  data_genes[which(data_genes$contig==contig),]$start=end-data_genes[which(data_genes$contig==contig),]$start
  data_genes[which(data_genes$contig==contig),]$end=end-data_genes[which(data_genes$contig==contig),]$end
  #change the coordinates for the centromeres
  data_centromere[which(data_centromere$contig==contig),]$start=end-data_centromere[which(data_centromere$contig==contig),]$start
  data_centromere[which(data_centromere$contig==contig),]$end=end-data_centromere[which(data_centromere$contig==contig),]$end
  #change the coordinates for the synteny databases
  nuc1[which(nuc1$chrom_lag==contig),]$start_lag=end-nuc1[which(nuc1$chrom_lag==contig),]$start_lag
  nuc1[which(nuc1$chrom_lag==contig),]$end_lag=end-nuc1[which(nuc1$chrom_lag==contig),]$end_lag
  nuc2[which(nuc2$chrom1==contig),]$start1=end-nuc2[which(nuc2$chrom1==contig),]$start1
  nuc2[which(nuc2$chrom1==contig),]$end1=end-nuc2[which(nuc2$chrom1==contig),]$end1

}

###CHANGE THE CONTIG NAMES###
new_names <- c(
  "tig00000123"="a1-HD",
  "tig00000006"="a1-PR1",
  "tig00000008"="a1-PR2",
  "tig00000064"="a1-PR3",
  "tig00000163"="a1-PR4",
  "tig00000022"="a1-tig022",
  "tig00000036"="a1-tig036",  
  "tig00000185"="a2-HD",
  "tig00000008"="a2-PR1",
  "tig00000533"="a2-PR2",
  "tig00000531"="a2-PR3",
  "tig00000007"="a2-tig007",
  "tig00005667"="a2-PR4",
  "tig00000010"="a2-tig010"
)
d=mutate(d,chrom_lag=new_names[chrom_lag])
d=mutate(d,chrom1=new_names[chrom1])
nuc1=mutate(nuc1,chrom_lag=new_names[chrom_lag])
nuc2=mutate(nuc2,chrom1=new_names[chrom1])
contigs=mutate(contigs,Chr=new_names[Chr])
# data_centromere=mutate(data_centromere,contig=new_names[contig])
data_genes=mutate(data_genes,contig=new_names[contig])


#make a matrix with the contigs start and end to initialize the circos
nb_contig=nrow(contigs)
m=matrix(c(rep(0, nb_contig), c(contigs$End)), ncol=2)

# png(file = paste("/Users/eliselucotte/desktop/server0/Assembly/ortholog_reconstruction/3.circos/A1vsA2/circos_",species1,"-",species2,".png",sep=''))#, width=200, height=200)

pdf(file = paste("/Users/eliselucotte/desktop/server0/Assembly/ortholog_reconstruction/Scorzo/circos/A1vsA2/circos_A1_vs_A2_stratacolor.pdf",sep=''))#, width=200, height=200)

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
contig_color=c("lightgrey","lightgrey","lightgrey")
circos.track(ylim=c(0, 1), panel.fun=function(x, y) {
  chr=CELL_META$sector.index
  xlim=CELL_META$xlim
  ylim=CELL_META$ylim
  circos.text(mean(xlim), max(ylim)+1, chr, cex=0.5, col='black', 
              facing="bending.inside", niceFacing=TRUE)
  
}, bg.col=contig_color, bg.border="grey40", track.height=0.06)

## Trace genomes x axis
##----------------------
#size of the biggest genome rounded
max_size=ceiling(max(c(contig1$End/10^6,contig2$End/10^6)))
brk <- seq(0,max_size,1)*10^6

circos.track(track.index = get.current.track.index(), panel.fun=function(x, y) {
  circos.axis(h="top", major.at=brk, labels=round(brk/10^6, 1), labels.cex=0.4,
              col=col_text, labels.col=col_text, lwd=0.7, labels.facing="clockwise")
}, bg.border=F)

##Plot the positions of the genes HD and PR
##-----------------------------------------
for(i in data_genes$contig)
{print(i)
  circos.genomicRect(data_genes[which(data_genes$contig==i),], sector.index=i, track.index=1, ytop = 1, ybottom = 0,col='red',border='red')}

##Plot the positions of the centromeres
##-------------------------------------
for(i in data_centromere$contig)
{print(i)
  circos.genomicRect(data_centromere[which(data_centromere$contig==i),], sector.index=i, track.index=1, ytop = 1, ybottom = 0,col='yellow',border='grey40')}


##Genomic rearrangements
##------------------------
##to color the inversion differently
#rcols <- scales::alpha(ifelse(sign(nuc1$start_lag - nuc1$end_lag) != sign(nuc2$start1 - nuc2$end1), "#f46d43", "#66c2a5"), alpha=1)
rcols=scales::alpha(ifelse(d$chrom_lag=='a1-HD',"blue","purple"),alpha=1)
circos.genomicLink(nuc1, nuc2, col=rcols, border=NA)



##END OF THE CIRCOS CODE#
##----------------------#
dev.off()





#dev.copy2pdf(file = paste("/Users/eliselucotte/desktop/server3/Assembly/sexchrom/MvDp1065/synteny_plots/sexchrom_MvSv_orto_MvDp/CIRCOS/circos_",MT,".pdf",sep=''), width=400, height=400)



# circos.genomicTrack(data=tgenes_all, panel.fun=function(x,y){
#   circos.genomicPoints(start, yaxis, pch = 16, col = "red")
# }, track.height=0.08, bg.border=F)

# circos.genomicLabels(data_genes, labels.column=4, cex=1, col=col_text, line_lwd=0.5, line_col="grey80", 
#                      side="inside", connection_height=0.05, labels_height=0.04)

# circos.genomicTrack(tgenes_all,numeric.column = 4, 
#                     panel.fun = function(region, value, ...) {
#                       circos.genomicPoints(region, value,...)
#                     })

## trace genomes x axis
##size of the biggest genome rounded
# max_size=ceiling(max(c(contig1$End/10^6,contig2$End/10^6)))
# brk <- seq(0,max_size,1)*10^6
# 
# circos.track(track.index = get.current.track.index(), panel.fun=function(x, y) {
#   circos.axis(h="top", major.at=brk, labels=round(brk/10^6, 1), labels.cex=0.4, 
#               col=col_text, labels.col=col_text, lwd=0.7, labels.facing="clockwise")
# }, bg.border=F)




