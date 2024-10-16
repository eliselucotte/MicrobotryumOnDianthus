# A script to analyse the TE proportion in the mycrobotryum genomes.
library(dplyr)
library(viridis)
library(ggplot2)
library(cowplot)
library(ggnewscale)
library(patchwork)

ThemeSobr=  theme( #Theme for the plot
    panel.border = element_blank(),  
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    text = element_text(size=10, face="bold"),
    axis.line = element_line(colour = "black", size=0.5),
    legend.spacing.y= unit(0, 'cm'),
    axis.title = element_text(face="bold"),
    plot.title=element_text(size=12, face="bold",hjust=0.5, vjust=2),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    axis.text = element_text(size = 10),
    legend.background = element_blank(),
    panel.grid.major.x = element_line(size = 0.2, color="grey50", linetype="dashed"),
    strip.text = element_text(face="bold", size=14)
)

SLScaffoldPR=c("MvDcarth-C212_tig00000115","MvDcarth-C212_tig00000126","MvDcarth-C212_tig00000121", "MvDcarth-C212_tig00000118",
               "MvDelt-1554_tig00000104","MvDelt-1554_tig00000105", "MvDelt-1554_tig00000103",
               "MvDp-1065-A1_tig00000080","MvDp-1065-A1_tig00000084","MvDp-1065-A1_tig00000082",
               "MvDp-1065-A2_tig00000068") #Name of the scaffold linked to PR
SLScaffoldHD=c("MvDcarth-C212_tig00000068", "MvDelt-1554_tig00000059","MvDp-1065-A1_tig00000040","MvDp-1065-A2_tig00000031") #name of the scaffold linked to HD
sexchrom=c(SLScaffoldHD,SLScaffoldPR)
setwd('/Users/eliselucotte/desktop/server1/Assembly/sexchrom/TE_Paul/retedata')
strata=read.table("Position_strata_Mv-sup-1065.txt", header=T, stringsAsFactors = F)

##A1
file="MvDp-1065-A1.gff3"
gff=read.table(file,stringsAsFactors = F) #Read the file
colnames(gff)=c("Scaffold", "Name", "V3", "Start", "End", "V6", "V7", "V8","Annotation" ) 

gff$size=gff$End-gff$Start #Size of the TE
gff$Pos="Autosomal"
MTchrom=c("Mv-sup-1065-A1_tig00000084","Mv-sup-1065-A1_tig00000080","Mv-sup-1065-A1_tig00000040","Mv-sup-1065-A2_tig00000031","Mv-sup-1065-A2_tig00000064")
dchrom <- c(
  "MvDp-1065-A1_tig00000084"="Mv-sup-1065-A1_tig00000084",
  "MvDp-1065-A1_tig00000080"="Mv-sup-1065-A1_tig00000080",
  "MvDp-1065-A1_tig00000040"="Mv-sup-1065-A1_tig00000040",
  "MvDp-1065-A2_tig00000031"="Mv-sup-1065-A2_tig00000031",
  "MvDp-1065-A2_tig00000068"="Mv-sup-1065-A2_tig00000068")

gff= gff %>% mutate(gff,Scaffold2=ifelse(Scaffold %in% sexchrom, dchrom[Scaffold], Scaffold))

for (i in 1:nrow(gff))
{
  if (gff$Scaffold[i] %in% c("MvDp-1065-A1_tig00000080","MvDp-1065-A1_tig00000084","MvDp-1065-A1_tig00000082","MvDp-1065-A2_tig00000068"))
  {
    gff$Pos[i]="Black PR" #By default, the region in the PR scaffold are considered to be part of the black stratum
  }
  if (gff$Scaffold[i] %in% c("MvDp-1065-A1_tig00000040","MvDp-1065-A2_tig00000031"))
  {
    gff$Pos[i]="HD chromosome" #By default, the region in the PR scaffold are considered to be part of the black stratum
  }
  if (gff$Scaffold2[i] %in% strata$contig)
  {
    y=which(strata$contig==gff$Scaffold2[i])
    for (Strata in y) #All other region belong to other strata or the PAR
    {
      if (strata$start[Strata] < gff$Start[i] & strata$end[Strata] > gff$End[i])
      {
        gff$Pos[i]=strata$strata[Strata]
      }
    }
  }
}
gff[gff$Pos=="Autosomal",]$Pos= gff[gff$Pos=="Autosomal",]$Scaffold
gff$Pos=gsub("MvDp-1065-A2_tig00000", "Autosome ",gff$Pos)
gff$Pos=gsub("MvDp-1065-A1_tig00000", "Autosome ",gff$Pos)

gffSum=gff %>% group_by(Scaffold2, Pos) %>% summarise(SumTE=sum(size), StartSc=min(Start), EndSc=max(End)) #Summarise by scaffold. Here, the size of the scaffold is determined by taking the end porition of the latest TE. This is suboptimal (there is still some dna sites after the last TE in the scaffold), but in practice very close to the reality
gffSum=as.data.frame(gffSum)
for (i  in 1:nrow(strata))
{
  if (strata$strata[i] %in% gffSum$Pos) # redefine the start and end of strata
  {
    y=which(gffSum$Pos==strata$strata[i])
    gffSum$StartSc[y]=strata$start[i]
    gffSum$EndSc[y]=strata$end[i]
  }
  else  #Fill the gffSum table with strata that are empty of TE
  {
    if (strata$contig[i] %in% gffSum$Scaffold)
    { 
      contig=strata$contig[i]
      strate=strata$strata[i]
      sizeTe=0
      start=strata$start[i]
      end=strata$end[i]
      gffSum[nrow(gffSum)+1,]= c(contig, strate, sizeTe, start, end)
    }
  }
}

gffSum$EndSc=as.numeric(gffSum$EndSc)
gffSum$StartSc=as.numeric(gffSum$StartSc)
gffSum$SumTE=as.numeric(gffSum$SumTE)
gffSum$Size=gffSum$EndSc - gffSum$StartSc #Size of the scaffols
gffSum[(grepl("dark_ruby", as.character(gffSum$Pos))),]$Pos="Dark ruby"
gffSum[(grepl("turquoise", as.character(gffSum$Pos))),]$Pos="Turquoise"
gffSum[(grepl("light_ruby", as.character(gffSum$Pos))),]$Pos="Light ruby"
gffSum[(grepl("PAR", as.character(gffSum$Pos))),]$Pos="PAR"
gffSum[(grepl("blue", as.character(gffSum$Pos))),]$Pos="Blue"
gffSum[(grepl("blackHD", as.character(gffSum$Pos))),]$Pos="Black HD"
gffSum=gffSum %>% group_by(Pos) %>% summarise(SumTE=sum(SumTE), Size=sum(Size)) #merge the strata being on diff chromosomes
gffSum$TEProp=gffSum$SumTE / gffSum$Size #overall proportion of TE. This is calculated only to order scaffold by TE proportion
colnames(gffSum)=c("Strata","SumTE","Size","TEProp")
gffSum$Age=1
gffSum[gffSum$Strata=="Black PR",]$Age=10
gffSum[gffSum$Strata=="Black HD",]$Age=3
gffSum[gffSum$Strata=="HD chromosome",]$Age=3
gffSum[gffSum$Strata=="Dark ruby",]$Age=8
gffSum[gffSum$Strata=="Light ruby",]$Age=7
gffSum[gffSum$Strata=="Turquoise",]$Age=6
gffSum[gffSum$Strata=="Blue",]$Age=5
gffSum[gffSum$Strata=="PAR",]$Age=4

gffSum$Strata <- factor(gffSum$Strata, levels = gffSum$Strata[order(gffSum$Age, gffSum$TEProp)]) #Modify factor order for nice plotting (sorted by TE proportion)

Unclass=c("unclassified:unclassified:unclassified", "?:?:?") #Change the annotation 
gff[gff$Annotation %in% Unclass,]$Annotation="Unclassified"
gff[gff$Annotation=="?:?:Simple",]$Annotation="Simple repeat"
gff[gff$Annotation=="?:?:Satellite",]$Annotation="Satellite"
gff[gff$Annotation=="?:?:rRNA",]$Annotation="rRNA"

ColLab=c("white",scales::viridis_pal(begin=0.0, end=0.8, option = "A", direction = 1)(2)) #Color of text label
Col=c("#729b5f", #Color of Annotation
      "#4f7dd6",
      "#87ce3c",
      "#68a0c6",
      "#71b453",
      "#6e9fa8",
      "#57c8a3",
      "#759681",
      "#d9b442",
      "#b352cd",
      "#c68343",
      "#d5488b",
      "#97785c",
      "#d94c2e",
      "#9f6da2",
      "#c2585d",
      "#daa9ba",
      "#417ab1",
      "#0d1724",
      "#467d8b",
      "#1b3344",
      "#3b586f")
gff[(grepl("dark_ruby", as.character(gff$Pos))),]$Pos="Dark ruby"
gff[(grepl("turquoise", as.character(gff$Pos))),]$Pos="Turquoise"
gff[(grepl("light_ruby", as.character(gff$Pos))),]$Pos="Light ruby"
gff[(grepl("PAR", as.character(gff$Pos))),]$Pos="PAR"
gff[(grepl("blue", as.character(gff$Pos))),]$Pos="Blue"
gff[(grepl("blackHD", as.character(gff$Pos))),]$Pos="Black HD"

gffSum2=gff %>% group_by(Pos,Annotation) %>% summarise(SumTE=sum(size)) #Genome span by each categories of TE
colnames(gffSum2)=c("Strata","Annotation","SumTE")
gffSum2= gffSum2 %>% inner_join(gffSum, by="Strata",)

gffSum2$TEProp=gffSum2$SumTE.x /gffSum2$Size #Proportion of each categories of TE
gffSum2$Age=1
gffSum2[gffSum2$Strata=="Black PR",]$Age=10
gffSum2[gffSum2$Strata=="Black HD",]$Age=3
gffSum2[gffSum2$Strata=="Dark ruby",]$Age=8
gffSum2[gffSum2$Strata=="Light ruby",]$Age=7
gffSum2[gffSum2$Strata=="Turquoise",]$Age=6
gffSum2[gffSum2$Strata=="Blue",]$Age=5
gffSum2[gffSum2$Strata=="PAR",]$Age=4
gffSum2[gffSum2$Strata=="HD chromosome",]$Age=3
gffSum2$Strata <- factor(gffSum2$Strata, levels = gffSum$Strata[order(gffSum$Age, gffSum$TEProp)]) #Reorder factor for nice plotting
# gffSum2$ColorLabel= ifelse(gffSum2$Scaffold %in% SLScaffoldPR, "PR", ifelse(gffSum2$Scaffold %in% SLScaffoldHD, "HD", "")) #Color of label
# gffSum$ColorLabel= ifelse(gffSum$Scaffold %in% SLScaffoldPR, "PR", ifelse(gffSum$Scaffold %in% SLScaffoldHD, "HD", ""))
gffSum2Sub=gffSum2[!(grepl("Autosome", as.character(gffSum2$Strata)) & gffSum2$Size<300000),]
gffSumSub=gffSum[!(grepl("Autosome", as.character(gffSum$Strata)) & gffSum$Size<300000),]
# gffSum2Sub[(gffSum2Sub$Pos=="Autosomal" &  gffSum2Sub$Size<300000),]

# gffSum2Sub=gffSum2Sub[order(gffSum$TEProp),]
# labelColor= ifelse(gffSum2Sub$Pos == "BlackPR", "black", ifelse(gffSum2Sub$Pos == "Autosomal", "grey50", "red"))
# gffSum2Sub$labelColor=labelColor
# 
# labelColor= ifelse(gffSum2Sub$Strata == "MvDp-1065-A2_tig00000068_BlackPR", "black", ifelse(gffSum2Sub$Strata == "MvDp-1065-A2_tig00000023_Autosomal", "grey50", "red"))

size=2

base=ggplot(gffSum2Sub) #Only consider scaffold larger than 300000 bp
PlotA1=base+geom_bar(aes(x=Strata, y=TEProp, fill=Annotation), #Plot of scaffold proportions
                   stat='identity', position="stack", color="black", size=0.05)+
  coord_flip()+
  scale_fill_manual("Annotation", values=Col)+
  #geom_point(aes(x=Scaffold, y=0.75, color=ColorLabel), shape=8)+
  # geom_text(aes(x=Scaffold, y=0.75, label=ColorLabel, color=ColorLabel), size=2)+
  scale_color_manual("", values=ColLab, guide=F)+
  ylab("Proportion of TEs")+xlab("Region")+ThemeSobr
# Plot
# save_plot(paste0(file,'cleanV3.png'), Plot, base_aspect_ratio = 2)
# save_plot(paste0(file,'cleanV3.pdf'), Plot, base_aspect_ratio = 2)


base=ggplot(gffSum2Sub)
PlotBpA1=base+ #Plot of scaffold size (in bp) and TE proportion. To do so, plot a white bar that represent scaffold size and plot over it the bar that represent TE proportion
  geom_bar(data=gffSumSub,aes(x=Strata, y=Size, fill="Non TE/repeat (e.g. genes, intergenic sequences)"),
           stat='identity', color="black", size=0.1)+
  scale_fill_manual("", values="white")+
  new_scale_fill()+
  geom_bar(aes(x=Strata, y=SumTE.x, fill=Annotation),
           stat='identity', color="black", size=0.05)+ #Plot TE proportion
  # geom_text(data=gffSum[gffSum$ScSize>300000,], aes(x=Scaffold, y=max(ScSize) + 0.10*max(ScSize), label=ColorLabel, color=ColorLabel), size=2)+
  scale_color_manual("", values=ColLab, guide=F)+
  coord_flip()+
  scale_fill_manual("", values=Col)+
  ylab(" Size spanned (bp)")+xlab("Region")+ThemeSobr+ylim(0,max(gffSum2$Size)+0.15*max(gffSum2$Size))
# save_plot(paste0(file,'.GenomeSize.cleanV3.png'), PlotBp, base_aspect_ratio = 2)
# save_plot(paste0(file,'.GenomeSize.cleanV3.pdf'), PlotBp, base_aspect_ratio = 2)
write.table(gffSum, paste0(file,'.summary.txt'),  quote = F, row.names = F, sep="\t")

#### SummaryForPaper
SexAutoSum=gffSum
SexAutoSum$Pos="Mating-type chromosome PR"
SexAutoSum[grepl("Autosome", SexAutoSum$Strata),]$Pos="Autosome"
SexAutoSum[grepl("HD", SexAutoSum$Strata),]$Pos="Mating-type chromosome HD"
SexAutoSum= SexAutoSum %>% group_by(Pos) %>% summarise(SumTE=sum(SumTE), Size=sum(Size)) %>% mutate(TEProp=SumTE/Size)
write.table(SexAutoSum, paste0(file,'.summarySexAuto.txt'),  quote = F, row.names = F, sep="\t")


##A2

file="MvDp-1065-A2.gff3"

gff=read.table(file,stringsAsFactors = F) #Read the file
colnames(gff)=c("Scaffold", "Name", "V3", "Start", "End", "V6", "V7", "V8","Annotation" ) 

gff$size=gff$End-gff$Start #Size of the TE
gff$Pos="Autosomal"
MTchrom=c("Mv-sup-1065-A1_tig00000084","Mv-sup-1065-A1_tig00000080","Mv-sup-1065-A1_tig00000040","Mv-sup-1065-A2_tig00000031","Mv-sup-1065-A2_tig00000064")
dchrom <- c(
  "MvDp-1065-A1_tig00000084"="Mv-sup-1065-A1_tig00000084",
  "MvDp-1065-A1_tig00000080"="Mv-sup-1065-A1_tig00000080",
  "MvDp-1065-A1_tig00000040"="Mv-sup-1065-A1_tig00000040",
  "MvDp-1065-A2_tig00000031"="Mv-sup-1065-A2_tig00000031",
  "MvDp-1065-A2_tig00000068"="Mv-sup-1065-A2_tig00000068")

gff= gff %>% mutate(gff,Scaffold2=ifelse(Scaffold %in% sexchrom, dchrom[Scaffold], Scaffold))

for (i in 1:nrow(gff))
{
  if (gff$Scaffold[i] %in% c("MvDp-1065-A1_tig00000080","MvDp-1065-A1_tig00000084","MvDp-1065-A1_tig00000082","MvDp-1065-A2_tig00000068"))
  {
    gff$Pos[i]="Black PR" #By default, the region in the PR scaffold are considered to be part of the black stratum
  }
  if (gff$Scaffold[i] %in% c("MvDp-1065-A1_tig00000040","MvDp-1065-A2_tig00000031"))
  {
    gff$Pos[i]="HD chromosome" #By default, the region in the PR scaffold are considered to be part of the black stratum
  }
  if (gff$Scaffold2[i] %in% strata$contig)
  {
    y=which(strata$contig==gff$Scaffold2[i])
    for (Strata in y) #All other region belong to other strata or the PAR
    {
      if (strata$start[Strata] < gff$Start[i] & strata$end[Strata] > gff$End[i])
      {
        gff$Pos[i]=strata$strata[Strata]
      }
    }
  }
}
gff[gff$Pos=="Autosomal",]$Pos= gff[gff$Pos=="Autosomal",]$Scaffold
gff$Pos=gsub("MvDp-1065-A2_tig00000", "Autosome ",gff$Pos)
gff$Pos=gsub("MvDp-1065-A1_tig00000", "Autosome ",gff$Pos)

gffSum=gff %>% group_by(Scaffold2, Pos) %>% summarise(SumTE=sum(size), StartSc=min(Start), EndSc=max(End)) #Summarise by scaffold. Here, the size of the scaffold is determined by taking the end porition of the latest TE. This is suboptimal (there is still some dna sites after the last TE in the scaffold), but in practice very close to the reality
gffSum=as.data.frame(gffSum)
for (i  in 1:nrow(strata))
{
  if (strata$strata[i] %in% gffSum$Pos) # redefine the start and end of strata
  {
    y=which(gffSum$Pos==strata$strata[i])
    gffSum$StartSc[y]=strata$start[i]
    gffSum$EndSc[y]=strata$end[i]
  }
  else  #Fill the gffSum table with strata that are empty of TE
  {
    if (strata$contig[i] %in% gffSum$Scaffold)
    { 
      contig=strata$contig[i]
      strate=strata$strata[i]
      sizeTe=0
      start=strata$start[i]
      end=strata$end[i]
      gffSum[nrow(gffSum)+1,]= c(contig, strate, sizeTe, start, end)
    }
  }
}

gffSum$EndSc=as.numeric(gffSum$EndSc)
gffSum$StartSc=as.numeric(gffSum$StartSc)
gffSum$SumTE=as.numeric(gffSum$SumTE)
gffSum$Size=gffSum$EndSc - gffSum$StartSc #Size of the scaffols
gffSum[(grepl("dark_ruby", as.character(gffSum$Pos))),]$Pos="Dark ruby"
gffSum[(grepl("turquoise", as.character(gffSum$Pos))),]$Pos="Turquoise"
gffSum[(grepl("light_ruby", as.character(gffSum$Pos))),]$Pos="Light ruby"
gffSum[(grepl("PAR", as.character(gffSum$Pos))),]$Pos="PAR"
gffSum[(grepl("blue", as.character(gffSum$Pos))),]$Pos="Blue"
gffSum[(grepl("blackHD", as.character(gffSum$Pos))),]$Pos="Black HD"
gffSum=gffSum %>% group_by(Pos) %>% summarise(SumTE=sum(SumTE), Size=sum(Size)) #merge the strata being on diff chromosomes
gffSum$TEProp=gffSum$SumTE / gffSum$Size #overall proportion of TE. This is calculated only to order scaffold by TE proportion
colnames(gffSum)=c("Strata","SumTE","Size","TEProp")
gffSum$Age=1
gffSum[gffSum$Strata=="Black PR",]$Age=10
gffSum[gffSum$Strata=="Black HD",]$Age=3
gffSum[gffSum$Strata=="HD chromosome",]$Age=3
gffSum[gffSum$Strata=="Dark ruby",]$Age=8
gffSum[gffSum$Strata=="Light ruby",]$Age=7
gffSum[gffSum$Strata=="Turquoise",]$Age=6
gffSum[gffSum$Strata=="Blue",]$Age=5
gffSum[gffSum$Strata=="PAR",]$Age=4

gffSum$Strata <- factor(gffSum$Strata, levels = gffSum$Strata[order(gffSum$Age, gffSum$TEProp)]) #Modify factor order for nice plotting (sorted by TE proportion)

Unclass=c("unclassified:unclassified:unclassified", "?:?:?") #Change the annotation 
gff[gff$Annotation %in% Unclass,]$Annotation="Unclassified"
gff[gff$Annotation=="?:?:Simple",]$Annotation="Simple repeat"
gff[gff$Annotation=="?:?:Satellite",]$Annotation="Satellite"
gff[gff$Annotation=="?:?:rRNA",]$Annotation="rRNA"

ColLab=c("white",scales::viridis_pal(begin=0.0, end=0.8, option = "A", direction = 1)(2)) #Color of text label
Col=c("#729b5f", #Color of Annotation
      "#4f7dd6",
      "#87ce3c",
      "#68a0c6",
      "#71b453",
      "#6e9fa8",
      "#57c8a3",
      "#759681",
      "#d9b442",
      "#b352cd",
      "#c68343",
      "#d5488b",
      "#97785c",
      "#d94c2e",
      "#9f6da2",
      "#c2585d",
      "#daa9ba",
      "#417ab1",
      "#0d1724",
      "#467d8b",
      "#1b3344",
      "#3b586f")
gff[(grepl("dark_ruby", as.character(gff$Pos))),]$Pos="Dark ruby"
gff[(grepl("turquoise", as.character(gff$Pos))),]$Pos="Turquoise"
gff[(grepl("light_ruby", as.character(gff$Pos))),]$Pos="Light ruby"
gff[(grepl("PAR", as.character(gff$Pos))),]$Pos="PAR"
gff[(grepl("blue", as.character(gff$Pos))),]$Pos="Blue"
gff[(grepl("blackHD", as.character(gff$Pos))),]$Pos="Black HD"

gffSum2=gff %>% group_by(Pos,Annotation) %>% summarise(SumTE=sum(size)) #Genome span by each categories of TE
colnames(gffSum2)=c("Strata","Annotation","SumTE")
gffSum2= gffSum2 %>% inner_join(gffSum, by="Strata",)

gffSum2$TEProp=gffSum2$SumTE.x /gffSum2$Size #Proportion of each categories of TE
gffSum2$Age=1
gffSum2[gffSum2$Strata=="Black PR",]$Age=10
gffSum2[gffSum2$Strata=="Black HD",]$Age=3
gffSum2[gffSum2$Strata=="Dark ruby",]$Age=8
gffSum2[gffSum2$Strata=="Light ruby",]$Age=7
gffSum2[gffSum2$Strata=="Turquoise",]$Age=6
gffSum2[gffSum2$Strata=="Blue",]$Age=5
gffSum2[gffSum2$Strata=="PAR",]$Age=4
gffSum2[gffSum2$Strata=="HD chromosome",]$Age=3
gffSum2$Strata <- factor(gffSum2$Strata, levels = gffSum$Strata[order(gffSum$Age, gffSum$TEProp)]) #Reorder factor for nice plotting
# gffSum2$ColorLabel= ifelse(gffSum2$Scaffold %in% SLScaffoldPR, "PR", ifelse(gffSum2$Scaffold %in% SLScaffoldHD, "HD", "")) #Color of label
# gffSum$ColorLabel= ifelse(gffSum$Scaffold %in% SLScaffoldPR, "PR", ifelse(gffSum$Scaffold %in% SLScaffoldHD, "HD", ""))
gffSum2Sub=gffSum2[!(grepl("Autosome", as.character(gffSum2$Strata)) & gffSum2$Size<300000),]
gffSumSub=gffSum[!(grepl("Autosome", as.character(gffSum$Strata)) & gffSum$Size<300000),]
# gffSum2Sub[(gffSum2Sub$Pos=="Autosomal" &  gffSum2Sub$Size<300000),]

# gffSum2Sub=gffSum2Sub[order(gffSum$TEProp),]
# labelColor= ifelse(gffSum2Sub$Pos == "BlackPR", "black", ifelse(gffSum2Sub$Pos == "Autosomal", "grey50", "red"))
# gffSum2Sub$labelColor=labelColor
# 
# labelColor= ifelse(gffSum2Sub$Strata == "MvDp-1065-A2_tig00000068_BlackPR", "black", ifelse(gffSum2Sub$Strata == "MvDp-1065-A2_tig00000023_Autosomal", "grey50", "red"))

base=ggplot(gffSum2Sub) #Only consider scaffold larger than 300000 bp
PlotA2=base+geom_bar(aes(x=Strata, y=TEProp, fill=Annotation), #Plot of scaffold proportions
                   stat='identity', position="stack", color="black", size=0.05)+
  coord_flip()+
  scale_fill_manual("Annotation", values=Col)+
  #geom_point(aes(x=Scaffold, y=0.75, color=ColorLabel), shape=8)+
  # geom_text(aes(x=Scaffold, y=0.75, label=ColorLabel, color=ColorLabel), size=2)+
  scale_color_manual("", values=ColLab, guide=F)+
  ylab("Proportion of TEs")+xlab("Region")+ThemeSobr

# save_plot(paste0(file,'cleanV3.png'), Plot, base_aspect_ratio = 2)
# save_plot(paste0(file,'cleanV3.pdf'), Plot, base_aspect_ratio = 2)
# 

base=ggplot(gffSum2Sub)
PlotBpA2=base+ #Plot of scaffold size (in bp) and TE proportion. To do so, plot a white bar that represent scaffold size and plot over it the bar that represent TE proportion
  geom_bar(data=gffSumSub,aes(x=Strata, y=Size, fill="Non TE/repeat (e.g. genes, intergenic sequences)"),
           stat='identity', color="black", size=0.1)+
  scale_fill_manual("", values="white")+
  new_scale_fill()+
  geom_bar(aes(x=Strata, y=SumTE.x, fill=Annotation),
           stat='identity', color="black", size=0.05)+ #Plot TE proportion
  # geom_text(data=gffSum[gffSum$ScSize>300000,], aes(x=Scaffold, y=max(ScSize) + 0.10*max(ScSize), label=ColorLabel, color=ColorLabel), size=2)+
  scale_color_manual("", values=ColLab, guide=F)+
  coord_flip()+
  scale_fill_manual("", values=Col)+
  ylab(" Size spanned (bp)")+xlab("Region")+ThemeSobr+ylim(0,max(gffSum2$Size)+0.15*max(gffSum2$Size))
# save_plot(paste0(file,'.GenomeSize.cleanV3.png'), PlotBp, base_aspect_ratio = 2)
# save_plot(paste0(file,'.GenomeSize.cleanV3.pdf'), PlotBp, base_aspect_ratio = 2)
write.table(gffSum, paste0(file,'.summary.txt'),  quote = F, row.names = F, sep="\t")

#### SummaryForPaper
SexAutoSum=gffSum
SexAutoSum$Pos="Mating-type chromosome PR"
SexAutoSum[grepl("Autosome", SexAutoSum$Strata),]$Pos="Autosome"
SexAutoSum[grepl("HD", SexAutoSum$Strata),]$Pos="Mating-type chromosome HD"
SexAutoSum= SexAutoSum %>% group_by(Pos) %>% summarise(SumTE=sum(SumTE), Size=sum(Size)) %>% mutate(TEProp=SumTE/Size)
write.table(SexAutoSum, paste0(file,'.summarySexAuto.txt'),  quote = F, row.names = F, sep="\t")



patchworkplot=PlotA1+PlotA2+PlotBpA1+PlotBpA2+  plot_layout(nrow=2, ncol=2, guides = 'collect')+plot_annotation(tag_levels = 'A')& theme(plot.tag = element_text(size = 20))

ggsave(plot=patchworkplot,file='PatchworkPlot.pdf',device='pdf', width=15,height=10)

patchworkplot2=PlotA1+PlotBpA1+plot_layout(nrow=2, ncol=1, guides = 'collect')+plot_annotation(tag_levels = 'A')& theme(plot.tag = element_text(size = 20))
patchworkplot2
ggsave(plot=patchworkplot2,file='PatchworkPlot2.pdf',device='pdf', width=10,height=10)
