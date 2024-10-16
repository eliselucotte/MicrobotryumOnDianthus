###PLOT THE DS AND THE STRATA ON THE ACTUAL ORDER 
library(ggplot2)
library(patchwork)

setwd('/Users/eliselucotte/desktop/server1/Assembly/ortholog_reconstruction/Scorzo/dSplot/')

d3=read.table('database_for_strata_M-scorzo_allcontig.txt', as.is=T, header=T)

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
#inverse like in the synteny plot
to_inv_A2=c("MscorzoA2_tig00000008","MscorzoA2_tig00000533","MscorzoA2_tig00000531","MscorzoA2_tig00005667","MscorzoA2_tig00005663")

to_inv_A1=c("MscorzoA1_tig00000008")

for (contig in to_inv_A2)
{
  end=max(d3[which(d3$contig_A2==contig),]$rank_A2)
  #change the coordinates for the synteny databases
  d3[which(d3$contig_A2==contig),]$rank_A2=end-d3[which(d3$contig_A2==contig),]$rank_A2
}
for (contig in to_inv_A1)
{
  end=max(d3[which(d3$contig_A1==contig),]$rank_A1)
  #change the coordinates for the synteny databases
  d3[which(d3$contig_A1==contig),]$rank_A1=end-d3[which(d3$contig_A1==contig),]$rank_A1
}

rank_PR2=max(d3[which(d3$contig_A1=="MscorzoA1_tig00000008"),]$rank_A1)
#add the rank of PR2 to PR1
d3[which(d3$contig_A1=="MscorzoA1_tig00000006"),]$rank_A1=d3[which(d3$contig_A1=="MscorzoA1_tig00000006"),]$rank_A1+rank_PR2
#add the rank of PR1 to PR3
rank_PR1=max(d3[which(d3$contig_A1=="MscorzoA1_tig00000006"),]$rank_A1)
d3[which(d3$contig_A1=="MscorzoA1_tig00000064"),]$rank_A1=d3[which(d3$contig_A1=="MscorzoA1_tig00000064"),]$rank_A1+rank_PR1

rank_PR2=max(d3[which(d3$contig_A2=="MscorzoA2_tig00000008"),]$rank_A2)
#add the rank of PR2 to PR1
d3[which(d3$contig_A2=="MscorzoA2_tig00000533"),]$rank_A2=d3[which(d3$contig_A2=="MscorzoA2_tig00000533"),]$rank_A2+rank_PR2
#add the rank of PR1 to PR3
rank_PR1=max(d3[which(d3$contig_A2=="MscorzoA2_tig00000533"),]$rank_A2)
d3[which(d3$contig_A2=="MscorzoA2_tig00000531"),]$rank_A2=d3[which(d3$contig_A2=="MscorzoA2_tig00000531"),]$rank_A2+rank_PR1


plot_A1PR=ggplot(data=d3[which(d3$contig_A1 %in% c("MscorzoA1_tig00000006",'MscorzoA1_tig00000008','MscorzoA1_tig00000064')),], aes(x=rank_A1, y=dS))+  
  #geom_rect(mapping=aes(xmin=1248, xmax=1449, ymin=0, ymax=0.3), alpha=0.01, fill="lightgoldenrod1")+
  geom_point(aes(color=as.factor(colors)))+
  scale_colour_identity()+
  theme_bw(base_size=20)+
  theme(plot.title = element_text(hjust=0.5))+
  labs(title=bquote("PR "*a[1]), x='gene rank')+
  ylim(0,0.3)
  #geom_vline(aes(xintercept= 2234), col='red')


plot_A2PR=ggplot(data=d3[which(d3$contig_A2 %in% c("MscorzoA2_tig00000533","MscorzoA2_tig00000008","MscorzoA2_tig00000531")),], aes(x=rank_A2, y=dS))+  
  #geom_rect(mapping=aes(xmin=1248, xmax=1822, ymin=0, ymax=0.3), alpha=0.01, fill="lightgoldenrod1")+
  geom_point(aes(color=as.factor(colors)))+
  scale_colour_identity()+
  theme_bw(base_size=20)+
  theme(plot.title = element_text(hjust=0.5))+
  labs(title=bquote("PR "*a[2]), x='gene rank')+
  ylim(0,0.3)
  #geom_vline(aes(xintercept= 2952), col='red')

(plot_A1PR+plot_A2PR)+plot_layout(nrow=2)+plot_annotation(tag_levels = 'A')


ggsave("dS_plot_actual_order_Mscorzo_PR_panels.pdf", device="pdf", height=10, width=15)
ggsave("dS_plot_actual_order_Mscorzo_PR_panels.png", device="png", height=10, width=15)

