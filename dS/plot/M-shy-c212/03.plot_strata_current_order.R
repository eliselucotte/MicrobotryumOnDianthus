###PLOT THE DS AND THE STRATA ON THE ACTUAL ORDER 
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(stringr)
setwd('/Users/eliselucotte/desktop/server1/Assembly/ortholog_reconstruction/')
d3=read.table('./2.make_DS_plot/find_strata/Mv-cat-C212/database_for_strata_Mv-cat-C212.txt', as.is=T, header=T)
d4=d3%>%select(c("gene_A1","contig_A1","start_A1","end_A1","rank_A1","gene_A2","contig_A2","start_A2","end_A2","rank_A2","dS","strata_lag","colors","strata_name"))

#A1
#h2tg000005l	5364995	5464323
#h2tg000018l	345244	540144
# h2tg000005l	6142584	6143931	PR
# h2tg000018l	870353	871691	HD1
# h2tg000018l	868588	870103	HD2
#telomeres
#h2tg000005l	13270000	13280000
#h2tg000005l	370000	380000

#A2
# h1tg000006l	7626814	7924938
# h1tg000006l	5370586	5371891	PR
# h1tg000040l	58917	60376	HD1
# h1tg000040l	60626	62141	HD2


#get the rank for PR
d4[which(d4$contig_A1=='Mv-cat-C212-A1_h2tg000005l' & d4$start_A1>6142584-372152-100000 & d4$start_A1<6142584-372152+100000),]
#1713
d4[which(d4$contig_A2=='Mv-cat-C212-A2_h1tg000006l' & d4$start_A2>5370586-40000 & d4$start_A2<5370586+15000),]
#1424

#get the rank for HD
d4[which(d4$contig_A1=='Mv-cat-C212-A1_h2tg000018l' & d4$start_A1>870353-10000 & d4$start_A1<870353+10000),]
#245
d4[which(d4$contig_A2=='Mv-cat-C212-A2_h1tg000040l' & d4$start_A2>58917-10000 & d4$start_A2<58917+10000),]
#18


#get ranks for centromere
d4[which(d4$contig_A1=='Mv-cat-C212-A1_h2tg000005l' & d4$start_A1>5364995-372152-400000 & d4$start_A1<5364995-372152+400000),]
d4[which(d4$contig_A1=='Mv-cat-C212-A1_h2tg000005l' & d4$start_A1>5464323-372152-10000 & d4$start_A1<5464323-372152+1000000),]
#1572-1697

		
d4[which(d4$contig_A2=='Mv-cat-C212-A2_h1tg000006l' & d4$start_A2>7626814-100000 & d4$start_A2<7626814+100000),]
d4[which(d4$contig_A2=='Mv-cat-C212-A2_h1tg000006l' & d4$start_A2>7924938-200000 & d4$start_A2<7924938+200000),]
#2013-2146

plot_A1PR=ggplot(data=d4[which(d4$contig_A1 %in% c('Mv-cat-C212-A1_h2tg000005l')),], aes(x=rank_A1, y=dS))+  
  geom_rect(mapping=aes(xmin=1572, xmax=1697, ymin=0, ymax=0.3), alpha=0.01, fill="lightgoldenrod1")+
  geom_point(aes(color=as.factor(colors)))+
  scale_colour_identity()+
  theme_bw(base_size=20)+
  theme(plot.title = element_text(hjust=0.5))+
  labs(title=bquote("PR "*a[1]), x='gene rank')+
  ylim(0,0.3)+
  geom_vline(aes(xintercept= 1713), col='red')



plot_A2PR=ggplot(data=d4[which(d4$contig_A2 %in% c('Mv-cat-C212-A2_h1tg000006l')),], aes(x=rank_A2, y=dS))+  
  geom_rect(mapping=aes(xmin=2013, xmax=2146, ymin=0, ymax=0.3), alpha=0.01, fill="lightgoldenrod1")+
  geom_point(aes(color=as.factor(colors)))+
  scale_colour_identity()+
  theme_bw(base_size=20)+
  theme(plot.title = element_text(hjust=0.5))+
  labs(title=bquote("PR "*a[2]), x='gene rank')+
  ylim(0,0.3)+
  geom_vline(aes(xintercept= 1424), col='red')


# plot_A1PR/plot_A2PR
# ggsave("./2.make_DS_plot/find_strata/Mv-cat-C212/dS_plot_actual_order_Mv-cat-C212.pdf", device="pdf", height=10, width=15)

plot_A1HD=ggplot(data=d4[which(d4$contig_A1 %in% c('Mv-cat-C212-A1_h2tg000018l')),], aes(x=rank_A1, y=dS))+  
  geom_point(aes(color=as.factor(colors)))+
  scale_colour_identity()+
  theme_bw(base_size=20)+
  theme(plot.title = element_text(hjust=0.5))+
  labs(title=bquote("HD "*a[1]), x='gene rank')+
  ylim(0,0.3)+
  geom_vline(aes(xintercept= 245), col='red')


plot_A2HD=ggplot(data=d4[which(d4$contig_A2 %in% c('Mv-cat-C212-A2_h1tg000040l')),], aes(x=rank_A2, y=dS))+  
  geom_point(aes(color=as.factor(colors)))+
  scale_colour_identity()+
  theme_bw(base_size=20)+
  theme(plot.title = element_text(hjust=0.5))+
  labs(title=bquote("HD "*a[2]), x='gene rank')+
  ylim(0,0.3)+
  geom_vline(aes(xintercept= 18), col='red')


# plot_A1HD/plot_A2HD
# ggsave("./2.make_DS_plot/find_strata/Mv-cat-C212/dS_plot_actual_order_Mv-cat-C212_HD.pdf", device="pdf", height=10, width=15)

(plot_A1HD+plot_A1PR+plot_A2HD+plot_A2PR)+plot_layout(ncol = 2, nrow=2,width=c(0.21,0.79))+plot_annotation(tag_levels = 'A')

ggsave("./2.make_DS_plot/find_strata/Mv-cat-C212/dS_plot_actual_order_Mv-cat-C212_PRHD_panels.pdf", device="pdf", height=10, width=15)
ggsave("./2.make_DS_plot/find_strata/Mv-cat-C212/dS_plot_actual_order_Mv-cat-C212_PRHD_panels.png", device="png", height=10, width=15)



####---------------------------------------------------------------------------------------------
##Make a file with the COORDINATES of the strata

new_tab=NULL
###PARS
PAR=d4[which(d4$strata_name=='PARs'),]
PAR[order(PAR$rank_A1),]
PARA1a=PAR[which(PAR$rank_A1>=128 & PAR$rank_A1<214),]
PARA1b=PAR[which(PAR$rank_A1>=3951),]
new_tab=rbind.data.frame(new_tab,c("PARA1a","A1","Mv-cat-C212-h2tg000005l",min(PARA1a$start_A1),max(PARA1a$end_A1),"grey",min(PARA1a$rank_A1),max(PARA1a$rank_A1)))
new_tab=rbind.data.frame(new_tab,c("PARA1b","A1","Mv-cat-C212-h2tg000005l",min(PARA1b$start_A1),max(PARA1b$end_A1),"grey",min(PARA1b$rank_A1),max(PARA1b$rank_A1)))

PAR[order(PAR$rank_A2),]
PARA2a=PAR[which(PAR$rank_A2>=35 & PAR$rank_A2<119),]
PARA2b=PAR[which(PAR$rank_A2>=3561),]
new_tab=rbind.data.frame(new_tab,c("PARA2a","A2","Mv-cat-C212-h1tg000006l",min(PARA2a$start_A2),max(PARA2a$end_A2),"grey",min(PARA2a$rank_A2),max(PARA2a$rank_A2)))
new_tab=rbind.data.frame(new_tab,c("PARA2b","A2","Mv-cat-C212-h1tg000006l",min(PARA2b$start_A2),max(PARA2b$end_A2),"grey",min(PARA2b$rank_A2),max(PARA2b$rank_A2)))

## amber
amber=d4[which(d4$strata_name=='amber'),]
amber[order(amber$rank_A1),]
amberA1=amber[which(amber$rank_A1>=2 & amber$rank_A1 <=103),]
new_tab=rbind.data.frame(new_tab,c("amberA1","A1","Mv-cat-C212-h2tg000005l",min(amberA1$start_A1),max(amberA1$end_A1),"amber",min(amberA1$rank_A1),max(amberA1$rank_A1)))

amber[order(amber$rank_A2),]
amberA2=amber[which(amber$rank_A2>=934 & amber$rank_A2 <=1022),]
new_tab=rbind.data.frame(new_tab,c("amberA2","A2","Mv-cat-C212-h1tg000006l",min(amberA2$start_A2),max(amberA2$end_A2),"amber",min(amberA2$rank_A2),max(amberA2$rank_A2)))

# ##black PR
PR=d4[which(d4$contig_A2 %in% c('Mv-cat-C212-A2_h1tg000006l')),]

blackPR=PR[which(PR$strata_name=='black'),]
View(blackPR[order(blackPR$rank_A1),])
blackPRA1=blackPR[which(blackPR$rank_A1>=214 & blackPR$rank_A1 <=3924),]
new_tab=rbind.data.frame(new_tab,c("blackPRA1","A1",'Mv-cat-C212-A1_h2tg000005l',min(blackPRA1$start_A1),max(blackPRA1$end_A1),"black_PR",min(blackPRA1$rank_A1),max(blackPRA1$rank_A1)))

View(blackPR[order(blackPR$rank_A2),])
blackPRA2=blackPR[which(blackPR$rank_A2>=125 & blackPR$rank_A2 <=3544),]
new_tab=rbind.data.frame(new_tab,c("blackPRA2","A2",'Mv-cat-C212-A2_h1tg000006l',min(blackPRA2$start_A1),max(blackPRA2$end_A1),"black_PR",min(blackPRA2$rank_A2),max(blackPRA2$rank_A2)))

###BLUE HD
HD=d4[which(d4$contig_A2 %in% c('Mv-cat-C212-A2_h1tg000040l')),]
blueHD=HD[which(HD$strata_name=='blue'),]
blueHD[order(blueHD$rank_A1),]
blueA1=blueHD[which(blueHD$rank_A1>=225 & blueHD$rank_A1 <=259),]
new_tab=rbind.data.frame(new_tab,c("blueA1","A1","Mv-cat-C212-A1_h2tg000018l",min(blueA1$start_A1),max(blueA1$end_A1),"blue",min(blueA1$rank_A1),max(blueA1$rank_A1)))

blueHD[order(blueHD$rank_A2),]
blueA2=blueHD[which(blueHD$rank_A2>=5 & blueHD$rank_A2 <=36),]
new_tab=rbind.data.frame(new_tab,c("blueA2","A2","Mv-cat-C212-A2_h1tg000040l",min(blueA2$start_A2),max(blueA2$end_A2),"blue",min(blueA2$rank_A2),max(blueA2$rank_A2)))


# ##dark ruby
# DR=d4[which(d4$strata_name=='dark_ruby'),]
# DR[order(DR$rank_A1),]
# DRA1=DR[which(DR$rank_A1>=2 & DR$rank_A1 <=103),]
# new_tab=rbind.data.frame(new_tab,c("dark_rubyA1","A1","Mv-cat-C212-h2tg000005l",min(DRA1$start_A1),max(DRA1$end_A1),"amber",min(DRA1$rank_A1),max(DRA1$rank_A1)))
# 
# 


colnames(new_tab)=c('strata','MT','contig','start','end','color','rankmin','rankmax')

write.table(new_tab, file=('./2.make_DS_plot/find_strata/Mv-cat-C212/Position_strata_Mv-cat-C212.txt'),row.names = F, quote = F)

tab_A1=new_tab[which(new_tab$MT=='A1' & new_tab$contig %in% c('Mv-cat-C212-h2tg000005l','Mv-cat-C212-h2tg000005l')),]
plot_A1=ggplot()+
  geom_point(data=d4[which(d4$contig_A1=='Mv-cat-C212-A1_h2tg000005l'),], aes(x=rank_A1, y=dS,color=as.factor(colors)))+
  scale_color_identity()+
  geom_rect(data=tab_A1, aes(xmin=as.numeric(as.character(tab_A1$rankmin)),xmax=as.numeric(as.character(tab_A1$rankmax)),ymin=0,ymax=0.3),alpha=0.2)+
  theme_bw()+
  ylim(0,0.3)
plot_A1

tab_A2=new_tab[which(new_tab$MT=='A2' & new_tab$contig %in% c('Mv-cat-C212-h1tg000006l')),]
plot_A2=ggplot()+  
  geom_point(data=d4[which(d4$contig_A2 %in% c('Mv-cat-C212-A2_h1tg000006l')),], aes(x=rank_A2, y=dS,color=as.factor(colors)))+
  scale_colour_identity()+
  geom_rect(data=tab_A2, aes(xmin=as.numeric(as.character(tab_A2$rankmin)),xmax=as.numeric(as.character(tab_A2$rankmax)),ymin=0,ymax=0.3),alpha=0.2)+
  theme_bw()+
  ylim(0,0.3)
plot_A2
plot_A1/plot_A2

ggsave("./2.make_DS_plot/find_strata/Mv-cat-C212/dS_plot_actual_order_Mv-cat-C212_with_rectangles_PR.pdf", device="pdf", height=10, width=15)


tab_A1=new_tab[which(new_tab$MT=='A1' & new_tab$contig %in% c('Mv-cat-C212-A1_h2tg000018l')),]
plot_A1=ggplot()+
  geom_point(data=d4[which(d4$contig_A1 %in% c('Mv-cat-C212-A1_h2tg000018l')),], aes(x=rank_A1, y=dS,color=as.factor(colors)))+
  scale_color_identity()+
  geom_rect(data=tab_A1, aes(xmin=as.numeric(as.character(tab_A1$rankmin)),xmax=as.numeric(as.character(tab_A1$rankmax)),ymin=0,ymax=0.3),alpha=0.2)+
  theme_bw()+
  ylim(0,0.3)
plot_A1

tab_A2=new_tab[which(new_tab$MT=='A2' & new_tab$contig %in% c('Mv-cat-C212-A2_h1tg000040l')),]
plot_A2=ggplot()+  
  geom_point(data=d4[which(d4$contig_A2 %in% c('Mv-cat-C212-A2_h1tg000040l')),], aes(x=rank_A2, y=dS,color=as.factor(colors)))+
  scale_colour_identity()+
  geom_rect(data=tab_A2, aes(xmin=as.numeric(as.character(tab_A2$rankmin)),xmax=as.numeric(as.character(tab_A2$rankmax)),ymin=0,ymax=0.3),alpha=0.2)+
  theme_bw()+
  ylim(0,0.3)

plot_A1/plot_A2
ggsave("./2.make_DS_plot/find_strata/Mv-cat-C212/dS_plot_actual_order_Mv-cat-C212_with_rectangles_HD.pdf", device="pdf", height=10, width=15)





