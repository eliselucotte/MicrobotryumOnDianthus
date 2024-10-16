###PLOT THE DS AND THE STRATA ON THE ACTUAL ORDER 
library(ggplot2)
library(patchwork)

setwd('/Users/eliselucotte/desktop/server1/Assembly/ortholog_reconstruction/')

d3=read.table('./2.make_DS_plot/Mv-sup-1065/database_for_strata_Mv-sup-1065.txt', as.is=T, header=T)

#positions centromeres
# tig00000080	6288000	7152000
# tig00000040	293000	657000
#tig00000080	3309490	3310840	PR
#tig00000040	845924	847447	HD2

#A2
#tig00000031	295000	786000
#tig00000068	5763000	6792000
# tig00000068	10857912	10859232	PR
# tig00000031	983597	984924	HD1

#get the rank for PR
d3[which(d3$contig_A1=='Mv-sup-1065-A1_tig00000080' & d3$start_A1>3200000 & d3$start_A1<3400000),]

d3[which(d3$contig_A2=='Mv-sup-1065-A2_tig00000068' & d3$start_A2>10840000 & d3$start_A2<10880000),]

#get the rank for HD
d3[which(d3$contig_A1=='Mv-sup-1065-A1_tig00000040' & d3$start_A1>840000 & d3$start_A1<850000),]

d3[which(d3$contig_A2=='Mv-sup-1065-A2_tig00000031' & d3$start_A2>983597-10000 & d3$start_A2<983597+10000),]

#get ranks for centromere
d3[which(d3$contig_A1=='Mv-sup-1065-A1_tig00000080' & d3$start_A1>6280000 & d3$start_A1<6290000),]
d3[which(d3$contig_A1=='Mv-sup-1065-A1_tig00000080' & d3$start_A1>7150000 & d3$start_A1<7160000),]

d3[which(d3$contig_A2=='Mv-sup-1065-A2_tig00000068' & d3$start_A2>5750000 & d3$start_A2<5800000),]
d3[which(d3$contig_A2=='Mv-sup-1065-A2_tig00000068' & d3$start_A2>6780000 & d3$start_A2<6799000),]

	
d3[which(d3$contig_A1=='Mv-sup-1065-A1_tig00000040' & d3$start_A1>293000-10000 & d3$start_A1<293000+10000),]
d3[which(d3$contig_A1=='Mv-sup-1065-A1_tig00000040' & d3$start_A1>657000-10000 & d3$start_A1<657000+10000),]

d3[which(d3$contig_A2=='Mv-sup-1065-A2_tig00000031' & d3$start_A2>295000-10000 & d3$start_A2<295000+10000),]
d3[which(d3$contig_A2=='Mv-sup-1065-A2_tig00000031' & d3$start_A2>786000-10000 & d3$start_A2<786000+10000),]

	

plot_A1PR=ggplot(data=d3[which(d3$contig_A1 %in% c('Mv-sup-1065-A1_tig00000084','Mv-sup-1065-A1_tig00000080')),], aes(x=rank_A1, y=dS))+  
  geom_rect(mapping=aes(xmin=1248, xmax=1449, ymin=0, ymax=0.3), alpha=0.01, fill="lightgoldenrod1")+
  geom_point(aes(color=as.factor(colors)))+
  scale_colour_identity()+
  theme_bw(base_size=20)+
  theme(plot.title = element_text(hjust=0.5))+
  labs(title=bquote("PR "*a[1]), x='gene rank')+
  ylim(0,0.3)+
  geom_vline(aes(xintercept= 2234), col='red')


plot_A2PR=ggplot(data=d3[which(d3$contig_A2 %in% c('Mv-sup-1065-A2_tig00000068')),], aes(x=rank_A2, y=dS))+  
  geom_rect(mapping=aes(xmin=1248, xmax=1822, ymin=0, ymax=0.3), alpha=0.01, fill="lightgoldenrod1")+
  geom_point(aes(color=as.factor(colors)))+
  scale_colour_identity()+
  theme_bw(base_size=20)+
  theme(plot.title = element_text(hjust=0.5))+
  labs(title=bquote("PR "*a[2]), x='gene rank')+
  ylim(0,0.3)+
  geom_vline(aes(xintercept= 2952), col='red')


# plot_A1PR/plot_A2PR
# 
# ggsave("./2.make_DS_plot/find_strata/Mv-sup-1065/dS_plot_actual_order_Mv-sup-1065.pdf", device="pdf", height=10, width=15)

plot_A1HD=ggplot(data=d3[which(d3$contig_A1 %in% c('Mv-sup-1065-A1_tig00000040')),], aes(x=rank_A1, y=dS))+  
  geom_rect(mapping=aes(xmin=79, xmax=167, ymin=0, ymax=0.3), alpha=0.01, fill="lightgoldenrod1")+
  geom_point(aes(color=as.factor(colors)))+
  scale_colour_identity()+
  theme_bw(base_size=20)+
  theme(plot.title = element_text(hjust=0.5))+
  labs(title=bquote("HD "*a[1]), x='gene rank')+
  ylim(0,0.3)+
  geom_vline(aes(xintercept= 225), col='red')

plot_A2HD=ggplot(data=d3[which(d3$contig_A2 %in% c('Mv-sup-1065-A2_tig00000031')),], aes(x=rank_A2, y=dS))+  
  geom_rect(mapping=aes(xmin=84, xmax=204, ymin=0, ymax=0.3), alpha=0.01, fill="lightgoldenrod1")+
  geom_point(aes(color=as.factor(colors)))+
  scale_colour_identity()+
  theme_bw(base_size=20)+
  theme(plot.title = element_text(hjust=0.5))+
  labs(title=bquote("HD "*a[2]), x='gene rank')+
  ylim(0,0.3)+
  geom_vline(aes(xintercept= 266), col='red')

# 
# plot_A1HD/plot_A2HD
# ggsave("./2.make_DS_plot/find_strata/Mv-sup-1065/dS_plot_actual_order_Mv-sup-1065_HD.pdf", device="pdf", height=10, width=15)
# 
# plot_A1HD+plot_A1PR+plot_layout(ncol = 2, nrow=1,width=c(0.21,0.79))
# ggsave("./2.make_DS_plot/find_strata/Mv-sup-1065/dS_plot_actual_order_Mv-sup-1065_A1_PRHD.pdf", device="pdf", height=5, width=15)
# 
# plot_A2HD+plot_A2PR+plot_layout(ncol = 2, nrow=1,width=c(0.21,0.79))
# ggsave("./2.make_DS_plot/find_strata/Mv-sup-1065/dS_plot_actual_order_Mv-sup-1065_A2_PRHD.pdf", device="pdf", height=5, width=15)


(plot_A1HD+plot_A1PR+plot_A2HD+plot_A2PR)+plot_layout(ncol = 2, nrow=2,width=c(0.21,0.79))+plot_annotation(tag_levels = 'A')

ggsave("./2.make_DS_plot/find_strata/Mv-sup-1065/dS_plot_actual_order_Mv-sup-1065_A2_PRHD_panels.pdf", device="pdf", height=10, width=15)
ggsave("./2.make_DS_plot/find_strata/Mv-sup-1065/dS_plot_actual_order_Mv-sup-1065_A2_PRHD_panels.png", device="png", height=10, width=15)


 ##Make a file for Paul

new_tab=NULL
###PARS
PAR=d3[which(d3$strata_name=='PARs'),]
PAR[order(PAR$rank_A1),]
#PARA1a=PAR[which(PAR$rank_A1>=10 & PAR$rank_A1<29),]
PARA1b=PAR[which(PAR$rank_A1>=3073),]
#new_tab=rbind.data.frame(new_tab,c("PARA1a","A1","Mv-sup-1065-A1_tig00000084",min(PARA1a$start_A1),max(PARA1a$end_A1),"grey",min(PARA1a$rank_A1),max(PARA1a$rank_A1)))
new_tab=rbind.data.frame(new_tab,c("PARA1b","A1","Mv-sup-1065-A1_tig00000080",min(PARA1b$start_A1),max(PARA1b$end_A1),"grey",min(PARA1b$rank_A1),max(PARA1b$rank_A1)))

PAR[order(PAR$rank_A2),]
#PARA2a=PAR[which(PAR$rank_A2>=11 & PAR$rank_A2<28),]
PARA2b=PAR[which(PAR$rank_A2>=3636),]
#new_tab=rbind.data.frame(new_tab,c("PARA2a","A2","Mv-sup-1065-A2_tig00000068",min(PARA2a$start_A2),max(PARA2a$end_A2),"grey",min(PARA2a$rank_A2),max(PARA2a$rank_A2)))
new_tab=rbind.data.frame(new_tab,c("PARA2b","A2","Mv-sup-1065-A2_tig00000068",min(PARA2b$start_A2),max(PARA2b$end_A2),"grey",min(PARA2b$rank_A2),max(PARA2b$rank_A2)))

##light ruby
##2 genes in the middle are actually PARs
## Mv-sup-1065-A2_tig00000068_g7145 Mv-sup-1065-A2_tig00000068_g7146
light_ruby=d3[which(d3$strata_name=='light_ruby'),]
light_ruby[order(light_ruby$rank_A1),]
temp=d3[which(d3$rank_A1>=15 & d3$rank_A1<=85 & d3$contig_A1=='Mv-sup-1065-A1_tig00000084'),]

temp[order(temp$rank_A1),]
light_rubyA1=light_ruby[which(light_ruby$rank_A1>=29 & light_ruby$rank_A1 <=49),]
new_tab=rbind.data.frame(new_tab,c("light_rubyA1","A1","Mv-sup-1065-A1_tig00000084",min(light_rubyA1$start_A1),max(light_rubyA1$end_A1),"indianred1",min(light_rubyA1$rank_A1),max(light_rubyA1$rank_A1)))

temp=d3[which(d3$rank_A2>=0 & d3$rank_A2<=70 & d3$contig_A2=='Mv-sup-1065-A2_tig00000068'),]
temp[order(temp$rank_A2),]
light_rubyA2=light_ruby[which(light_ruby$rank_A2>=28 & light_ruby$rank_A2 <=55),]
new_tab=rbind.data.frame(new_tab,c("light_rubyA2","A2","Mv-sup-1065-A2_tig00000068",min(light_rubyA2$start_A2),max(light_rubyA2$end_A2),"indianred1",min(light_rubyA2$rank_A2),max(light_rubyA2$rank_A2)))



## dark ruby
dark_ruby=d3[which(d3$strata_name=='dark_ruby'),]
dark_ruby=dark_ruby[order(dark_ruby$rank_A1),]
temp=d3[which(d3$rank_A1>=500 & d3$rank_A1 <=1000),]
temp[order(temp$rank_A1),]
dark_rubyA1a=dark_ruby[which(dark_ruby$rank_A1>=574 & dark_ruby$rank_A1 <=577),]
dark_rubyA1b=dark_ruby[which(dark_ruby$rank_A1>=723 & dark_ruby$rank_A1 <=731),]
dark_rubyA1c=dark_ruby[which(dark_ruby$rank_A1>=784 & dark_ruby$rank_A1 <=811),]
new_tab=rbind.data.frame(new_tab,c("dark_rubyA1a","A1","Mv-sup-1065-A1_tig00000084",min(dark_rubyA1a$start_A1),max(dark_rubyA1a$end_A1),"firebrick",min(dark_rubyA1a$rank_A1),max(dark_rubyA1a$rank_A1)))
new_tab=rbind.data.frame(new_tab,c("dark_rubyA1b","A1","Mv-sup-1065-A1_tig00000084",min(dark_rubyA1b$start_A1),max(dark_rubyA1b$end_A1),"firebrick",min(dark_rubyA1b$rank_A1),max(dark_rubyA1b$rank_A1)))
new_tab=rbind.data.frame(new_tab,c("dark_rubyA1c","A1","Mv-sup-1065-A1_tig00000084",min(dark_rubyA1c$start_A1),max(dark_rubyA1c$end_A1),"firebrick",min(dark_rubyA1c$rank_A1),max(dark_rubyA1c$rank_A1)))


dark_ruby[order(dark_ruby$rank_A2),]
temp=d3[which(d3$rank_A2>=50 & d3$rank_A2 <=100 & d3$contig_A2=='Mv-sup-1065-A2_tig00000068'),]
temp[order(temp$rank_A2),]
dark_rubyA2a=dark_ruby[which(dark_ruby$rank_A2>=57 & dark_ruby$rank_A2 <=71),]
new_tab=rbind.data.frame(new_tab,c("dark_rubyA2a","A2","Mv-sup-1065-A2_tig00000068",min(dark_rubyA2a$start_A2),max(dark_rubyA2a$end_A2),"firebrick",min(dark_rubyA2a$rank_A2),max(dark_rubyA2a$rank_A2)))


## turquoise
turquoise=d3[which(d3$strata_name=='turquoise'),]
turquoise[order(turquoise$rank_A1),]
turquoiseA1a=turquoise[which(turquoise$rank_A1>=191 & turquoise$rank_A1 <=245),]
turquoiseA1b=turquoise[which(turquoise$rank_A1>=733 & turquoise$rank_A1 <=765),]
new_tab=rbind.data.frame(new_tab,c("turquoiseA1a","A1","Mv-sup-1065-A1_tig00000084",min(turquoiseA1a$start_A1),max(turquoiseA1a$end_A1),"turquoise",min(turquoiseA1a$rank_A1),max(turquoiseA1a$rank_A1)))
new_tab=rbind.data.frame(new_tab,c("turquoiseA1b","A1","Mv-sup-1065-A1_tig00000084",min(turquoiseA1b$start_A1),max(turquoiseA1b$end_A1),"turquoise",min(turquoiseA1b$rank_A1),max(turquoiseA1b$rank_A1)))

turquoise[order(turquoise$rank_A2),]
turquoiseA2b=turquoise[which(turquoise$rank_A2>=2898 & turquoise$rank_A2 <=2922),]
turquoiseA2c=turquoise[which(turquoise$rank_A2>=3618 & turquoise$rank_A2 <=3626),]
new_tab=rbind.data.frame(new_tab,c("turquoiseA2b","A2","Mv-sup-1065-A2_tig00000068",min(turquoiseA2b$start_A2),max(turquoiseA2b$end_A2),"turquoise",min(turquoiseA2b$rank_A2),max(turquoiseA2b$rank_A2)))
new_tab=rbind.data.frame(new_tab,c("turquoiseA2c","A2","Mv-sup-1065-A2_tig00000068",min(turquoiseA2c$start_A2),max(turquoiseA2c$end_A2),"turquoise",min(turquoiseA2c$rank_A2),max(turquoiseA2c$rank_A2)))

# ##black PR
# PR=d3[which(d3$contig_A2 %in% c('Mv-sup-1065-A2_tig00000068')),]
# 
# blackPR=PR[which(PR$strata_name=='black'),]
# View(blackPR[order(blackPR$rank_A1),])
# blackPRA1=blackPR[which(blackPR$rank_A1>=829 & blackPR$rank_A1 <=3051),]
# new_tab=rbind.data.frame(new_tab,c("blackPRA1","A1","Mv-sup-1065-A1_tig00000040",min(blackPRA1$start_A1),max(blackPRA1$end_A1)))
# 
# View(blackPR[order(blackPR$rank_A2),])
# blackHDA2=blackPR[which(blackPR$rank_A2>=296 & blackPR$rank_A2 <=348),]
# new_tab=rbind.data.frame(new_tab,c("blackHDA2","A2","Mv-sup-1065-A2_tig00000031",min(blackHDA2$start_A1),max(blackHDA2$end_A1)))
# 
# 


##black HD
HD=d3[which(d3$contig_A2 %in% c('Mv-sup-1065-A2_tig00000031')),]
blackHD=HD[which(HD$strata_name=='black'),]
View(blackHD[order(blackHD$rank_A1),])
blackHDA1=blackHD[which(blackHD$rank_A1>=255 & blackHD$rank_A1 <=303),]
new_tab=rbind.data.frame(new_tab,c("blackHDA1","A1","Mv-sup-1065-A1_tig00000040",min(blackHDA1$start_A1),max(blackHDA1$end_A1),"black",min(blackHDA1$rank_A1),max(blackHDA1$rank_A1)))

View(blackHD[order(blackHD$rank_A2),])
blackHDA2=blackHD[which(blackHD$rank_A2>=296 & blackHD$rank_A2 <=348),]
new_tab=rbind.data.frame(new_tab,c("blackHDA2","A2","Mv-sup-1065-A2_tig00000031",min(blackHDA2$start_A1),max(blackHDA2$end_A1),"black",min(blackHDA2$rank_A2),max(blackHDA2$rank_A2)))


blueHD=HD[which(HD$strata_name=='blue'),]
blueHD[order(blueHD$rank_A1),]
blueA1=blueHD[which(blueHD$rank_A1>=208 & blueHD$rank_A1 <=254),]
new_tab=rbind.data.frame(new_tab,c("blueA1","A1","Mv-sup-1065-A1_tig00000040",min(blueA1$start_A1),max(blueA1$end_A1),"blue",min(blueA1$rank_A1),max(blueA1$rank_A1)))

blueHD[order(blueHD$rank_A2),]
blueA2=blueHD[which(blueHD$rank_A2>=249 & blueHD$rank_A2 <=295),]
new_tab=rbind.data.frame(new_tab,c("blueA2","A2","Mv-sup-1065-A2_tig00000031",min(blueA2$start_A2),max(blueA2$end_A2),"blue",min(blueA2$rank_A2),max(blueA2$rank_A2)))

colnames(new_tab)=c('strata','MT','contig','start','end','color','rankmin','rankmax')

write.table(new_tab, file=('./2.make_DS_plot/find_strata/Mv-sup-1065/Position_strata_Mv-sup-1065.txt'),row.names = F, quote = F)


##Plot the regions with rectangles
tab_A1=new_tab[which(new_tab$MT=='A1' & new_tab$contig %in% c('Mv-sup-1065-A1_tig00000084','Mv-sup-1065-A1_tig00000080')),]
plot_A1=ggplot()+
  geom_point(data=d3[which(d3$contig_A1 %in% c('Mv-sup-1065-A1_tig00000084','Mv-sup-1065-A1_tig00000080')),], aes(x=rank_A1, y=dS,color=as.factor(colors)))+
  scale_color_identity()+
  geom_rect(data=tab_A1, aes(xmin=as.numeric(as.character(tab_A1$rankmin)),xmax=as.numeric(as.character(tab_A1$rankmax)),ymin=0,ymax=0.3),alpha=0.2)+
  theme_bw()+
  ylim(0,0.3)
plot_A1

tab_A2=new_tab[which(new_tab$MT=='A2' & new_tab$contig %in% c('Mv-sup-1065-A2_tig00000068')),]
plot_A2=ggplot()+  
  geom_point(data=d3[which(d3$contig_A2 %in% c('Mv-sup-1065-A2_tig00000068')),], aes(x=rank_A2, y=dS,color=as.factor(colors)))+
  scale_colour_identity()+
  geom_rect(data=tab_A2, aes(xmin=as.numeric(as.character(tab_A2$rankmin)),xmax=as.numeric(as.character(tab_A2$rankmax)),ymin=0,ymax=0.3),alpha=0.2)+
  theme_bw()+
  ylim(0,0.3)

plot_A1/plot_A2

ggsave("./2.make_DS_plot/find_strata/Mv-sup-1065/dS_plot_actual_order_Mv-sup-1065_with_rectangles_PR.pdf", device="pdf", height=10, width=15)


tab_A1=new_tab[which(new_tab$MT=='A1' & new_tab$contig %in% c('Mv-sup-1065-A1_tig00000040')),]
plot_A1=ggplot()+
  geom_point(data=d3[which(d3$contig_A1 %in% c('Mv-sup-1065-A1_tig00000040')),], aes(x=rank_A1, y=dS,color=as.factor(colors)))+
  scale_color_identity()+
  geom_rect(data=tab_A1, aes(xmin=as.numeric(as.character(tab_A1$rankmin)),xmax=as.numeric(as.character(tab_A1$rankmax)),ymin=0,ymax=0.3),alpha=0.2)+
  theme_bw()+
  ylim(0,0.3)
plot_A1

tab_A2=new_tab[which(new_tab$MT=='A2' & new_tab$contig %in% c('Mv-sup-1065-A2_tig00000031')),]
plot_A2=ggplot()+  
  geom_point(data=d3[which(d3$contig_A2 %in% c('Mv-sup-1065-A2_tig00000031')),], aes(x=rank_A2, y=dS,color=as.factor(colors)))+
  scale_colour_identity()+
  geom_rect(data=tab_A2, aes(xmin=as.numeric(as.character(tab_A2$rankmin)),xmax=as.numeric(as.character(tab_A2$rankmax)),ymin=0,ymax=0.3),alpha=0.2)+
  theme_bw()+
  ylim(0,0.3)

plot_A1/plot_A2
ggsave("./2.make_DS_plot/find_strata/Mv-sup-1065/dS_plot_actual_order_Mv-sup-1065_with_rectangles_HD.pdf", device="pdf", height=10, width=15)
