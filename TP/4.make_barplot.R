library(ggplot2)
library(tidyr)
library(dplyr)
library(wesanderson)
library(patchwork)
setwd('/Users/eliselucotte/desktop/server1/Assembly/transpecific_polymorphism/')
d=read.table('counts_strata_TP.txt', header=F, as.is=T, sep=' ')
colnames(d)=c('Type_gene','Stratum','Count')

#d=gather(d, key='Type_gene','value'="Count", Total, Full, Part, No,Unresolved)

d=d[-which(d$Stratum=='in_HD'),]

d=mutate(d, type2=ifelse(Type_gene=='species','No','Yes'))


rcols=wes_palette("Darjeeling1")
#rcols=wes_palette("Zissou1")
size=20
plot_count=ggplot(data=d,aes(x=as.factor(Stratum), fill=Type_gene)) +
  geom_bar(data=subset(d,type2=="Yes"),aes(y=Count),stat='identity')+
  geom_bar(data=subset(d,type2=="No"),aes(y=Count*-1),stat='identity')+
  scale_y_continuous(breaks=seq(-200,200,20),labels=abs(seq(-200,200,20))) + 
  scale_fill_manual("Type_gene", breaks = c("species", "partial", "full"), labels=c("Species tree", "Part TP", "Full TP"),
                    values=c(rcols[5], rcols[3],rcols[4]), name='gene type')+
  theme_bw(base_size=size)+
  scale_x_discrete(breaks=c('PARs','light_ruby','dark_ruby','orange','purple','black','turquoise'), labels=c('PARs','light ruby','dark ruby','orange','purple','black','turquoise'),limits=c('PARs','light_ruby','dark_ruby','orange','purple','black','turquoise'))+
  labs(y='Gene count', x="Stratum",)+
  coord_flip()
plot_count
#rcols[2],rcols[5]
ggsave('gene_count_instratum.pdf',plot=plot_count,device='pdf',height=10,width=10)
ggsave('gene_count_instratum.png',plot=plot_count,device='png',height=10,width=10)



dico_tot=c()
for (stratum in unique(d$Stratum))
{
  temp=d[which(d$Stratum==stratum),]
  dico_tot[stratum]=sum(temp$Count)
}

d=mutate(d, perc=as.numeric(Count)/dico_tot[Stratum])

rcols=wes_palette("Darjeeling1")
# rcols=wes_palette("Zissou1")
size=30

plot_perc=ggplot(data=d,aes(x=as.factor(Stratum), fill=Type_gene)) +
  geom_bar(data=subset(d,type2=="Yes"),aes(y=perc*100),stat='identity')+
  geom_bar(data=subset(d,type2=="No"),aes(y=perc*-100),stat='identity')+
  geom_text(aes(x=Stratum, y=-105, label=dico_tot[Stratum]),size=7)+
  scale_y_continuous(breaks=seq(-100,100,20),labels=abs(seq(-100,100,20)))+
  scale_fill_manual("Type_gene", breaks = c("perc_No", "perc_Part", "perc_Full"), labels=c("Species tree", "Part TP", "Full TP"),
                    values=c(rcols[5], rcols[2], rcols[3],rcols[4]),
                    name='Topology')+
  theme_bw(base_size=size)+
  scale_x_discrete(breaks=c('PARs','light_ruby','dark_ruby','orange','purple','black','turquoise'),
                   labels=c('PARs','light ruby','dark ruby','orange','purple','black','turquoise'),
                   limits=c('PARs','light_ruby','dark_ruby','orange','purple','black','turquoise'))+
  labs(y='Percentage of genes', x="Stratum",)+
  coord_flip()
plot_perc
ggsave('gene_perc_instratum.pdf',plot=plot_perc,device='pdf',height=10,width=20)
ggsave('gene_perc_instratum.png',plot=plot_perc,device='png',height=10,width=20)


patchworkplot=plot_count+plot_perc+ plot_layout(nrow=1, ncol=2, guides = 'collect')+plot_annotation(tag_levels = 'A')& theme(plot.tag = element_text(size = 20))
write.table(d, file='table_count_percentage.txt',quote=F, row.names = F)

patchworkplot
ggsave('gene_percandcount_instratum.png',plot=patchworkplot,device='png',height=10,width=20)
