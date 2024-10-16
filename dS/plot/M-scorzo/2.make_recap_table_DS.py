#Make a big table with all the DS to make a plot in R
import os
import sys
d_id={}
species="M-scorzo"

d_gene={}
bed=open('../0.data/bed/Mv-lag-1253-A1.bed','r')
for line in bed:
	line=line.split()
	if line[7]=='gene':
		gene=line[3]
		#print(gene)
		chrom=line[0]
		start=line[1]
		end=line[2]
		d_gene[gene]=[chrom,start,end] #/!\ it doesn't necessarily takes the right transcript but it's not important here
bed.close()

d_gene['Mv-lag-1253-A1_C32_PR']=['Mv-lag-1253-A1_C32','343080','344431']
d_gene['Mv-lag-1253-A1_MC12_PR']=['Mv-lag-1253-A1_MC12_PR','483002','484218']

f4=open(f'Orthologues_Mv-lag-1253-A1_{species}_A1_A2.txt','r')
f4.readline()
d_OG={}
for line in f4:
	line=line.split()
	d_OG[line[4].split('.')[0]]=line[1].split('.')[0]
f4.close()
#print(d_gene)


#print(d_id)
w1=open(f'table_all_DS.txt','w')
header=['name_lag','chrom_lag','start_lag', 'end_lag','seq1','chrom1','real_id1','seq2','chrom2','real_id2','dS','dS_se']
print('\t'.join(header),file=w1)
#Ds SE_Ds Dn SE_Dn Gene1 Gene2
f2=open('results_YN.txt','r')
for line in f2:
	line=line.split()
	geneA1=line[5].split('.')[0]
	geneA2=line[4].split('.')[0]
	chrom_1=geneA1.split('_')[1]
	chrom_2=geneA2.split('_')[1]
	try:
		name_lag=d_OG[geneA1]
	except KeyError:
		continue
	chrom_lag=d_gene[name_lag][0]
	dS=line[0]
	dS_se=line[1]
	if chrom_lag in['Mv-lag-1253-A1_MC03','Mv-lag-1253-A1_MC12','Mv-lag-1253-A1_MC16']:
		start_lag=d_gene[name_lag][1]
		end_lag=d_gene[name_lag][2]
		new_line=[name_lag, chrom_lag, start_lag, end_lag,line[0],chrom_1,geneA1,line[1],chrom_2,geneA2, dS, dS_se]
		print(new_line)
		print('\t'.join(new_line), file=w1)
