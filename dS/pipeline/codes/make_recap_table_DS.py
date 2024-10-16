#Make a big table with all the DS to make a plot in R
import os
import sys
path=sys.argv[1]
species=sys.argv[2]

names=set()
f1=open(f'{path}/paml/file_not_empty.txt','r')
for line in f1:
	line=line.strip()
	names.add(line)

d_id={}
for name in names:
	d_id[name]={}
	f2=open(f'{path}/alignement/correspondance_id/{name}_correspondance_id.txt','r')
	for line in f2:
		line=line.split()
		d_id[name][line[2]]=line[1]

d_gene={}
bed=open('/home/elise/Assembly/ortholog_reconstruction/0.data/GTF_final/Mv-lag-1253-A1_fixed.bed','r')
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

f4=open(f'{path}/Orthologues_Mv-lag-1253-A1_{species}_A1_A2.txt','r')
f4.readline()
d_OG={}
for line in f4:
	line=line.split()
	OG=line[0]
	d_OG[OG]=line[1].split('.')[0]
f4.close()
#print(d_gene)


#print(d_id)
w1=open(f'{path}table_all_DS.txt','w')
header=['ortho_gp','name_lag','chrom_lag','start_lag', 'end_lag','seq1','chrom1','real_id1','seq2','chrom2','real_id2','S','N','t','kappa','omega','dN','dN_SE','dS','dS_se']
print('\t'.join(header),file=w1)
for name in names: 
	f2=open(f'{path}/paml/results/{name}/{name}_DNDStable_Yang_Nielsen.txt','r')
	f2.readline()
	for line in f2:
		line=line.split()
		real_id1=d_id[name][line[0]]
		real_id2=d_id[name][line[1]]
		chrom_1=real_id1.split('_')[1]
		chrom_2=real_id2.split('_')[1]
		name_lag=d_OG[name]
		chrom_lag=d_gene[name_lag][0]
		if chrom_lag in['Mv-lag-1253-A1_MC03','Mv-lag-1253-A1_MC12','Mv-lag-1253-A1_MC16']:
			start_lag=d_gene[name_lag][1]
			end_lag=d_gene[name_lag][2]
			new_line=[name,name_lag, chrom_lag, start_lag, end_lag,line[0],chrom_1,real_id1,line[1],chrom_2,real_id2]+line[2:8]+line[9:11]+[line[12]]
			print(new_line)
			print('\t'.join(new_line), file=w1)
