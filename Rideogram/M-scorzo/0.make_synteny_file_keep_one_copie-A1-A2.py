#Make a big table with all the DS to make a plot in R
import os
import sys
species=sys.argv[1]
path_ortholog=f'/home/elise/Assembly/ortholog_reconstruction/Scorzo/0.data/OrthoFinder/Results_Dec18_1/Orthologues/Orthologues_{species}-A1'


d_gene={}
bed=open(f'/home/elise/Assembly/ortholog_reconstruction/Scorzo/0.data/bed/{species}-A1.bed','r')
for line in bed:
	line=line.split()
	gene=line[3].split('.')[0]
	#print(gene)
	chrom=line[0]
	start=line[1]
	end=line[2]
	d_gene[gene]=[chrom,start,end] #/!\ it doesn't necessarily takes the right transcript but it's not important here
bed.close()

#print(d_gene)

d_gene2={}
bed2=open(f'/home/elise/Assembly/ortholog_reconstruction/Scorzo/0.data/bed/{species}-A2.bed','r')
for line in bed2:
	line=line.split()
	#print(line)
	gene=line[3].split('.')[0]
	#print(gene)
	chrom=line[0]
	start=line[1]
	end=line[2]
	d_gene2[gene]=[chrom,start,end] #/!\ it doesn't necessarily takes the right transcript but it's not important here
bed2.close()

#d_gene2['Mv-cat-C212-A2_h1tg000006l_g6063']=['Mv-cat-C212-A2_h1tg000006l_g6063','5370586','5371891']

#print(d_gene['Mv-lag-1253-A1_MC02_g3802'])
#print(d_gene2['Mv-sup-1065-A1_tig00000041_g4167'])
#print(d_gene2)

d={}
f1=open(f'{path_ortholog}/{species}-A1__v__{species}-A2.tsv','r')
f1.readline()
for line in f1:
	line=line.split('\t')
	OG=line[0]
	Mvlag=line[1].strip()
	Mvsup=line[2].strip()
	if len(Mvlag.split(', ')) !=1 or len(Mvsup.split(', ')) !=1:
		continue
	Mvlag=Mvlag.split(' CDS')[0]
	Mvsup=Mvsup.split(' CDS')[0]
	contig_Mvlag=Mvlag.split('_')[1]
	gene_Mvlag=Mvlag.split('.')[0]
	contig_Mvsup=Mvsup.split('_')[1]
	gene_Mvsup=Mvsup.split('.')[0]
	#print(gene_Mvsup)
	clean_gene_Mvlag=gene_Mvlag.split('.')[0]
	#print(gene_Mvsup)
	try:
		d[clean_gene_Mvlag+'_'+OG].append([gene_Mvlag,contig_Mvlag,d_gene[gene_Mvlag][1],d_gene[gene_Mvlag][2],gene_Mvsup,contig_Mvsup,d_gene2[gene_Mvsup][1],d_gene2[gene_Mvsup][2]])
	except KeyError:
		try:
			d[clean_gene_Mvlag+'_'+OG]=[]
			d[clean_gene_Mvlag+'_'+OG].append([gene_Mvlag,contig_Mvlag,d_gene[gene_Mvlag][1],d_gene[gene_Mvlag][2],gene_Mvsup,contig_Mvsup,d_gene2[gene_Mvsup][1],d_gene2[gene_Mvsup][2]])
		except KeyError:
			print('no key', gene_Mvsup)

f1.close()


#print(d_id)
w1=open(f'synteny_{species}-A1_{species}-A2.txt','w')
header=['ortho_gp','name_lag','chrom_lag','start_lag','end_lag','name1','chrom1','start1','end1']
print('\t'.join(header),file=w1)
for OG in d.keys():
	for i in range(0,len(d[OG])):
		#print(d[OG][i])
		new_line=[OG]+d[OG][i]
		#print(new_line)
		if len(new_line)!=9:
			print('problem',new_line)
			continue
		print('\t'.join(new_line), file=w1)