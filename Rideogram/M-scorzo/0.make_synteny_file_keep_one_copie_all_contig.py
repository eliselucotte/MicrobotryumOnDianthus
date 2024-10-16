#Make a big table with all the DS to make a plot in R
import os
import sys
species=sys.argv[1]
MT=sys.argv[2]
path_ortholog=f'/home/elise/Assembly/ortholog_reconstruction/Scorzo/0.data/OrthoFinder/Results_Dec18_1/Orthologues/Orthologues_Mv-lag-1253-{MT}/'
d_gene={}
bed=open(f'/home/elise/Assembly/ortholog_reconstruction/Scorzo/0.data/bed/Mv-lag-1253-{MT}.bed','r')
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

d_gene['Mv-lag-1253-A1_MC12_PR']=['Mv-lag-1253-A1_MC12','343080','344431']
d_gene['Mv-lag-1253-A2_MC12_PR']=['Mv-lag-1253-A2_MC12','483002','484218']

#print(d_gene)

d_gene2={}
bed2=open(f'/home/elise/Assembly/ortholog_reconstruction/Scorzo/0.data/bed/M-scorzo-{MT}.bed','r')
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
d_gene2['MscorzoA2_tig00000008_MP006917']=['MscorzoA2_tig00000008','251815','253276']

print(d_gene2['MscorzoA2_tig00000008_MP006917'])
#print(d_gene2['Mv-sup-1065-A1_tig00000041_g4167'])
#print(d_gene2)

d={}
f1=open(f'{path_ortholog}/Mv-lag-1253-{MT}__v__{species}.tsv','r')
f1.readline()
for line in f1:
	line=line.split('\t')
	OG=line[0]
	Mvlag=line[1].strip()
	Mvsup=line[2].strip()
	if len(Mvlag.split(', ')) !=1 or len(Mvsup.split(', ')) !=1:
		continue
	contig_Mvlag=Mvlag.split('_')[1]
	gene_Mvlag=Mvlag.split('.')[0]
	contig_Mvsup=Mvsup.split('_')[1]
	gene_Mvsup=Mvsup.split('.')[0]
	#print(gene_Mvsup)
	clean_gene_Mvlag=gene_Mvlag.split('.')[0]
	#print(gene_Mvsup)
	#print(d_gene2[gene_Mvsup])
	try:
		d[clean_gene_Mvlag+'_'+OG].append([gene_Mvlag,contig_Mvlag,d_gene[gene_Mvlag][1],d_gene[gene_Mvlag][2],gene_Mvsup,contig_Mvsup,d_gene2[gene_Mvsup][1],d_gene2[gene_Mvsup][2]])
	except KeyError:
		try:
			d[clean_gene_Mvlag+'_'+OG]=[]
			d[clean_gene_Mvlag+'_'+OG].append([gene_Mvlag,contig_Mvlag,d_gene[gene_Mvlag][1],d_gene[gene_Mvlag][2],gene_Mvsup,contig_Mvsup,d_gene2[gene_Mvsup][1],d_gene2[gene_Mvsup][2]])
		except KeyError:
			print('no key', gene_Mvsup, gene_Mvlag)

f1.close()


#print(d_id)
w1=open(f'synteny_Mv-lag-1253-{MT}_{species}.txt','w')
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