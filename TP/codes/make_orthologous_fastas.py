#Create multi sequences fasta files to run alignments on macse per gene
from Bio import SeqIO
import sys
path='/home/elise/Assembly/ortholog_reconstruction/0.data/GTF_final/'
#read the orthologous file
d_ortho={}
f1=open(f'data_strata_datation_Mv-sup-1065.txt','r')
f1.readline()
for line in f1:
	line=line.split()
	print(line)
	ortho=line[0]
	d_ortho[ortho]={}
	d_ortho[ortho]['Mv-lag-1253-A1']=[line[1]]
	d_ortho[ortho]['Mv-sup-1065-A1']=[line[2]]
	d_ortho[ortho]['Mv-sup-1065-A2']=[line[3]]
	d_ortho[ortho]['Mv-cat-C212-A1']=line[4].split(',')
	d_ortho[ortho]['Mv-cat-C212-A2']=line[5].split(',')
	d_ortho[ortho]['Mv-lag-1253-A2']=line[6].split(',')

print(d_ortho)
#creates fasta files with the fasta sequences for each orthologous group
path='/home/elise/Assembly/ortholog_reconstruction/0.data/GTF_final/'
list_ortho=list(d_ortho.keys())
d_set={}
for ortho in list_ortho:
	d_set[ortho]={}
	w1=open(f'./fasta_ortho/{ortho}.fa', 'w')
	for species in d_ortho[ortho].keys():
		d_set[ortho][species]=set()
		#print(species)
		if species in ['Mv-lag-1253-A1','Mv-lag-1253-A2']:
			file=f'{path}{species}.cds.fa'
		else:
			file=f'{path}{species}_final.cds.fa'
		for seq_ortho in SeqIO.parse(file, "fasta"):
			gene=seq_ortho.id
			gene=gene.split('.')[0]
			#print(gene)
			if gene in d_ortho[ortho][species] and gene not in d_set[ortho][species]:
				print('>'+seq_ortho.id, file=w1)
				print(seq_ortho.seq, file=w1)
				d_set[ortho][species].add(gene)
