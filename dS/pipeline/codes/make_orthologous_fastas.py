#Create multi sequences fasta files to run alignments on macse per gene
from Bio import SeqIO
import sys
name_species=sys.argv[1]
#read the orthologous file
d_ortho={}
f1=open(f'Orthologues_Mv-lag-1253-A1_{name_species}_A1_A2.txt','r')
f1.readline()
for line in f1:
	line=line.split()
	ortho=line[0]
	d_ortho[ortho]={}
	d_ortho[ortho][f'{name_species}-A1']=line[6]
	d_ortho[ortho][f'{name_species}-A2']=line[9]

#print(d_ortho)
#creates fasta files with the fasta sequences for each orthologous group
path='/home/elise/Assembly/ortholog_reconstruction/0.data/GTF_final/'
list_ortho=list(d_ortho.keys())
for ortho in list_ortho:
	w1=open(f'./fasta_ortho/{ortho}.fa', 'w')
	for species in [f'{name_species}-A1',f'{name_species}-A2']:
		for seq_ortho in SeqIO.parse(f'{path}{species}_final.cds.fa', "fasta"):
			gene=seq_ortho.id
			gene=gene.split('_')[2]
			#print(gene)
			if gene in d_ortho[ortho][species]:
				print('>'+seq_ortho.id, file=w1)
				print(seq_ortho.seq, file=w1)
