MTchrom_A1=['h2tg000005l','h2tg000018l','tig00000040','tig00000080','tig00000084','MC12','MC16','MC03']
MTchrom_A2=['h1tg000006l','h1tg000040l','tig00000068', 'tig00000031','MC12','MC16','MC03']
#import sys
#path=sys.argv[1]
## get in a dictionnary all the orthologs Mv-sup - Mv-cat for A1
# you end up with a dictionnary d_ortho[each gene for Mvsup]=[list of genes in Mvcat]
d_ortho_A1={}
path='/home/elise/Assembly/ortholog_reconstruction/1.results/Results_Jun12_2/Orthologues/'
f1=open(f'{path}/Orthologues_Mv-sup-1065-A1_final.cds.prot_longest_transcript/Mv-sup-1065-A1_final.cds.prot_longest_transcript__v__Mv-cat-C212-A1_final.cds.prot_longest_transcript.tsv','r')
for line in f1:
	line=line.strip()
	line=line.split('\t')
	genes1=line[1]
	genes2=line[2]
	# print(genes1)
	# print(genes2)
	for gene in genes1.split(', '):
		#print(gene)
		gene=gene.split('.')[0]
		tig=gene.split('_')[1]
		if tig not in MTchrom_A1:
			continue
		d_ortho_A1[gene]=[]
		for gene2 in genes2.split(', '):
			#print(gene2)
			gene2=gene2.split('.')[0]
			tig=gene2.split('_')[1]
			if tig not in MTchrom_A1:
				continue
			else:
				d_ortho_A1[gene].append(gene2)
#print(d_ortho_A1)

## get in a dictionnary all the orthologs Mv-sup - Mv-cat for A2
d_ortho_A2={}
f2=open(f'{path}Orthologues_Mv-sup-1065-A2_final.cds.prot_longest_transcript/Mv-sup-1065-A2_final.cds.prot_longest_transcript__v__Mv-cat-C212-A2_final.cds.prot_longest_transcript.tsv','r')
for line in f2:
	line=line.strip()
	line=line.split('\t')
	genes1=line[1]
	genes2=line[2]
	# you end up with a dictionnary d_ortho[each gene for Mvsup]=[list of genes in Mvcat]
	for gene in genes1.split(', '):
		gene=gene.split('.')[0]
		tig=gene.split('_')[1]
		if tig not in MTchrom_A2:
			continue
		d_ortho_A2[gene]=[]
		for gene2 in genes2.split(', '):
			gene2=gene2.split('.')[0]
			tig=gene2.split('_')[1]
			if tig not in MTchrom_A2:
				continue
			else:
				d_ortho_A2[gene].append(gene2)
#print(d_ortho_A2)


## get in a dictionnary all the orthologs Mv-lag A1 - Mv-lag A2
d_ortho_lag_A1={}
f2=open(f'{path}Orthologues_Mv-lag-1253-A1_final.cds.prot_longest_transcript/Mv-lag-1253-A1_final.cds.prot_longest_transcript__v__Mv-lag-1253-A2_final.cds.prot_longest_transcript.tsv','r')
for line in f2:
	line=line.strip()
	line=line.split('\t')
	genes1=line[1]
	genes2=line[2]
	for gene in genes1.split(', '):
		gene=gene.split('.')[0]
		tig=gene.split('_')[1]
		if tig not in MTchrom_A2:
			continue
		d_ortho_lag_A1[gene]=[]
		for gene2 in genes2.split(', '):
			gene2=gene2.split('.')[0]
			tig=gene2.split('_')[1]
			if tig not in MTchrom_A2:
				continue
			else:
				d_ortho_lag_A1[gene].append(gene2)
#print(d_ortho_lag_A1)

##Read a first time the file, and add a 'NA' if there are no orthologs in MvCat
f0=open('/home/elise/Assembly/ortholog_reconstruction/2.make_DS_plot/Mv-sup-1065/database_for_strata_Mv-sup-1065.txt','r')
f0.readline()
for line in f0:
	line=line.split()
	#print(line)
	gene_lag=line[1]
	gene_A1=line[6]
	gene_A2=line[11]
	try:
		if d_ortho_A1[gene_A1]==[]:
			d_ortho_A1[gene_A1]=['NA']
	except KeyError:
		d_ortho_A1[gene_A1]=['NA']
	try:
		if d_ortho_A2[gene_A2]==[]:
			d_ortho_A2[gene_A2]=['NA']
	except KeyError:
		d_ortho_A2[gene_A2]=['NA']

	try:
		if d_ortho_lag_A1[gene_lag]==[]:
			d_ortho_lag_A1[gene_lag]=['NA']
	except KeyError:
		d_ortho_lag_A1[gene_lag]=['NA']
f0.close()

w1=open('../data_strata_datation_Mv-sup-1065.txt','w')
f0=open('/home/elise/Assembly/ortholog_reconstruction/2.make_DS_plot/Mv-sup-1065/database_for_strata_Mv-sup-1065.txt','r')
header=f0.readline()
new_header=['OG','gene_Mv-lag-1253-A1','gene_Mv-sup-1065-A1','gene_Mv-sup-1065-A2','gene_Mv-cat-C212-A1', 'gene_Mv-cat-C212-A2','gene_Mv-lag-1253-A2','strata']
print('\t'.join(new_header), file=w1)
for line in f0:
	line=line.split()
	#print(line)
	gene_lag=line[1]
	gene_A1=line[6]
	gene_A2=line[11]
	elt=(',').join(d_ortho_A1[gene_A1])
	elt2=(',').join(d_ortho_A2[gene_A2])
	elt3=(',').join(d_ortho_lag_A1[gene_lag])
	new_line=[line[0],line[1], line[6],line[11]]+[elt, elt2, elt3,line[-1]]
	print('\t'.join(new_line), file=w1)





