import sys

d={}
name_species="M-scorzo"
path='../0.data/OrthoFinder/Results_Dec18_1/Orthologues/Orthologues_Mv-lag-1253-A1'

list_chrom=['MC03','MC12','MC16']
for MAT in ['A1','A2']:
	f1=open(f'{path}/Mv-lag-1253-A1__v__{name_species}-{MAT}.tsv','r')
	f1.readline()
	for line in f1:
		line=line.split('\t')
		OG=line[0]
		Mvlag=line[1].strip()
		Mvsup=line[2].strip()
		if len(Mvlag.split(',')) != 1 or len(Mvsup.split(','))!=1:
			continue
		contig_Mvlag=Mvlag.split('_')[1]
		gene_Mvlag=Mvlag.split('_')[2]
		contig_Mvsup=Mvsup.split('_')[1]
		gene_Mvsup=Mvsup.split('_')[2]
		clean_gene_Mvlag=gene_Mvlag.split('.')[0]
		if contig_Mvlag not in list_chrom:
			continue
		try:
			d[clean_gene_Mvlag+'_'+OG].extend([Mvsup, contig_Mvsup, gene_Mvsup])
		except KeyError:
			d[clean_gene_Mvlag+'_'+OG]=[Mvlag,contig_Mvlag,gene_Mvlag]
			d[clean_gene_Mvlag+'_'+OG].extend([Mvsup, contig_Mvsup, gene_Mvsup])
	f1.close()
#print(d)

w1=open(f'Orthologues_Mv-lag-1253-A1_{name_species}_A1_A2.txt','w')
header=['gene','species1','contig1','gene1','species2','contig2','gene2','species3','contig3','gene3']
print('\t'.join(header), file=w1)
for gene in d.keys():
	if len(d[gene]) !=9:
		continue
	new_line=[gene]+d[gene]
	print('\t'.join(new_line), file=w1)
w1.close()
