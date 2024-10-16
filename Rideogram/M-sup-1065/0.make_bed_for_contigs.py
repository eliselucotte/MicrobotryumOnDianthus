path_genomes="/home/elise/Assembly/ortholog_reconstruction/0.data/00_export_data_to_freeze/1.Genomes/"
for species in ['Mv-sup-1065-A1','Mv-sup-1065-A2','Mv-lag-1253-A1','Mv-lag-1253-A2']:
	f0=open(f'{path_genomes}{species}.fa.fai','r')
	w1=open(f'{species}_contig.bed', 'w')
	print('\t'.join(['Chr','Start','End','fill',"species",'size','color']), file=w1)
	for line in f0:
		line=line.split()
		contig=line[0]
		new_line=[contig.split('_')[1], "0", line[1], '969696', species,'12','252525']
		print('\t'.join(new_line), file=w1)