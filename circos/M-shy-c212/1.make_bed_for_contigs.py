d_species={}
d_species['MvDcarth-C212-A1']='/home/elise/Assembly/fasta_files/MvDcarth-C212-A1.fasta.fai'
d_species['MvDcarth-C212-A2']='/home/elise/Assembly/fasta_files/MvDcarth-C212-A2.fasta.fai'
for MT in ['A1','A2']:
	old_species='MvDcarth-C212-'+MT
	f0=open(f'{d_species[old_species]}','r')
	new_species='Mv-cat-C212-'+MT
	w1=open(f'{new_species}_contig.bed', 'w')
	print('\t'.join(['Chr','Start','End','fill',"species",'size','color']), file=w1)
	for line in f0:
		line=line.split()
		contig=line[0]
		new_line=[contig, "0", line[1], '969696', new_species,'12','252525']
		print('\t'.join(new_line), file=w1)