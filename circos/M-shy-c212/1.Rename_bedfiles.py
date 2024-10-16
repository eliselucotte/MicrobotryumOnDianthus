for MT in ['A1','A2']:
	w1=open(f'Mv-sup-1065-{MT}_contig.bed','w')

	f1=open(f'MvDp-1065-20210906-{MT}_contig.bed','r')
	header=f1.readline()
	print('\t'.join(header.split()), file=w1)

	for line in f1:
		line=line.split()
		print('\t'.join(line[0:4]+[f'Mv-sup-1065-{MT}']+line[5:]), file=w1)
	w1.close()
	f1.close()


	w1=open(f'Mv-lag-1253-{MT}_contig.bed','w')
	f1=open(f'MvSv-1253-{MT}-R1-20151026_contig.bed','r')
	header=f1.readline()
	print('\t'.join(header.split()), file=w1)

	for line in f1:
		line=line.split()
		print('\t'.join(line[0:4]+[f'Mv-lag-1253-{MT}']+line[5:]), file=w1)
	w1.close()
	f1.close()

