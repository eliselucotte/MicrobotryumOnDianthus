path='/home/elise/Assembly/ortholog_reconstruction/1.results/Results_Jun12_2/Orthogroups/'

list_species=['Mv-dio-1303-A1','Mv-dio-1303-A2','Mv-lag-1253-A1','Mv-lag-1253-A2','Mv-lyc-1064-A1','Mv-lyc-1064-A2']
wanted_species=set(list_species)


#get a set of the single copy orthologues
SC_OG=set()
f1=open(f'{path}Orthogroups_SingleCopyOrthologues.txt','r')
for line in f1:
	line=line.strip()
	SC_OG.add(line)
f1.close()

w1=open('Orthogroups_filtered.txt','w')

f0=open(f'{path}Orthogroups.txt','r')
for line in f0:
	set_species=set()
	line=line.strip()
	line=line.split(' ')
	OG=line[0].split(':')[0]
	#print(len(line))
	if OG not in SC_OG:
		continue
	for i in range(1,len(line)):
		print(line[i])
		if line[i][0:4]=='CDS':
			continue
		species=line[i].split('_')[0]
		set_species.add(species)
	if set_species.issuperset(wanted_species):
		print(' '.join(line), file=w1)

f0.close()