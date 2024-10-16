f0=open('/home/elise/Assembly/ortholog_reconstruction/2.make_DS_plot/find_strata/datation_strata/data_strata_datation_Mv-sup-1065.txt','r')
f0.readline()
d_strata={}
for line in f0:
	line=line.split()
	OG=line[0]
	strata=line[-3]
	try:
		d_strata[strata].add(OG)
	except KeyError:
		d_strata[strata]=set()
		d_strata[strata].add(OG)
f0.close()

w1=open('launch_move_trees.sh','w')
print(f'mkdir ./trees/', file=w1)
for strata in d_strata.keys():
	print(f'mkdir ./trees/'+strata, file=w1)
	for OG in d_strata[strata]:
		print(f'cp ./alignement/'+OG+'_NT_phylip.fa.iqtree ./trees/'+strata+'/'+OG+'_NT_phylip.fa.iqtree', file=w1)
w1.close()