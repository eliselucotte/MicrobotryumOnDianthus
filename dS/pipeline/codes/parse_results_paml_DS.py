##Parse the output from yn00 of PAML to put them into separate files
import os
import sys
path=sys.argv[1]+'/paml/'
list_files=os.listdir(path+"results/")
#print(list_files)
w2=open(f'{path}file_NA.txt','w')
w4=open(f'{path}file_empty.txt','w')
w3=open(f'{path}file_not_empty.txt','w')
for name in list_files:
	#print(name)
	table=[]
	'/home/elise/Assembly/ortholog_reconstruction/2.make_DS_plot/Mv-cat-C212//paml/results/file_NA.txt/yn00_orthogp.file_NA.txt.out'

	f1=open(f'{path}results/{name}/yn00_orthogp.{name}.out', 'r')
	w1=open(f'{path}results/{name}/{name}_DNDStable_Yang_Nielsen.txt', 'w')
	file=f1.readlines()
	i=0
	flag1='NA'
	flag2='NA'
	for line in file:
		line=line.split()
		if len(line)==0:
			i+=1
			continue
		if line[0]=='seq.':
			flag1=i
			i+=1
		elif line[0]=='(C)':
			flag2=i
			i+=1
		else:
			i+=1
			continue
	#print(flag1, flag2)
	if flag1=='NA' or flag2=='NA':
		print(name, file=w2)
	elif len(file[flag1:flag2])==4:
		print(name, file=w4)

	else:
		print(name, file=w3)
		for line in file[flag1:flag2]:
			line=line.split()
			if len(line)==0:
				continue
			print('\t'.join(line), file=w1)

