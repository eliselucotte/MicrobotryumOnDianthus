import subprocess
import shlex
from collections import defaultdict
import sys
alignment_path=sys.argv[1]
path=sys.argv[2]
ortho=sys.argv[3]

template=open(f'{path}/yn00_template.ctl','r')

infile=f'{ortho}_NT_phylip.fa'
outfile=f'yn00_orthogp.{ortho}.out'
ctl_file=f'yn00_{ortho}.ctl'

w1=open(f'{path}/results/{ortho}/{ctl_file}','w')
i=0
for line in template:
	line=line.strip()
	if i==0:
		new_line=line.replace('PATH/SEQFILE', f'{alignment_path}{infile}')
		print(new_line, file=w1)
	elif i==1:
		new_line=line.replace('PATH/OUTFILE', f'{path}{outfile}')
		print(new_line, file=w1)
	else:
		print(line, file=w1)
	i+=1