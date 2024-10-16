
import sys
from Bio import SeqIO
path=sys.argv[1]
#ortholog='MDag00362'
ortholog=sys.argv[2]
d_id={}
i=1
w1=open(f'{path}alignement/{ortholog}_NT_newid.fa','w')
w2=open(f'{path}alignement/correspondance_id/{ortholog}_correspondance_id.txt','w')
#MDag00362_NT_noexcl.fa
#for seq_ortho in SeqIO.parse(f'{path}{ortholog}.fa_edited.fa_NT.fa', 'fasta'):
for seq_ortho in SeqIO.parse(f'{path}alignement/{ortholog}_NT_noexcl.fa', 'fasta'):
	#print(seq_ortho)
	d_id[seq_ortho.id]=i
	seq_ortho.id=i
	#print(seq_ortho.id)
	print('>'+str(seq_ortho.id), file=w1)
	print(seq_ortho.seq, file=w1)
	i+=1

for name in d_id.keys():
	print('\t'.join([ortholog,name, str(d_id[name])]), file=w2)




