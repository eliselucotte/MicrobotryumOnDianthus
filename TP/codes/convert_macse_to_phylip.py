import sys
from Bio import SeqIO
path=sys.argv[1]
ortholog=sys.argv[2]
records=SeqIO.parse(f'{path}/alignement/{ortholog}_NT_newid.fa', 'fasta')
SeqIO.write(records, f'{path}/alignement/{ortholog}_NT_phylip.fa', "phylip")