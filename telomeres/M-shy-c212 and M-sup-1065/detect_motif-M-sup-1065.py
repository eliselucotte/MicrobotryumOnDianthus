from collections import defaultdict
from Bio import SeqIO
import re
species='MvDp-1065'
contigs={}
contigs['A1']=['tig00000040','tig00000080','tig00000084']
contigs['A2']=['tig00000031','tig00000068']

d_motif={}
d_motif['CCCTAA']='_inv_compl'
d_motif['GGGATT']='_inv'
d_motif['AATCCC']='_compl'
d_motif['TTAGGG']=''
window=1000
for motif in d_motif.keys():
	pattern=re.compile(motif)
	for MAT in ["A1","A2"]:
		dico_match=defaultdict(lambda:defaultdict(int))
		w1=open(f'{species}-{MAT}-count_telomeres_motifs{d_motif[motif]}_{window}bp.txt','w')
		print('\t'.join(['contig','interval','count']), file=w1)

		path=f'/home/elise/Assembly/fasta_files/'
		for record in SeqIO.parse(f'{path}{species}-{MAT}.contigs.fasta', 'fasta'):
			if record.id in contigs[MAT]:
				print(record.id)
				#print(str(record.seq))
				list_match=[]
				for match in pattern.finditer(str(record.seq)):
					list_match.append(match.start())
				#print(list_match)
				for i in range(0,len(record.seq),window):
					for match in list_match:
						if int(match)>=i and int(match)<i+window:
							dico_match[record.id][i]+=1
		for contig in dico_match.keys():
			for interval in dico_match[contig].keys():
				print('\t'.join([contig,str(interval),str(dico_match[contig][interval])]), file=w1)

