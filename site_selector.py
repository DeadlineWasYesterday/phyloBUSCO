import pandas as pd
from Bio import SeqIO

pd.read_csv('single_copy.list', header = None)

pd.read_csv('mi30448.csv')


df = pd.merge(pd.read_csv('single_copy.list', header = None),pd.read_csv('mi30448.csv'), left_on = 0, right_on = 'gene')


selected = df.sort_values('mean')[:1800]['site'].values


t = []
for r in SeqIO.parse("stock_alignments/super.afa", "fasta"):
	ws = ''
	for i in selected:
		ws = ws + r.seq[i]
	
	t.append(('>' + r.id,ws))


pd.DataFrame(t).to_csv('clean_alignment.afa', sep = '\n', index = 0, header = None)
