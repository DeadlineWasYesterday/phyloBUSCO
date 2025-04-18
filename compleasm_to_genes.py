import sys,os,gzip,csv,urllib
import numpy as np
import pandas as pd
from Bio import SeqIO

lin = sys.argv[1]
ualnsites = int(sys.argv[2])
wd = os.getcwd()
max_gene_len = int(sys.argv[3])

print("Lineage: ", lin)
print("Number of sites to align: ", ualnsites)
print("Max gene length: ", max_gene_len)
print("Working directory: ", wd)

if not os.path.isdir('unaligned'):
    os.mkdir(wd+'/unaligned')



def export_gene(df,gene):
    td = df[df['Gene'] == gene][[1,'s']]
    td[1] = '>' + td[1]
    td.to_csv(wd+'/unaligned/{0}.fa'.format(gene),header = None, sep ='\n', index = 0)


def select_from_duplicates(di):
    df = di.drop(columns = ['Codons','s'])
    df = df.sort_values(['Sequence', 'Gene Start'])
    #consider both forward and reverse
    df['d'] = 'f'
    dfr = df.copy()[::-1]
    dfr['d'] = 'r'
    df2 = pd.concat([df,dfr], axis = 0)
    df2['Sequence'] = df2['d']+df2['Sequence']
    df = df2.copy()
    df = df.reset_index(drop = 1)

    df['gr'] = df['Gene']+'_'+df['Gene'].shift(-1)
    df['sr'] = df['Sequence'].shift(-1)

    #gene block count at termination
    cnt = df['gr'].value_counts().reset_index()
    cnt = cnt.rename(columns = {'count':'mc'})
    tmp = df.index #index fuckery
    df = pd.merge(df,cnt, how = 'left')
    df.index = tmp

    #remove end of contig
    df.loc[(df['Sequence'] != df['sr']),'tcr'] = 0

    #remove singleton gene blocks
    df.loc[(df['Sequence'] == df['sr']) & (df['tcr'] == 1), 'tcr'] = 0
    
    #join orientations
    dff = df[df['d'] == 'f']
    dfr = df[df['d'] == 'r']
    di['fc'] = dff['tcr'].values
    di['rc'] = dfr['tcr'].values
    di['tc'] = di['fc'] + di['rc']
    
    #fraction * identity
    di['fi'] = di['Fraction'].astype(float) * di['Identity'].astype(float)
    
    return di.sort_values([1,'Gene','tc','fi'],ascending = [1,1,0,0]).drop_duplicates([1,'Gene'])




#read stock taxa
a1 = pd.read_csv('StockTaxa.txt',sep = '\t', header = None)
#make sure no duplicates
assert sum(a1[0].duplicated()) == 0
assert sum(a1[1].duplicated()) == 0
#replace spaces with underscores
a1[1] = a1[1].apply(lambda x: x.replace(' ','_'))


#download metadata and alignments
fname = '{0}.ssv.gz'.format(sys.argv[1][:2])
if not os.path.isfile(fname):
    print("Downloading stock alignment file.")
    print("https://ava.genome.arizona.edu/~mdalam/leg/{0}".format(fname))
    urllib.request.urlretrieve("https://ava.genome.arizona.edu/~mdalam/leg/{0}".format(fname),fname)

fname = 'j1a.tsv'.format(sys.argv[1][:2])
if not os.path.isfile(fname):
    print("Downloading stock metadata.")
    urllib.request.urlretrieve("https://ava.genome.arizona.edu/~mdalam/leg/j1a.tsv".format(fname),fname)


#find assembly identifiers
a2 = pd.read_csv('j1a.tsv', sep='\t')
assemblies = pd.merge(a1,a2,left_on=0,right_on='on',how = 'left')['Assembly'].values
bassemblies = [bytes(s,'utf-8') for s in assemblies]



#extract stock BUSCO data
t = []
bcomma = bytes(',','utf-8')
with gzip.open('vi.ssv.gz','r') as fin:
    header_line = [s.decode('utf-8') for s in fin.readline().strip('\n'.encode('utf-8')).split(bcomma)]
    for line in fin:
        if line.split(bcomma)[14] in bassemblies:
            t.append([s.decode('utf-8') for s in line.split(bcomma)])

d1 = pd.DataFrame(t,columns = header_line).iloc[:,1:] #there's a stupid index column in the stock
d1 = pd.merge(d1,a2[['on','Assembly']], how = 'left')
d1 = pd.merge(d1,a1,left_on='on',right_on=0,how='left')



#load user BUSCO data
a3 = pd.read_csv('UserTaxaCm.txt', header = None, sep = '\t')
#make sure no duplicates
assert sum(a3[0].duplicated()) == 0
assert sum(a3[1].duplicated()) == 0
#replace spaces with underscores
a3[1] = a3[1].apply(lambda x: x.replace(' ','_'))



d2 = pd.DataFrame()
for p,i in a3.values:
    td = pd.read_csv(p+'/'+lin+'_odb10/full_table.tsv', sep = '\t')
    td[1] = i
    d2 = pd.concat([d2,td], axis = 0)

#only single and duplicated
d2 = d2[(d2['Status'] == 'Single') | (d2['Status'] == 'Duplicated')].reset_index(drop = 1)
#primary key to merge sequences
d2['GG'] = d2[1] + '__' + d2['Best gene'] + '|' + d2['Sequence'] + ':' + d2['Gene Start'].astype(int).astype(str) + '-'+ d2['Gene End'].astype(int).astype(str)
#remove alternate codon configurations
d2 = d2.drop_duplicates('GG')

#load sequences
t = []
for p,i in a3.values:
    td = d2[d2[1] == i]
    vals = td['GG'].values
    for r in SeqIO.parse(p+'/'+lin+'_odb10/translated_protein.fasta', 'fasta'):
        if i + '__' + r.id in vals:
            t.append((i + '__' + r.id, str(r.seq)))

dfa = pd.DataFrame(t, columns = ['GG', 's']).drop_duplicates('GG')
d2 = pd.merge(d2,dfa)

#join stock and user
d2 = pd.concat([d1,d2],axis=0).reset_index(drop = 1)
#in case there are duplicate chromosome names
d2['Sequence'] = d2[1]+'_'+d2['Sequence']

#select best duplicate
d2 = select_from_duplicates(d2)

#filter out genes in < 80% assemblies
d2 = d2[d2['Gene'].isin(d2['Gene'].value_counts()[d2['Gene'].value_counts() > 0.8 * d2[1].nunique()].index)]

#select best genes. length x tc
d2['Length'] = d2['Length'].astype(float)
d3 = d2.groupby('Gene')[['Length','tc']].mean()
#filter out genes too long to align on smol computers
d3 = d3[d3['Length'] < max_gene_len]
d3['ltc'] = d3['Length']*d3['tc']
d3 = d3.sort_values('ltc',ascending = 0)
d3['LCum'] = d3['Length'].cumsum()
genes = d3[d3['LCum'] < ualnsites].index


#export unaligned
d2['s'] = d2['s'].apply(lambda x: x.replace('\n',''))
for gene in genes:
    export_gene(d2,gene)


#for gene in d2['Gene'].unique():
#    export_gene(d2,gene)


