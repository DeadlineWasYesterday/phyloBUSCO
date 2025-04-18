import sys,os,gzip,csv
import numpy as np
import pandas as pd
from Bio import SeqIO



def pf(f):
    f['gaps'] = f.apply(lambda r: sum(r == '-'), axis=1)
    f = f[f['gaps'] < (f.shape[1] - 2)*0.20] #####
    f['unc'] = f.iloc[:,0:f.shape[1]-2].apply(lambda r: r.loc[r != '-'].nunique(), axis = 1)
    return f

def es(adir): #export sites
    t = []
    for i in os.listdir(adir):
        for r in SeqIO.parse(adir+i, 'fasta'):
            t.append((i, r.id, str(r.seq)))
            
    d1 = pd.DataFrame(t)
    #d1 = d1[d1[0].isin(d1[0].value_counts()[d1[0].value_counts() > d1[1].nunique() * 0.9].index)] #90% gene incidence 
    #d1 = d1[d1[1].isin(d1[1].value_counts()[d1[1].value_counts() > d1[0].nunique() * 0.9].index)] #90% assembly incidence
    d1 = d1.drop_duplicates([1,2]) ##
    d1 = d1.pivot(index=0,columns=1,values=2).reset_index() #long to wide
    #fill nulls
    for i,r in d1.iterrows():
        gl = len(r.dropna().iloc[1])
        d1.loc[i,:] = d1.loc[i,:].fillna('-'*gl)
    
    #gene index
    t2 = []
    for i,r in d1.iterrows():
        t2 = t2 + [i]*len(r.iloc[1]) #gene id
    
    #seq to char
    t = []
    for c in d1.columns[1:]:
        t1 = []
        for g in d1.loc[:,c]:
            for char in g:
                t1.append(char) 
        t.append(t1)
        
    ti = list(map(list, zip(*t))) #transpose list
    
    f1 = pd.DataFrame(ti)
    f1['gi'] = t2
    
    return d1,f1


def ef(d,f,ss,sl,p):
    t = f[(f['unc'] >= ss) & (f['unc'] <= sl)].iloc[:,:-3].apply(lambda c: ''.join(list(c))).reset_index()
    t[1] = t[0].apply(lambda x: x.count('-'))
    t = t[t[1] <= 0.25*len(t[0].iloc[0])]
    t['index'] = t['index'].astype(int)
    t2 = pd.DataFrame(d.columns[1:]).rename(columns = {1:'Tx'}).reset_index()
    t = pd.merge(t,t2)
    t['Tx'] = '>'+t['Tx']
    t[['Tx',0]].to_csv(p+'{0}_{1}.afa'.format(ss,sl), sep = '\n', index = 0, header = None)







d,f = es(sys.argv[1]+'/')
f = pf(f)


ef(d,f,int(sys.argv[2]),int(sys.argv[3]),sys.argv[4])



