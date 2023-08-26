#anchors a blast alignment to an existing alignment file

#python anchor.py gene_name.fasta

import sys
from Bio import SeqIO

gene_name = sys.argv[1].split('.')[0]

f = open(gene_name + ".fasta.bout", "r")
temp = f.read()[:-1]
f.close()


query_id= temp.split('\t')[0]
subject_id= temp.split('\t')[1]
sstart= int(temp.split('\t')[2]) - 1

query_seq = temp.split('\t')[3]
subject_seq = temp.split('\t')[4]


#remove subject gaps on blast alignment
c=0
for i in range(len(subject_seq)):
    if subject_seq[c] == '-':
        query_seq = query_seq[:c] + query_seq[(c+1):]
        subject_seq = subject_seq[:c] + subject_seq[(c+1):]
        c-=1
    c+=1

#read original subject
for read in SeqIO.parse("../stock_alignments/{0}.fasta".format(gene_name), 'fasta'):
    if read.id == subject_id:
        original_subject = str(read.seq)
        break

#append original
query_seq = '-' * sstart + query_seq + '-' * (len(original_subject) - sstart - len(query_seq))
subject_seq = original_subject[:sstart] + subject_seq + original_subject[(sstart + len(subject_seq)):]

#read aligned original
for read in SeqIO.parse("../stock_alignments/{0}_sfirst.efa".format(gene_name), 'fasta'):
    if read.id == subject_id:
        stock_seq = str(read.seq)
        break

    
#emit
string = ''
c = 0
for i in range(len(stock_seq)):
    if stock_seq[i] == '-':
        string = string + '-'
    else:
        string = string + query_seq[c]
        c+=1
        
        
#append query to alignment
f = open("../stock_alignments/{0}_sfirst.efa".format(gene_name), 'a')
f.write(">" + query_id + "\n")
f.write(string + '\n')
f.close()
