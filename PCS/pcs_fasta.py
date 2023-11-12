
import os
import pandas as pd
df = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/hg38-galVar1_outputs/pcs_sequences/all/hg38_galVar1_all_pcs_seqs.tsv", sep='\t')
a = df.drop(['second.seqnames','second.start','second.end','second.sequence','second.strand'],axis=1)
b = df.drop(['first.seqnames','first.start','first.end','first.sequence','first.strand'],axis=1)

result1 = []
for i in a.iterrows():
    result1.append('>'+ str(i[1][0]) + ':' + str(i[1][1]) + '-' + str(i[1][2]))#0 means seqnames, 1 means first.start
    result1.append(i[1][3])
    
final_a = pd.DataFrame(result1)
final_a.to_csv(r'/bucket/MillerU/Abrar/hg38-galVar1_outputs/pcs_sequences/all/hg38_pcs.fasta', sep = '\n', index = None,  header = None)

result2 = []
for i in b.iterrows():
    result2.append('>'+ str(i[1][0]) + ':' + str(i[1][1]) + '-' + str(i[1][2]))#0 means seqnames, 1 means second.start
    result2.append(i[1][3])
    
final_b = pd.DataFrame(result2)
final_b.to_csv(r'/bucket/MillerU/Abrar/hg38-galVar1_outputs/pcs_sequences/all/galVar1_pcs.fasta', sep = '\n', index = None,  header = None)