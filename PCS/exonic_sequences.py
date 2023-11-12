
import pandas as pd 

A = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/hg19-mm10_outputs/pcs_sequences/all/hg19_mm10_all_pcs_seqs_222222.tsv", sep='\t')
A = A.drop(A.columns[[4,5,6,7,8,9,10]],axis = 1)
A = A.sort_values('first_start').reset_index(drop=True)
print(A.head())
B = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/hg19-mm10_outputs/hg19-mm10_filtered_pcs/hg19_1_all_1_identical_seqs.csv")
B = B.drop(B.columns[[3,4]],axis = 1)
B = B.sort_values('first_start').reset_index(drop=True)
print(B.head())
out = pd.merge_asof(B, A, on=['first_end'], direction='forward', suffixes=('','_y'))

out.loc[:,['first_start','first_end']] = out.loc[:,['first_start','first_end']].sub(out.first_start_y, axis=0)

out = out.assign(sequences=out.apply(lambda row: 
          row['first_sequence'][row['first_start']:row['first_end']+1], 
          axis=1)).drop(['first_sequence','first_start_y','first_seqnames_y'], axis=1)

out.update(B)
print(out)
out.to_csv("/bucket/.mabuya/MillerU/Abrar/hg19-mm10_outputs/hg19-mm10_filtered_pcs/exonic_sequences_hg19-mm10.tsv", sep='\t', index = False)