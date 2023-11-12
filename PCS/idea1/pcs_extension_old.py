


import pandas as pd
gg = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/hg38-gorGor6_outputs/pcs_sequences/all/hg38_gorGor6_all_pcs_seqs.tsv", sep='\t')['first.sequence'].tolist()
temp = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/hg38-mm39_outputs/pcs_sequences/all/hg38_mm39_all_pcs_seqs.tsv", sep='\t')
cutoff = 30    ## lower cutoff of Mouse PCS
temp = temp[temp['pcs.width']>=cutoff]
mm = temp['first.sequence'].tolist() 

'''
indices = []
#added_mm_pcs = []
for s in mm:
  for idx, c in enumerate(gg):
    if s in c:
      #added_mm_pcs.append(s)
      indices.append(idx)
'''

indices = [idx for s in mm for idx,c in enumerate(gg)  if s in c]      
      
main = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/hg38-gorGor6_outputs/pcs_sequences/all/hg38_gorGor6_all_pcs_seqs.tsv", sep='\t')

print(main)

#df2 = main.iloc[indices]
#df2['matched.mm.pcs'] = added_mm_pcs
#df2 = df2[~df2.index.duplicated(keep='first')]
#df2.to_csv("/bucket/.mabuya/MillerU/Abrar/hg38-gorGor6_outputs/pcs_sequences/all/hg38_gorGor6_mm39_pcs_{}.tsv".format(cutoff), sep='\t', index=False)

df_main = main.iloc[indices].drop_duplicates()
pcs_filter = df_main['pcs.width'].tolist()

serial_count=[]
for i in range(0, int(max(pcs_filter)+5)):
    serial_count.append(pcs_filter.count(i))
    
lmer = [x for x in range(int(max(pcs_filter)+5))]
# dictionary of lists  
dict = {'PCS length (L)': lmer, 'Number of PCS of Length L': serial_count}  
        
df3 = pd.DataFrame(dict)
df3.to_csv(r'/bucket/MillerU/Abrar/PCS/PCS_hg38-gg6-mm39_{}.csv'.format(cutoff), index = False) 


