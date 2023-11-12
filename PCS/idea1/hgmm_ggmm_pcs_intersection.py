
import pandas as pd
from functools import reduce

############## BASE GENOME PAIR = hg38/gorGor6 and others are compared against it ################
cutoff = 20
hg_mm = pd.read_csv('/bucket/.mabuya/MillerU/Abrar/hg19-mm10_outputs/pcs_sequences/all/hg19_mm10_all_pcs_seqs.tsv', sep='\t')
hg_mm = hg_mm[hg_mm['pcs.width']>=cutoff]
print('done')
gg_mm = pd.read_csv('/bucket/.mabuya/MillerU/Abrar/gorGor6-mm10_outputs/pcs_sequences/all/gorGor6_mm10_all_pcs_seqs.tsv', sep='\t')['first.sequence']
print('done')

data_frames = [hg_mm, gg_mm]# the PCs based on which others will be filtered -> B

df_merged = reduce(lambda  left,right: pd.merge(left,right,on=['first.sequence'],
                                            how='inner'), data_frames)

df_merged.to_csv('/bucket/.mabuya/MillerU/Abrar/PCS/idea1/common_pcs_intersection_hgmm_ggmm.tsv', sep='\t')

pcs_filter = df_merged['pcs.width'].tolist()

serial_count=[]
for i in range(0, int(max(pcs_filter)+5)):
    serial_count.append(pcs_filter.count(i))
    
lmer = [x for x in range(int(max(pcs_filter)+5))]
# dictionary of lists  
dict = {'PCS length (L)': lmer, 'Number of PCS of Length L': serial_count}  
        
df3 = pd.DataFrame(dict)

df3.to_csv("/bucket/.mabuya/MillerU/Abrar/PCS/idea1/common_intersection_L{}_hgmm_ggmm.csv".format(cutoff), index = False) 




