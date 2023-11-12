
import pandas as pd
from functools import reduce

############## BASE GENOME PAIR = hg38/gorGor6 and others are compared against it ################
cutoff = 20

gg = pd.read_csv('/bucket/.mabuya/MillerU/Abrar/hg38-gorGor6_outputs/pcs_sequences/all/hg38_gorGor6_all_pcs_seqs.tsv', sep='\t')
gg = gg[gg['pcs.width']>=cutoff]
print('done')
pon = pd.read_csv('/bucket/.mabuya/MillerU/Abrar/hg38-ponAbe3_outputs/pcs_sequences/all/hg38_ponAbe3_all_pcs_seqs.tsv', sep='\t')['first.sequence']
print('done')
nom = pd.read_csv('/bucket/.mabuya/MillerU/Abrar/hg38-nomLeu3_outputs/pcs_sequences/all/hg38_nomLeu3_all_pcs_seqs.tsv', sep='\t')['first.sequence']
print('done')
rhe = pd.read_csv('/bucket/.mabuya/MillerU/Abrar/hg38-rheMac10_outputs/pcs_sequences/all/hg38_rheMac10_all_pcs_seqs.tsv', sep='\t')['first.sequence']
print('done')
panpan = pd.read_csv('/bucket/.mabuya/MillerU/Abrar/hg38-panPan3_outputs/pcs_sequences/all/hg38_panPan3_all_pcs_seqs.tsv', sep='\t')['first.sequence']
print('done')
pantro = pd.read_csv('/bucket/.mabuya/MillerU/Abrar/hg38-panTro6_outputs/pcs_sequences/all/hg38_panTro6_all_pcs_seqs.tsv', sep='\t')['first.sequence']
print('done')
papanu = pd.read_csv('/bucket/.mabuya/MillerU/Abrar/hg38-papAnu4_outputs/pcs_sequences/all/hg38_papAnu4_all_pcs_seqs.tsv', sep='\t')['first.sequence']
print('done')
mac = pd.read_csv('/bucket/.mabuya/MillerU/Abrar/hg38-macFas5_outputs/pcs_sequences/all/hg38_macFas5_all_pcs_seqs.tsv', sep='\t')['first.sequence']
print('done')
chl = pd.read_csv('/bucket/.mabuya/MillerU/Abrar/hg38-chlSab2_outputs/pcs_sequences/all/hg38_chlSab2_all_pcs_seqs.tsv', sep='\t')['first.sequence']
print('done')
nas = pd.read_csv('/bucket/.mabuya/MillerU/Abrar/hg38-nasLar1_outputs/pcs_sequences/all/hg38_nasLar1_all_pcs_seqs.tsv', sep='\t')['first.sequence']
print('done')


data_frames = [gg,pon,nom,rhe,panpan,pantro,papanu,mac,chl,nas]# the PCs based on which others will be filtered -> B

df_merged = reduce(lambda  left,right: pd.merge(left,right,on=['first.sequence'],
                                            how='inner'), data_frames)

df_merged.to_csv('/bucket/.mabuya/MillerU/Abrar/PCS/primate_common_pcs_intersection_10_genomes.tsv', sep='\t')

pcs_filter = df_merged['pcs.width'].tolist()

serial_count=[]
for i in range(0, int(max(pcs_filter)+5)):
    serial_count.append(pcs_filter.count(i))
    
lmer = [x for x in range(int(max(pcs_filter)+5))]
# dictionary of lists  
dict = {'PCS length (L)': lmer, 'Number of PCS of Length L': serial_count}  
        
df3 = pd.DataFrame(dict)

df3.to_csv("/bucket/.mabuya/MillerU/Abrar/PCS/PCS_common_intersection_primate_L{}_10_genomes.csv".format(cutoff), index = False) 




