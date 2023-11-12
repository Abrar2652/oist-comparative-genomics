
import pandas as pd
from functools import reduce

############## BASE GENOME PAIR = hg38/mm39 and others are compared against it ################
cutoff = 20
mm = pd.read_csv('/bucket/.mabuya/MillerU/Abrar/hg38-mm39_outputs/pcs_sequences/all/hg38_mm39_all_pcs_seqs.tsv', sep='\t')
mm = mm[mm['pcs.width']>=cutoff]
print('done')
mic = pd.read_csv('/bucket/.mabuya/MillerU/Abrar/hg38-micMur3_outputs/pcs_sequences/all/hg38_micMur3_all_pcs_seqs.tsv', sep='\t')['first.sequence']
print('done')
tar = pd.read_csv('/bucket/.mabuya/MillerU/Abrar/hg38-tarSyr2_outputs/pcs_sequences/all/hg38_tarSyr2_all_pcs_seqs.tsv', sep='\t')['first.sequence']
print('done')
ory = pd.read_csv('/bucket/.mabuya/MillerU/Abrar/hg38-oryCun2_outputs/pcs_sequences/all/hg38_oryCun2_all_pcs_seqs.tsv', sep='\t')['first.sequence']
print('done')
equ = pd.read_csv('/bucket/.mabuya/MillerU/Abrar/hg38-equCab3_outputs/pcs_sequences/all/hg38_equCab3_all_pcs_seqs.tsv', sep='\t')['first.sequence']
print('done')
fel = pd.read_csv('/bucket/.mabuya/MillerU/Abrar/hg38-felCat9_outputs/pcs_sequences/all/hg38_felCat9_all_pcs_seqs.tsv', sep='\t')['first.sequence']
print('done')
ovi = pd.read_csv('/bucket/.mabuya/MillerU/Abrar/hg38-oviAri4_outputs/pcs_sequences/all/hg38_oviAri4_all_pcs_seqs.tsv', sep='\t')['first.sequence']
print('done')
sus = pd.read_csv('/bucket/.mabuya/MillerU/Abrar/hg38-susScr11_outputs/pcs_sequences/all/hg38_susScr11_all_pcs_seqs.tsv', sep='\t')['first.sequence']
print('done')
tri = pd.read_csv('/bucket/.mabuya/MillerU/Abrar/hg38-triMan1_outputs/pcs_sequences/all/hg38_triMan1_all_pcs_seqs.tsv', sep='\t')['first.sequence']
print('done')
orn = pd.read_csv('/bucket/.mabuya/MillerU/Abrar/hg38-ornAna2_outputs/pcs_sequences/all/hg38_ornAna2_all_pcs_seqs.tsv', sep='\t')['first.sequence']
print('done')
gal = pd.read_csv('/bucket/.mabuya/MillerU/Abrar/hg38-galGal6_outputs/pcs_sequences/all/hg38_galGal6_all_pcs_seqs.tsv', sep='\t')['first.sequence']
print('done')
tha = pd.read_csv('/bucket/.mabuya/MillerU/Abrar/hg38-thaSir1_outputs/pcs_sequences/all/hg38_thaSir1_all_pcs_seqs.tsv', sep='\t')['first.sequence']
print('done')
xen = pd.read_csv('/bucket/.mabuya/MillerU/Abrar/hg38-xenTro10_outputs/pcs_sequences/all/hg38_xenTro10_all_pcs_seqs.tsv', sep='\t')['first.sequence']
print('done')

data_frames = [mm,mic,tar,ory,equ,fel,ovi,sus,tri,orn,gal,tha,xen]# the PCs based on which others will be filtered -> B

df_merged = reduce(lambda  left,right: pd.merge(left,right,on=['first.sequence'],
                                            how='inner'), data_frames)

#df_merged.to_csv('/bucket/.mabuya/MillerU/Abrar/PCS/distant_common_pcs_merged.tsv', sep='\t')
df_merged.to_csv('distant_common_pcs_intersection.tsv', sep='\t')

pcs_filter = df_merged['pcs.width'].tolist()

serial_count=[]
for i in range(0, int(max(pcs_filter)+5)):
    serial_count.append(pcs_filter.count(i))
    
lmer = [x for x in range(int(max(pcs_filter)+5))]
# dictionary of lists  
dict = {'PCS length (L)': lmer, 'Number of PCS of Length L': serial_count}  
        
df3 = pd.DataFrame(dict)
#df3.to_csv("/bucket/.mabuya/MillerU/Abrar/PCS/PCS_common_merged_distant.csv".format(cutoff), index = False) 
df3.to_csv("PCS_common_intersection_distant_L{}.csv".format(cutoff), index = False)









