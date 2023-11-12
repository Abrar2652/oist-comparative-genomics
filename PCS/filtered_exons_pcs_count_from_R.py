
import pandas as pd
import os

genome = 'mm10'

df = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/hg19-{0}_outputs/hg19-{0}_filtered_pcs/hg19_1_all_1_identical_seqs.csv".format(genome))

seq_count = df['width'].tolist()

print(seq_count[:10])


serial_count=[]
for i in range(0, int(max(seq_count)+5)):
    serial_count.append(seq_count.count(i))


lmer = [x for x in range(int(max(seq_count)+5))]
# dictionary of lists  
dict = {'PCS length (L)': lmer, 'Number of PCS of Length L': serial_count}  
       
df = pd.DataFrame(dict) 
    
# saving the dataframe 
df.to_csv(r'/bucket/MillerU/Abrar/PCS/hg19-{0}/Filtered_PCS_hg19{0}.csv'.format(genome), index = False) 



