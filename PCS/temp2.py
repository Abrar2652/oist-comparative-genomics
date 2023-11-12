


import os
import glob
import pandas as pd
import numpy as np
from numpy import array
from scipy.stats import entropy
from scipy.special import rel_entr, kl_div


genome_list = ['gorGor6']
names = ['Gorilla']
window_size = int(50000/1000) #(mb -> kb) [10000000, 5000000, 2500000, 1200000, 15000, 50000, 100000, 30000]
pinf = float('+inf')
ninf = float('-inf')

########## q_i = true distribution/reference ###############
###### (to avoid p=0) putting 0s in case of all infinity ###########
df3 = pd.read_csv("/bucket/MillerU/Abrar/PCS/hg38-gorGor6/PCS_hg38gorGor6.csv")
a = df3['Number of PCS of Length L'].values.tolist()
serial_count1 = a + [0] * (4000 - len(a))
q_i = serial_count1
print(q_i[:10])
total1 = 0
for ele in range(0, len(q_i)):
    total1 = total1 + q_i[ele]
q_i = [m/total1 for m in q_i]
print(q_i[:10])
                
######## p_i = windows ###############

for genome,name in zip(genome_list,names):
        chrom = []
        win = []
        entro = []
        for i in [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X','Y']:
            lines=[]
            path = '/bucket/MillerU/Abrar/PCS/hg38-{0}/relative_entropy/'.format(genome,window_size,i)
            #os.makedirs(path)
            with open("/bucket/MillerU/Abrar/hg38-{1}_outputs/hg38-{1}_sliding_pcs_{2}kb/chr{0}_{2}kb_min_window_sliding_pcs.txt".format(i, genome,window_size)) as file:
                while True:
                    # Get next line from file
                    line = file.readline().rstrip('\n')
                    # if line is empty
                    # end of file is reached
                    if not line:
                        break
                    elif line==("\n"):
                        continue
                    #print(line)
                    lines.append([line])
            lines_frame = pd.DataFrame(lines, columns=['PCS'])
            
            ########### Single window analysis ############
         
            for k in range(len(lines_frame)) :
                pcs=[]
                pcs.append(lines_frame.iloc[k].values.tolist())
                flat_pcs = [item for sublist in pcs for item in sublist]
                final_pcs = [int(k) for k in ','.join(flat_pcs).split(',')]        
                serial_count2 = []
                for j in range(0, 4000):
                    serial_count2.append(final_pcs.count(j))
                p_i = serial_count2
                total2 = 0
                for ele in range(0, len(p_i)):
                    total2 = total2 + p_i[ele]
                p_i = [m/total2 for m in p_i]
                
                chrom.append(i)
                win.append(k+1)
                x = kl_div(p_i, q_i)
                x[x == pinf] = 0.0
                x[x == ninf] = 0.0
                entro.append(sum(x))
                                               
                
rl_entro_df = pd.DataFrame(list(zip(chrom, win, entro)), columns =['chrom', 'window number', 'relative entropy'])
rl_entro_df.to_csv("/flash/MillerU/chr/sliding_relative_entropy_qdown_{1}kb_win_hggg_ref_hggg_idea3.csv".format(genome,window_size), index=False)
















        
        
        
    
    
    
