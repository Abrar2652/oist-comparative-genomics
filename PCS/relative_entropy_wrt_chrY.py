


import os
import glob
import pandas as pd
import numpy as np
from numpy import array
from scipy.stats import entropy
from scipy.special import rel_entr, kl_div


genome_list = ['gorGor6']
names = ['Gorilla']
window_size = int(10000000/1000) #(mb -> kb) [10000000, 5000000, 2500000, 1200000, 15000, 50000, 100000, 30000]

########## p_i ###############
for genome,name in zip(genome_list,names):
        for i in ['Y']:
            lines=[]
            with open("/bucket/MillerU/Abrar/hg38-{1}_outputs/hg38-{1}_tiled_pcs_{2}kb/chr{0}_{2}kb_min_window_tiled_pcs.txt".format(i, genome,window_size)) as file:
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
         
            for k in [14] : ## Set the best reference chromosome index if it follows -4 power law [1,2,3,12,14]
                pcs1=[]
                pcs1.append(lines_frame.iloc[k].values.tolist())
                flat_pcs1 = [item for sublist in pcs1 for item in sublist]
                final_pcs1 = [int(k) for k in ','.join(flat_pcs1).split(',')]        
                serial_count1 = []
                for j in range(0, 4000):
                    serial_count1.append(final_pcs1.count(j))
                p_i = np.array(serial_count1).astype(float) 
                p_i[np.where(p_i == 0)] = 0.000001
                
######## q_i ###############

for genome,name in zip(genome_list,names):
        chrom = []
        win = []
        entro = []
        for i in [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X','Y']:
            lines=[]
            path = '/bucket/MillerU/Abrar/PCS/hg38-{0}/relative_entropy/'.format(genome,window_size,i)
            #os.makedirs(path)
            with open("/bucket/MillerU/Abrar/hg38-{1}_outputs/hg38-{1}_tiled_pcs_{2}kb/chr{0}_{2}kb_min_window_tiled_pcs.txt".format(i, genome,window_size)) as file:
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
                q_i = np.array(serial_count2).astype(float) 
                q_i[np.where(q_i == 0)] = 0.000001
                chrom.append(i)
                win.append(k+1)
                entro.append(entropy(p_i, q_i))
                                
                #print(p_i,q_i,entropy(p_i, q_i))
                
                
rl_entro_df = pd.DataFrame(list(zip(chrom, win, entro)), columns =['chrom', 'window number', 'relative entropy'])
rl_entro_df.to_csv("/bucket/MillerU/Abrar/PCS/hg38-{0}/relative_entropy/relative_entropy_{1}kb_hggg.csv".format(genome,window_size), index=False)
















        
        
        
    
    
    
