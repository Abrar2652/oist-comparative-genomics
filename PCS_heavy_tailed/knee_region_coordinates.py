#!/usr/bin/env python
# coding: utf-8

## Red in case of Mouse, Blue in case of Gorilla
import os
import glob
import pandas as pd

##### FUNCTION FOR CREATING CHUNKS OF LISTS ####
def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]
        
############## TILED CHROMOSOMES ##################
lines2=[]
for i in [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X','Y']:
    with open("/bucket/MillerU/Abrar/hg38-gorGor3_outputs/hg38-gorGor3_tiled_chrom_100percent/chr{0}_15kb_min_window_tiled_chrom.txt".format(i)) as file:
        while True:
            # Get next line from file
            line_i = file.readline().rstrip('\n')
            # if line is empty
            # end of file is reached
            if not line_i:
                break
            elif line_i==("\n"):
                continue
            line = [i for i in line_i.split(",")]
            lines2.append(line)
                    
############## TILED COORDINATES ##################
lines1=[]
for i in [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X','Y']:
    with open("/bucket/MillerU/Abrar/hg38-gorGor3_outputs/hg38-gorGor3_tiled_coordinates_100percent/chr{0}_15kb_min_window_tiled_coordinates.txt".format(i)) as file:
        while True:
            # Get next line from file
            line_i = file.readline().rstrip('\n')
            # if line is empty
            # end of file is reached
            if not line_i:
                break
            elif line_i==("\n"):
                continue
            line_j = [int(i) for i in line_i.split(",")]
            line = list(chunks(line_j, 2))
            lines1.append(line)

#### TILED PCS #########
lines=[]
for i in [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X','Y']:
    with open("/bucket/MillerU/Abrar/hg38-gorGor3_outputs/hg38-gorGor3_tiled_pcs_100percent/chr{0}_15kb_min_window_tiled_pcs.txt".format(i)) as file:
        while True:
            # Get next line from file
            line_i = file.readline().rstrip('\n')
            # if line is empty
            # end of file is reached
            if not line_i:
                break
            elif line_i == ("\n"):
                continue
            line = [int(i) for i in line_i.split(",")]
            lines.append(line)


import pandas as pd
# kurtosis_100pc_500col.csv" in case of gorGor3
nash_kurtosis = pd.read_csv("/bucket/MillerU/Abrar/Moment_val_bins/hg38gorGor3/All_bin_30kb/quantile_kurtosis_100pc.csv") 
nash_kurtosis = nash_kurtosis[nash_kurtosis['0'] > 0] #nash_kurtosis['0'] != 0  #Manually convert excel to numeric (upto 9 decimal places) if any error is encountered

# In[92]:
species1 = 'hg38'
species2 = 'gorGor3'
percent = '100' #change accordingly

# In[40]:
cutoff = 6.25
backtrack_index = nash_kurtosis[nash_kurtosis['0'] > cutoff].index.tolist() #take decision while selecting this 'cutoff' #high kurtosis

# Making lines and lines1 of the same length by extracting their high kurtosis windows only
lines3 = [] # pcs
for i in backtrack_index:
    lines3.append(lines[i])
    
lines4 = [] # coordinates
for i in backtrack_index:
    lines4.append(lines1[i])
    
lines5 = [] # chroms
for i in backtrack_index:
    lines5.append(lines2[i])
    
## Collecting PCS those are in the range(440, 560) and storing their indices in respective windows to knee_indices
need = [x for x in range(370, 480+1)]   ####### CHANGE THIS RANGE FOR EACH CASES


knee_indices = []
for win in lines3:
    new = []
    for k, pcs in enumerate(win):
        if pcs in need:
            new.append(k)
    knee_indices.append(new)

# Finding the knee PCS coordinates from their relative indices stored in knee_indices
knee_coord = []
for i, x in enumerate(knee_indices):
    for j in x:
        knee_coord.append(lines4[i][j])

knee_chrom = []
for i, x in enumerate(knee_indices):
    for j in x:
        knee_chrom.append(lines5[i][j])
                
my_df = pd.DataFrame(knee_coord, columns=['Start','End'])
my_df['Chrom'] = knee_chrom
my_df = my_df[['Chrom', 'Start', 'End']]
my_df.to_csv('/bucket/.mabuya/MillerU/Abrar/PCS_heavy_tailed/knee_pcs_coordinates_hg38gorGor3.csv', index=False)














