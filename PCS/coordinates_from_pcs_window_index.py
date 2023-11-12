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
        
window_size  = 50       # [10000, 5000, 2500, 1200, 500, 100, 50]              
############## TILED COORDINATES ##################
#[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X','Y']
lines1=[]
for i in [10]:
    with open("/bucket/MillerU/Abrar/hg38-gorGor6_outputs/hg38-gorGor6_tiled_coordinates_{1}kb/chr{0}_{1}kb_min_window_tiled_coordinates.txt".format(i,window_size)) as file:
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
print(len(lines1)) # 5 for chrY 10MB

#backtrack_index = [0,1,2,3,4] # window index
#backtrack_index = [14,18,30,54,55,56,61,188,256,257,258,259,260,262,264,267,330,332,333,334,429,430,592,622,1268,2074,2075,2417,2502,2503,2504,2866,2870,2871,2872,2887,2893,2906,2907,2908,2909,2910,2922,2928,2961,2963,2971,2972,2981,2982,2887,2990,2991,2992,2993,3045,3047,3103,3146,3230,4058,4152,4553,4572,4574,4800,4931,4964,4966,4970]
#backtrack_index = [i-1 for i in backtrack_index]

# empty windows will be skipped
for i in range(len(lines1)):
    #print(i)    
    if(lines1[i] != [[0]]):
        chrom=[]
        for j in range(len(lines1[i])):
            chrom.append('chr10')
                      
        my_df = pd.DataFrame(lines1[i], columns=['start','end'])
        my_df['chrom'] = chrom
        my_df = my_df[['chrom', 'start', 'end']]
        #print(my_df.start.iloc[0],my_df.end.iloc[-1]) ## only for the special windows 
        my_df.to_csv("/bucket/.mabuya/MillerU/Abrar/PCS/hg38-gorGor6/all_window_comparisons/{0}kb/chr10_idx{1}_{0}kb_coordinates_hg38gorGor6.csv".format(window_size, i), index=False)
    else:
        continue













