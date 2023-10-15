

## Red in case of Mouse, Blue in case of Gorilla
import os
import glob
import pandas as pd

lines=[]
for i in [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X','Y']:
    with open("/bucket/MillerU/Abrar/hg19-gorGor3_outputs/hg19-gorGor3_tiled_pcs_100percent/chr{0}_30kb_min_window_tiled_pcs.txt".format(i)) as file:
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

import pandas as pd
# kurtosis_100pc_500col.csv" in case of gorGor3
nash_kurtosis = pd.read_csv("/bucket/MillerU/Abrar/Moment_val_bins/hg19gorGor3/All_bin_30kb/kurtosis_100pc_500col.csv") #in the Moment_val_bins folder, the kurtosis files are located
nash_kurtosis = nash_kurtosis[nash_kurtosis['0'] > 0] #nash_kurtosis['0'] != 0  #Manually convert excel to numeric (upto 9 decimal places) if any error is encountered

# In[92]:
species1 = 'hg19'
species2 = 'gorGor3'
percent = '100' #change accordingly

# In[40]:
# 13->Gorilla, 12.5->Mouse

cutoff1 = 8.75
backtrack_index1 = nash_kurtosis[nash_kurtosis['0'] > cutoff1].index.tolist() #take decision while selecting this 'cutoff' #high kurtosis
#backtrack_index1 = nash_kurtosis[nash_kurtosis['0'] < cutoff1].index.tolist() #take decision while selecting this 'cutoff' #low kurtosis

tail_pcs1=[]
for i in backtrack_index1:
    tail_pcs1.append(lines_frame.iloc[i].values.tolist())
flat_tail_pcs1 = [item for sublist in tail_pcs1 for item in sublist]
final_tail_pcs1 = [int(i) for i in ','.join(flat_tail_pcs1).split(',')]        

serial_count1 = []
for i in range(0, int(max(final_tail_pcs1)+5)):
    serial_count1.append(final_tail_pcs1.count(i))


cutoff2 = 1.25
#backtrack_index2 = nash_kurtosis[nash_kurtosis['0'] > cutoff2].index.tolist() #take decision while selecting this 'cutoff' #high kurtosis
backtrack_index2 = nash_kurtosis[nash_kurtosis['0'] < cutoff2].index.tolist() #take decision while selecting this 'cutoff' #low kurtosis

tail_pcs2=[]
for i in backtrack_index2:
    tail_pcs2.append(lines_frame.iloc[i].values.tolist())
flat_tail_pcs2 = [item for sublist in tail_pcs2 for item in sublist]
final_tail_pcs2 = [int(i) for i in ','.join(flat_tail_pcs2).split(',')]        

serial_count2 = []
for i in range(0, int(max(final_tail_pcs2)+5)):
    serial_count2.append(final_tail_pcs2.count(i))


cutoff3 = 10.625
backtrack_index3 = nash_kurtosis[nash_kurtosis['0'] > cutoff3].index.tolist() #take decision while selecting this 'cutoff' #high kurtosis
#backtrack_index3 = nash_kurtosis[nash_kurtosis['0'] < cutoff3].index.tolist() #take decision while selecting this 'cutoff' #low kurtosis

tail_pcs3=[]
for i in backtrack_index3:
    tail_pcs3.append(lines_frame.iloc[i].values.tolist())
flat_tail_pcs3 = [item for sublist in tail_pcs3 for item in sublist]
final_tail_pcs3 = [int(i) for i in ','.join(flat_tail_pcs3).split(',')]        

serial_count3 = []
for i in range(0, int(max(final_tail_pcs3)+5)):
    serial_count3.append(final_tail_pcs3.count(i))
        
    
################# LINEAR #########################
# PCS
# Importing necessary libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Preparing the data for the plot

dfx = pd.read_csv("/bucket/MillerU/Abrar/PCS/hg19-gorGor3/PCS_hg19gorGor3.csv")
x = [i for i in dfx.index]
y = dfx['Number of PCS of Length L']
m = [x for x in range(int(max(final_tail_pcs1)+5))][1:] #high kurtosis
n = serial_count1[1:]
q = [x for x in range(int(max(final_tail_pcs2)+5))][1:] #low kurtosis
r = serial_count2[1:]
a = [x for x in range(int(max(final_tail_pcs3)+5))][1:] #high kurtosis
b = serial_count3[1:]

# Resizing the figure
plt.figure(figsize=[15, 10])

# Plotting the graph with Log ticks at x and y axis using loglog
plt.plot(x, y, '.r', linewidth=2, label='Human/Gorilla whole genome PCS')
plt.plot(m ,n, '+g',linewidth=2, label='Frequency of heavy tailed kurtosis PCS (cutoff = {})'.format(cutoff1))
plt.plot(q ,r, '*b',linewidth=2, label='Frequency of left clustered kurtosis PCS (cutoff = {})'.format(cutoff2))
plt.plot(a, b, '.m', linewidth=2, label='Frequency of heavy tailed kurtosis PCS (cutoff = {})'.format(cutoff3))

plt.title('Length distribution comparisons of Perfectly Conserved L-mers from hg19/gorGor3 genome alignment for L>=1', fontsize=15)

plt.xlabel('PCS', fontsize=13)
plt.ylabel('Counts', fontsize=13)

plt.legend()
plt.savefig("/bucket/MillerU/Abrar/PCS_comparison/hg19gorGor3/PCS_comparison_hg19gorGor3_linear.png")



#####################  LOG-LOG #############################
# PCS

# Preparing the data for the plot

x = [i for i in dfx.index]
y = dfx['Number of PCS of Length L']
m = [x for x in range(int(max(final_tail_pcs1)+5))][1:] #high kurtosis
n = serial_count1[1:]
q = [x for x in range(int(max(final_tail_pcs2)+5))][1:] #low kurtosis
r = serial_count2[1:]
a = [x for x in range(int(max(final_tail_pcs3)+5))][1:] #high kurtosis
b = serial_count3[1:]

# Resizing the figure
plt.figure(figsize=[15, 10])

# Plotting the graph with Log ticks at x and y axis using loglog
plt.loglog(x, y, '.r', linewidth=2, label='Human/Gorilla whole genome PCS')
plt.loglog(m, n, '+g',linewidth=2, label='Frequency (log10) of heavy tailed kurtosis PCS (log10) (cutoff = {})'.format(cutoff1))
plt.loglog(q, r, '*b',linewidth=2, label='Frequency (log10) of left clustered kurtosis PCS (log10) (cutoff = {})'.format(cutoff2))
plt.loglog(a, b, '.m', linewidth=2, label='Frequency (log10) of heavy tailed kurtosis PCS (log10) (cutoff = {})'.format(cutoff3))

plt.title('Length distribution comparisons of Perfectly Conserved L-mers from hg19/gorGor3 genome alignment (log-log) for L>=1', fontsize=15)
plt.xlabel('PCS (log10)', fontsize=13)
plt.ylabel('Counts (log10)', fontsize=13)

plt.legend()
plt.savefig('/bucket/MillerU/Abrar/PCS_comparison/hg19gorGor3/PCS_comparison_hg19gorGor3_loglog_L1.png')



################### SEMILOG Y ################################
# PCS
# Preparing the data for the plot
x = [i for i in dfx.index]
y = dfx['Number of PCS of Length L']
m = [x for x in range(int(max(final_tail_pcs1)+5))][1:] #high kurtosis
n = serial_count1[1:]
q = [x for x in range(int(max(final_tail_pcs2)+5))][1:] #low kurtosis
r = serial_count2[1:]
a = [x for x in range(int(max(final_tail_pcs3)+5))][1:] #high kurtosis
b = serial_count3[1:]

plt.figure(figsize=[15, 10])

plt.semilogy(x, y, '.r', linewidth=2, base=10, label='Human/Gorilla whole genome PCS')
plt.semilogy(m, n, '+g', linewidth=2, base=10, label='Frequency (log10) of heavy tailed kurtosis PCS (cutoff = {})'.format(cutoff1))
plt.semilogy(q, r, '*b', linewidth=2, base=10, label='Frequency (log10) of left clustered kurtosis PCS (cutoff = {})'.format(cutoff2))
plt.semilogy(a, b, '.m', linewidth=2, base=10, label='Frequency (log10) of heavy tailed kurtosis PCS (cutoff = {})'.format(cutoff3))

plt.xlabel('PCS', fontsize=13)
plt.ylabel('Counts (log10)', fontsize=13)
plt.title("Length distribution comparisons of Perfectly Conserved L-mers from hg19/gorGor3 genome alignment (semi-log) for L>=1", fontsize=15)

plt.legend()
plt.savefig('/bucket/MillerU/Abrar/PCS_comparison/hg19gorGor3/PCS_comparison_hg19gorGor3_semilog.png')























