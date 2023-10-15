

## Red in case of Mouse, Blue in case of Gorilla
import os
import glob
import pandas as pd

lines=[]
for i in [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X','Y']:
    with open("/bucket/MillerU/Abrar/hg19-mm10_outputs/hg19-mm10_tiled_pcs_100percent_vanilla_L20/chr{0}_30kb_min_window_tiled_pcs.txt".format(i)) as file:
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
# kurtosis_100pc_500col.csv" in case of mm10
nash_kurtosis = pd.read_csv("/bucket/MillerU/Abrar/Moment_val_bins/hg19mm10/All_bin_30kb_vanilla_L20/quantile_kurtosis_100pc.csv") #in the Moment_val_bins folder, the kurtosis files are located
nash_kurtosis = nash_kurtosis[nash_kurtosis['0'] > 0] #nash_kurtosis['0'] != 0  #Manually convert excel to numeric (upto 9 decimal places) if any error is encountered

# In[92]:

species1 = 'hg19'
species2 = 'mm10'
percent = '100' #change accordingly

# In[40]:
# 13->Mouse, 12.5->Mouse

cutoff1 = 10
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



cutoff2 = 1.5
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


cutoff3 = 20
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
        


cutoff4 = 4
backtrack_index4 = nash_kurtosis[nash_kurtosis['0'] > cutoff4].index.tolist() #take decision while selecting this 'cutoff' #high kurtosis
#backtrack_index4 = nash_kurtosis[nash_kurtosis['0'] < cutoff4].index.tolist() #take decision while selecting this 'cutoff' #low kurtosis

tail_pcs4=[]
for i in backtrack_index4:
    tail_pcs4.append(lines_frame.iloc[i].values.tolist())
flat_tail_pcs4 = [item for sublist in tail_pcs4 for item in sublist]
final_tail_pcs4 = [int(i) for i in ','.join(flat_tail_pcs4).split(',')]        

serial_count4 = []
for i in range(0, int(max(final_tail_pcs4)+5)):
    serial_count4.append(final_tail_pcs4.count(i))
    
   
cutoff5 = 6
backtrack_index5 = nash_kurtosis[nash_kurtosis['0'] > cutoff5].index.tolist() #take decision while selecting this 'cutoff' #high kurtosis
#backtrack_index5 = nash_kurtosis[nash_kurtosis['0'] < cutoff5].index.tolist() #take decision while selecting this 'cutoff' #low kurtosis

tail_pcs5=[]
for i in backtrack_index5:
    tail_pcs5.append(lines_frame.iloc[i].values.tolist())
flat_tail_pcs5 = [item for sublist in tail_pcs5 for item in sublist]
final_tail_pcs5 = [int(i) for i in ','.join(flat_tail_pcs5).split(',')]        

serial_count5 = []
for i in range(0, int(max(final_tail_pcs5)+5)):
    serial_count5.append(final_tail_pcs5.count(i))
         
################# LINEAR #########################
# PCS
# Importing necessary libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Preparing the data for the plot

dfx = pd.read_csv("/bucket/MillerU/Abrar/PCS/hg19-mm10/PCS_hg19mm10.csv")
x = [i for i in dfx.index][20:]
y = dfx['Number of PCS of Length L'][20:]
m = [x for x in range(int(max(final_tail_pcs1)+5))][20:] #high kurtosis
n = serial_count1[20:]
q = [x for x in range(int(max(final_tail_pcs2)+5))][20:] #high kurtosis
r = serial_count2[20:]
a = [x for x in range(int(max(final_tail_pcs3)+5))][20:] #high kurtosis
b = serial_count3[20:]
s = [x for x in range(int(max(final_tail_pcs4)+5))][20:] #high kurtosis
t = serial_count4[20:]
c = [x for x in range(int(max(final_tail_pcs5)+5))][20:] #high kurtosis
d = serial_count5[20:]

# Resizing the figure
plt.figure(figsize=[10, 7])

# Plotting the graph with Log ticks at x and y axis using loglog
plt.plot(x, y, '.r', linewidth=2, label='Human/Mouse whole genome PCS')
plt.plot(m ,n, '+g',linewidth=2, label='Frequency of right heavy tailed kurtosis PCS (cutoff = {})'.format(cutoff1))
plt.plot(q ,r, '*b',linewidth=2, label='Frequency of left heavy tailed kurtosis PCS (cutoff = {})'.format(cutoff2))
plt.plot(a, b, '.m', linewidth=2, label='Frequency of right heavy tailed kurtosis PCS (cutoff = {})'.format(cutoff3))
plt.plot(s, t, 'gold', marker='*', linewidth=2, linestyle='None', label='Frequency of right heavy tailed kurtosis PCS (cutoff = {})'.format(cutoff4))
plt.plot(c, d, 'brown', marker='.', linewidth=2, linestyle='None', label='Frequency of right heavy tailed kurtosis PCS (cutoff = {})'.format(cutoff5))

plt.title('Heavy-tailed vanilla PCS length distribution comparisons of\n hg19/mm10 genome alignment for L>=20', fontsize=14)

plt.xlabel('PCS Length', fontsize=20)
plt.ylabel('Counts', fontsize=20)
plt.tick_params(axis='x', labelsize=15)
plt.tick_params(axis='y', labelsize=15)

plt.legend(fontsize=13, markerscale=3)
plt.savefig("/bucket/MillerU/Abrar/PCS_comparison/hg19mm10/PCS_comparison_hg19mm10_linear_L20.png")



#####################  LOG-LOG #############################
# PCS

# Preparing the data for the plot
from matplotlib.ticker import FormatStrFormatter,LogFormatter,ScalarFormatter,FuncFormatter,LogLocator,NullFormatter

##########################
class scalar(ScalarFormatter):
    def _formatSciNotation(self, s):
        # transform 1e+004 into 1e4, for example
        if self._useLocale:
            decimal_point = locale.localeconv()['decimal_point']
            positive_sign = locale.localeconv()['positive_sign']
        else:
            decimal_point = '.'
            positive_sign = '+'
        tup = s.split('e')
        try:
            significand = tup[0].rstrip('0').rstrip(decimal_point)
            sign = tup[1][0].replace(positive_sign, '')
            exponent = tup[1][1:].lstrip('0')
            if self._useMathText or self._usetex:
                if significand == '1' and exponent != '':
                    # reformat 1x10^y as 10^y
                    significand = ''
                if exponent:
                    exponent = '10^{%s%s}' % (sign, exponent)
                if significand and exponent:
                    return r'%s{\times}%s' % (significand, exponent)
                else:
                    return r'%s%s' % (significand, exponent)
            else:
                s = ('%se%s%s' % (significand, sign, exponent)).rstrip('e')
                return s
        except IndexError:
            return s

#########################
f = scalar(useOffset=False, useMathText=True)
g = lambda x,pos : "${}$".format(f._formatSciNotation('%1.10e' % x))
fmt = FuncFormatter(g)

class myformatter(LogFormatter):
    def _num_to_string(self, x, vmin, vmax):
        if x >= 1:
            s = fmt(x)
        elif x < 1 and x >= 0.001:
            s = fmt(x)
        return s
############################

x = [i for i in dfx.index][20:]
y = dfx['Number of PCS of Length L'][20:]
m = [x for x in range(int(max(final_tail_pcs1)+5))][20:] #high kurtosis
n = serial_count1[20:]
q = [x for x in range(int(max(final_tail_pcs2)+5))][20:] #high kurtosis
r = serial_count2[20:]
a = [x for x in range(int(max(final_tail_pcs3)+5))][20:] #high kurtosis
b = serial_count3[20:]
s = [x for x in range(int(max(final_tail_pcs4)+5))][20:] #high kurtosis
t = serial_count4[20:]
c = [x for x in range(int(max(final_tail_pcs5)+5))][20:] #high kurtosis
d = serial_count5[20:]

# Resizing the figure
fig, ax = plt.subplots(1, figsize=(20, 13))

# Plotting the graph with Log ticks at x and y axis using loglog
plt.loglog(x, y, '.r', linewidth=2, label='Human/Mouse whole genome PCS')
plt.loglog(m, n, '+g',linewidth=2, label='Frequency (log10) of right heavy tailed kurtosis PCS (log10) (cutoff = {})'.format(cutoff1))
plt.loglog(q, r, '*b',linewidth=2, label='Frequency (log10) of left heavy tailed kurtosis PCS (log10) (cutoff = {})'.format(cutoff2))
plt.loglog(a, b, '.m', linewidth=2, label='Frequency (log10) of right heavy tailed kurtosis PCS (log10) (cutoff = {})'.format(cutoff3))
plt.loglog(s, t, 'gold', marker='*', linewidth=2, linestyle='None', label='Frequency (log10) of right heavy tailed kurtosis PCS (log10) (cutoff = {})'.format(cutoff4))
plt.loglog(c, d, 'brown', marker='.', linewidth=2, linestyle='None', label='Frequency (log10) of right heavy tailed kurtosis PCS (log10) (cutoff = {})'.format(cutoff5))

plt.title('Heavy-tailed vanilla PCS length distribution comparisons of\n hg19/mm10 genome alignment (log-log) for L>=20', fontsize=26)
ax.set_xlabel('PCS Length (log10)', fontsize=40)
ax.set_ylabel('Counts (log10)', fontsize=40)

########## Majoring and Minoring #############

xfmt = myformatter(labelOnlyBase=False, minor_thresholds=(4, 0.2))
yfmt = myformatter(labelOnlyBase=False, minor_thresholds=(10, 0.05)) 

ax.set_xlim([2*(10**1), 10**4])
ax.set_ylim([0.5, 10**6])

locmin = LogLocator(base=10.0, subs=np.arange(2, 10) * .1, numticks=20)
ax.xaxis.set_minor_locator(locmin)
ax.xaxis.set_minor_formatter(xfmt)
ax.yaxis.set_minor_locator(locmin)
ax.yaxis.set_minor_formatter(yfmt)


###################################################

ax.tick_params(axis='both', which='major', labelsize=28)
ax.tick_params(axis='y', which='minor', labelsize=12)
ax.tick_params(axis='x', which='minor', labelsize=13)
ax.legend(fontsize=22, markerscale=4)
plt.savefig('/bucket/MillerU/Abrar/PCS_comparison/hg19mm10/PCS_comparison_hg19mm10_loglog_L20.png')



################### SEMILOG Y ################################
# PCS
# Preparing the data for the plot
x = [i for i in dfx.index][20:]
y = dfx['Number of PCS of Length L'][20:]
m = [x for x in range(int(max(final_tail_pcs1)+5))][20:] #high kurtosis
n = serial_count1[20:]
q = [x for x in range(int(max(final_tail_pcs2)+5))][20:] #high kurtosis
r = serial_count2[20:]
a = [x for x in range(int(max(final_tail_pcs3)+5))][20:] #high kurtosis
b = serial_count3[20:]
s = [x for x in range(int(max(final_tail_pcs4)+5))][20:] #high kurtosis
t = serial_count4[20:]
c = [x for x in range(int(max(final_tail_pcs5)+5))][20:] #high kurtosis
d = serial_count5[20:]

plt.figure(figsize=[10, 7])

plt.semilogy(x, y, '.r', linewidth=2, base=10, label='Human/Mouse whole genome PCS')
plt.semilogy(m, n, '+g', linewidth=2, base=10, label='Frequency (log10) of right heavy tailed kurtosis PCS (cutoff = {})'.format(cutoff1))
plt.semilogy(q, r, '*b', linewidth=2, base=10, label='Frequency (log10) of left heavy tailed kurtosis PCS (cutoff = {})'.format(cutoff2))
plt.semilogy(a, b, '.m', linewidth=2, base=10, label='Frequency (log10) of right heavy tailed kurtosis PCS (cutoff = {})'.format(cutoff3))
plt.semilogy(s, t, 'gold', marker='*', linewidth=2, linestyle='None', base=10, label='Frequency (log10) of right heavy tailed kurtosis PCS (cutoff = {})'.format(cutoff4))
plt.semilogy(c, d, 'brown', marker='.', linewidth=2, linestyle='None', base=10, label='Frequency (log10) of right heavy tailed kurtosis PCS (cutoff = {})'.format(cutoff5))

plt.xlabel('PCS Length', fontsize=20)
plt.ylabel('Counts (log10)', fontsize=20)
plt.tick_params(axis='x', labelsize=15)
plt.tick_params(axis='y', labelsize=15)

plt.title("Heavy-tailed vanilla PCS length distribution comparisons of\n hg19/mm10 genome alignment (semi-log) for L>=20", fontsize=14)

plt.legend(fontsize=13, markerscale=3)
plt.savefig('/bucket/MillerU/Abrar/PCS_comparison/hg19mm10/PCS_comparison_hg19mm10_semilog_L20.png')






















