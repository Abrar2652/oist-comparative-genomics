

import os
import glob
import pandas as pd


#####################  LOG-LOG #############################
# PCS

# Preparing the data for the plot
import numpy as np
import matplotlib.pyplot as plt
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
# Resizing the figure
fig, ax = plt.subplots(1, figsize=(20, 13))

############ LOOPING THROUGH ALL THE INDIVIDUAL PLOTS ##################

#genomes = ['gorGor6','mm39','galVar1','micMur3','ponAbe3','rheMac10','nomLeu3','saiBol1']
genomes = ['calJac4','danRer10','galGal6','ornAna2','equCab3','oryCun2','petMar3','tarSyr2','thaSir1','triMan1','xenTro10']
#colors = ['blue','red','gold','green','olive','black','brown','darkorange']
colors = ['blue','red','darkorange','black','firebrick','green','gold','turquoise','indigo','deeppink','olive']

for genome, color in zip(genomes, colors):
        lines=[]
        for i in [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X','Y']:
            with open("/bucket/MillerU/Abrar/hg38-{0}_outputs/hg38-{0}_tiled_pcs_15kb/chr{1}_15kb_min_window_tiled_pcs.txt".format(genome, i)) as file:
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
        # kurtosis_100pc_500col.csv" in case of galVar1
        nash_kurtosis = pd.read_csv("/bucket/MillerU/Abrar/Moment_val_bins/hg38{0}/All_bin_15kb/quantile_kurtosis_100pc.csv".format(genome)) #in the Moment_val_bins folder, the kurtosis files are located
        nash_kurtosis = nash_kurtosis[nash_kurtosis['0'] > 0] #nash_kurtosis['0'] != 0  #Manually convert excel to numeric (upto 9 decimal places) if any error is encountered
        
        # In[92]:
        species1 = 'hg38'
        species2 = genome
        percent = '100' #change accordingly
        
        # In[40]:
        
        cutoff1 = 5 #10
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
        
        

        
        # Plotting the graph with Log ticks at x and y axis using loglog
        m = [x for x in range(int(max(final_tail_pcs1)+5))][30:]
        n = serial_count1[30:]
        ax.loglog(m, n, '{}'.format(color), marker='.', alpha=0.5, linewidth=2, linestyle='None', label='hg38/{0} right heavy-tail cutoff = {1}'.format(genome, cutoff1))
   
   
        
ax.set_title('Right Heavy-tailed Kurtosis PCS Length Distribution Comparisons (log-log) for L>=30', fontsize=26)
ax.set_xlabel('PCS Length (log10)', fontsize=40)
ax.set_ylabel('Counts (log10)', fontsize=40)

########## Majoring and Minoring #############

xfmt = myformatter(labelOnlyBase=False, minor_thresholds=(4, 0.2))
yfmt = myformatter(labelOnlyBase=False, minor_thresholds=(10, 0.05)) 

ax.set_xlim([3*(10**1), 10**4])
ax.set_ylim([0.5, 10**6])

locmin = LogLocator(base=10.0, subs=np.arange(2, 10) * .1, numticks=20)
ax.xaxis.set_minor_locator(locmin)
ax.xaxis.set_minor_formatter(xfmt)
ax.yaxis.set_minor_locator(locmin)
ax.yaxis.set_minor_formatter(yfmt)
#for axis in ['bottom', 'left']:
#  ax.spines[axis].set_linewidth(2.5)
###################################################
ax.tick_params(axis='both', which='major', labelsize=28, width=2.5)
ax.tick_params(axis='y', which='minor', labelsize=12, width=2.5)
ax.tick_params(axis='x', which='minor', labelsize=13, width=2.5)
ax.legend(fontsize=22, markerscale=4)

plt.savefig('/bucket/MillerU/Abrar/PCS_comparison/Right_kurtosis_cutoff{0}_PCS_comparison_loglog.png'.format(cutoff1))
               
               

        
        
        




















