

import numpy as np
import matplotlib.pyplot as plt
import csv
import pandas as pd
import scipy
import lmfit
import os
import glob

# General Functions
def linear(x, m, c):
    return c / (x ** (m))
    
#####################  LOG-LOG SHADED #############################
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
#fig, ax = plt.subplots(2, figsize=(28, 30)) #28, 29

############ LOOPING THROUGH ALL THE INDIVIDUAL PLOTS ##################


genomes = ['gorGor6','mm39','galVar1','micMur3','ponAbe3','rheMac10','nomLeu3','saiBol1']
colors = ['blue','red','green','gold','olive','black','brown','darkorange']

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
        
        cutoff1 = 5
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
        m = [x for x in range(int(max(final_tail_pcs1)+5))][20:]
        n = serial_count1[20:]

        ax.loglog(m, n, '{}'.format(color), marker='.', linewidth=2, alpha=0.5, linestyle='None', label='hg38/{0} right heavy-tail cutoff = {1}'.format(genome, cutoff1))
           
        if (genome == 'mm39'):
              ######## SLOPE plot #########
              dict1 = {'PCS length (L)': m, 'Number of PCS of Length L': n}  
              df = pd.DataFrame(dict1)
              print(df.head())
              df = df[df['Number of PCS of Length L']!=0]
              x = df['PCS length (L)']
              x_fake = df['PCS length (L)']
              y = df['Number of PCS of Length L']
              y_fake = df['Number of PCS of Length L']
              x=np.array(x)
              x_fake=np.array(x_fake)
              y=np.array(y)
              y_fit = linear(x_fake, 4, 10**10.7)
              absError = np.log10(y_fit) - np.log10(y_fake)
              SE = np.square(absError) # squared errors
              MSE = np.mean(SE) # mean squared errors
              RMSE = np.sqrt(MSE) # Root Mean Squared Error, RMSE
              Rsquared = 1.0 - (np.var(absError) / np.var(np.log10(y_fake)))
              ax.annotate("Linear Curve:\nMSE: {:.3f}\nRMSE: {:.3f}\nR-squared: {:.3f}".format(MSE,RMSE,Rsquared),(9*10**2,6*10**1), fontsize=22, xytext=(900,6*10**1))
              ax.loglog(x, y_fit, 'b--', label="Fitted Curve: Slope = -4")
              
        if (genome == 'gorGor6'):
              ######## SLOPE plot #########  [-0.18336822  1.85262218 -6.90339853 20.08558534]
              dict1 = {'PCS length (L)': m, 'Number of PCS of Length L': n}  
              df = pd.DataFrame(dict1)
              print(df.head())
              df = df[df['Number of PCS of Length L']!=0]
              x = df['PCS length (L)']
              x_fake = df['PCS length (L)'][:1030]
              y_fake = df['Number of PCS of Length L'][:1030]
              logx, logy = np.log(x_fake), np.log(y_fake)
              p = np.polyfit(logx, logy, 3)
              y_fit = np.exp(np.polyval(p, logx)) #plyval->Compute polynomial values.
              absError = np.log(y_fit) - np.log(y_fake)
              SE = np.square(absError) # squared errors
              MSE = np.mean(SE) # mean squared errors
              RMSE = np.sqrt(MSE) # Root Mean Squared Error, RMSE
              Rsquared = 1.0 - (np.var(absError) / np.var(np.log10(y_fake)))
             
              ax.annotate("Polynomial Curve:\nMSE: {:.3f}\nRMSE: {:.3f}\nR-squared: {:.3f}".format(MSE,RMSE,Rsquared),(3*10**3,6*10**1), fontsize=22, xytext=(3000,6*10**1))
              ax.loglog(x_fake, y_fit, 'r--', label="Fitted Curve: Degree = 3")        
   
######################

ax.axvspan(20, 60, color="b", alpha=0.2)  ######## change

########################
        
ax.set_title('Right Heavy-tailed Kurtosis PCS Length Distribution Comparisons (log-log) for L>=20', fontsize=26)

#### for making the left and bottom spines bold ########
#for axis in ['bottom', 'left']:
#  ax.spines[axis].set_linewidth(2.5)
  
ax.set_xlabel('PCS Length (log10)', fontsize=40)
ax.set_ylabel('Counts (log10)', fontsize=40)

########## Majoring and Minoring #############

xfmt = myformatter(labelOnlyBase=False, minor_thresholds=(4, 0.2))
yfmt = myformatter(labelOnlyBase=False, minor_thresholds=(10, 0.05)) 

ax.set_xlim([2*(10**1), 10**4])
ax.set_ylim([0.5, 10**7])

locmin = LogLocator(base=10.0, subs=np.arange(2, 10) * .1, numticks=20)
ax.xaxis.set_minor_locator(locmin)
ax.xaxis.set_minor_formatter(xfmt)
ax.yaxis.set_minor_locator(locmin)
ax.yaxis.set_minor_formatter(yfmt)

###################################################
ax.tick_params(axis='both', which='major', labelsize=28, width=2.5)
ax.tick_params(axis='y', which='minor', labelsize=12, width=2.5)
ax.tick_params(axis='x', which='minor', labelsize=13, width=2.5)
ax.legend(fontsize=22, markerscale=4)

plt.savefig('/bucket/MillerU/Abrar/PCS_comparison/Right_kurtosis_cutoff{0}_PCS_blue_shaded.png'.format(cutoff1))

#################### Zoomed In #######################################

fig, ax = plt.subplots(1, figsize=(20, 13))
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
        
        cutoff1 = 5
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
        m = [x for x in range(int(max(final_tail_pcs1)+5))][20:]
        n = serial_count1[20:]

        xlims = [20, 60] ############# change
        ylims = [0.5, 10**7]
       
        ax.loglog(m, n, '{}'.format(color), marker='.', alpha=0.5, linewidth=2, markersize=20, linestyle='None', label='hg38/{0} right heavy-tail cutoff = {1}'.format(genome, cutoff1))
        
        if (genome == 'mm39'):
              ######## SLOPE plot #########
              dict1 = {'PCS length (L)': m, 'Number of PCS of Length L': n}  
              df = pd.DataFrame(dict1)
              df = df[df['Number of PCS of Length L']!=0]
              x = df['PCS length (L)']
              x_fake = df['PCS length (L)']
              y = df['Number of PCS of Length L']
              y_fake = df['Number of PCS of Length L']
              x=np.array(x)
              x_fake=np.array(x_fake)
              y=np.array(y)
              y_fit = linear(x_fake, 4, 10**10.7)
              absError = np.log10(y_fit) - np.log10(y_fake)
              SE = np.square(absError) # squared errors
              MSE = np.mean(SE) # mean squared errors
              RMSE = np.sqrt(MSE) # Root Mean Squared Error, RMSE
              Rsquared = 1.0 - (np.var(absError) / np.var(np.log10(y_fake)))
              ax.annotate("Linear Curve:\nMSE: {:.3f}\nRMSE: {:.3f}\nR-squared: {:.3f}".format(MSE,RMSE,Rsquared),(2.5*10**1,3), fontsize=22, xytext=(25,3))#10**6
              ax.loglog(x, y_fit, 'b--', label="Fitted Curve: Slope = -4")      
               
        if (genome == 'gorGor6'):
              ######## SLOPE plot #########
              dict1 = {'PCS length (L)': m, 'Number of PCS of Length L': n}  
              df = pd.DataFrame(dict1)
              df = df[df['Number of PCS of Length L']!=0]
              x = df['PCS length (L)']
              x_fake = df['PCS length (L)'][:1030]
              y = df['Number of PCS of Length L']
              y_fake = df['Number of PCS of Length L'][:1030]
          
              logx, logy = np.log(x_fake), np.log(y_fake)
              p = np.polyfit(logx, logy, 3)
              y_fit = np.exp(np.polyval(p, logx)) #plyval->Compute polynomial values.
              absError = np.log(y_fit) - np.log(y_fake)
              SE = np.square(absError) # squared errors
              MSE = np.mean(SE) # mean squared errors
              RMSE = np.sqrt(MSE) # Root Mean Squared Error, RMSE
              Rsquared = 1.0 - (np.var(absError) / np.var(np.log10(y_fake)))
              ax.annotate("Polynomial Curve:\nMSE: {:.3f}\nRMSE: {:.3f}\nR-squared: {:.3f}".format(MSE,RMSE,Rsquared),(3.5*10**1,3), fontsize=22, xytext=(35,3))#10**6
              ax.loglog(x_fake, y_fit, 'r--', label="Fitted Curve: Degree = 3")   
                              
        ax.set_title('PCS for 20<=L<=60 (zoomed in)', fontsize=32) # change

        xfmt = myformatter(labelOnlyBase=False, minor_thresholds=(4, 0.2))#3,4
        yfmt = myformatter(labelOnlyBase=False, minor_thresholds=(10, 0.6))   
    
        ax.set_xlabel('PCS Length (log10)', fontsize=40)
        ax.set_ylabel('Counts (log10)', fontsize=40)
        ax.set_xlim(xlims)
        ax.set_ylim(ylims)
    
        locmaj = LogLocator(base=10.0, subs=(1.0, ), numticks=100)
        ax.xaxis.set_major_locator(locmaj)
    
        locmin = LogLocator(base=10.0, subs=np.arange(2, 10) * .1, numticks=10000)
        locmin2 = LogLocator(base=10.0, subs=np.arange(2, 10) * .1, numticks=20)
        ax.yaxis.set_minor_locator(locmin2)
        ax.yaxis.set_minor_formatter(yfmt)
        ax.xaxis.set_minor_locator(locmin)
        ax.xaxis.set_minor_formatter(xfmt)
        #for axis in ['bottom', 'left']:
        #  ax.spines[axis].set_linewidth(2.5)
        ax.tick_params(axis='both', which='major', labelsize=29, width=2.5)
        ax.tick_params(axis='y', which='minor', labelsize=17, width=2.5)
        ax.tick_params(axis='x', which='minor', labelsize=19, width=2.5)
        plt.xticks(weight='bold')
        ax.grid()


plt.savefig('/bucket/MillerU/Abrar/PCS_comparison/Right_kurtosis_cutoff{0}_PCS_zoomed_in.png'.format(cutoff1))






















