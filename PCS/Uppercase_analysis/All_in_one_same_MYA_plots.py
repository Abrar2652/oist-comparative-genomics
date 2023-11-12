

from pandas._libs.algos import diff_2d
import os
import glob
import matplotlib.pyplot as plt
import csv
import pandas as pd
import scipy
import lmfit


df1 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/PCS/Uppercase_analysis/indel_terminated/vanilla_indel_terminated_PCS_hg38rheMac10.csv")
x1 = [i for i in df1.index][1:]
y1 = df1['Number of PCS of Length L'][1:]

df2 = pd.read_csv('/bucket/.mabuya/MillerU/Abrar/PCS/Uppercase_analysis/indel_terminated/vanilla_indel_terminated_PCS_hg38papAnu4.csv')
x2 = [i for i in df2.index][1:]
y2 = df2['Number of PCS of Length L'][1:]

df3 = pd.read_csv('/bucket/.mabuya/MillerU/Abrar/PCS/Uppercase_analysis/indel_terminated/vanilla_indel_terminated_PCS_hg38macFas5.csv')
x3 = [i for i in df3.index][1:]
y3 = df3['Number of PCS of Length L'][1:]

df4 = pd.read_csv('/bucket/.mabuya/MillerU/Abrar/PCS/Uppercase_analysis/indel_terminated/vanilla_indel_terminated_PCS_hg38chlSab2.csv')
x4 = [i for i in df4.index][1:]
y4 = df4['Number of PCS of Length L'][1:]

df5 = pd.read_csv('/bucket/.mabuya/MillerU/Abrar/PCS/Uppercase_analysis/indel_terminated/vanilla_indel_terminated_PCS_hg38nasLar1.csv')
x5 = [i for i in df5.index][1:]
y5 = df5['Number of PCS of Length L'][1:]

df6 = pd.read_csv('/bucket/.mabuya/MillerU/Abrar/PCS/Uppercase_analysis/indel_terminated/vanilla_indel_terminated_PCS_hg38cerAty1.csv')
x6 = [i for i in df6.index][1:]
y6 = df6['Number of PCS of Length L'][1:]

df7 = pd.read_csv('/bucket/.mabuya/MillerU/Abrar/PCS/Uppercase_analysis/indel_terminated/vanilla_indel_terminated_PCS_hg38manLeu1.csv')
x7 = [i for i in df7.index][1:]
y7 = df7['Number of PCS of Length L'][1:]

df8 = pd.read_csv('/bucket/.mabuya/MillerU/Abrar/PCS/Uppercase_analysis/indel_terminated/vanilla_indel_terminated_PCS_hg38macNem1.csv')
x8 = [i for i in df8.index][1:]
y8 = df8['Number of PCS of Length L'][1:]

df9 = pd.read_csv('/bucket/.mabuya/MillerU/Abrar/PCS/Uppercase_analysis/indel_terminated/vanilla_indel_terminated_PCS_hg38colAng1.csv')
x9 = [i for i in df9.index][1:]
y9 = df9['Number of PCS of Length L'][1:]

df10 = pd.read_csv('/bucket/.mabuya/MillerU/Abrar/PCS/Uppercase_analysis/indel_terminated/vanilla_indel_terminated_PCS_hg38rhiBie1.csv')
x10 = [i for i in df10.index][1:]
y10 = df10['Number of PCS of Length L'][1:]

df11 = pd.read_csv('/bucket/.mabuya/MillerU/Abrar/PCS/Uppercase_analysis/indel_terminated/vanilla_indel_terminated_PCS_hg38rhiRox1.csv')
x11 = [i for i in df11.index][1:]
y11 = df11['Number of PCS of Length L'][1:]

df12 = pd.read_csv('/bucket/.mabuya/MillerU/Abrar/PCS/Uppercase_analysis/indel_terminated/vanilla_indel_terminated_PCS_hg38gorGor6.csv')
x12 = [i for i in df12.index][1:]
y12 = df12['Number of PCS of Length L'][1:]

df13 = pd.read_csv('/bucket/.mabuya/MillerU/Abrar/PCS/Uppercase_analysis/indel_terminated/vanilla_indel_terminated_PCS_hg38mm39.csv')
x13 = [i for i in df13.index][1:]
y13 = df13['Number of PCS of Length L'][1:]

# Resizing the figure
plt.figure(figsize=[10, 7])
plt.xlim(0.7,3*10**3)
plt.ylim(0.5,10**8)

# Plotting the graph with Log ticks at x and y axis using loglog

plt.loglog(x1, y1, 'crimson', marker='.', linewidth=2, linestyle='None', label='Human/Rhesus')
plt.loglog(x2, y2, 'yellow', marker='.', linewidth=2, linestyle='None', label='Human/Baboon')
plt.loglog(x3, y3, 'olive', marker='.', linewidth=2, linestyle='None', label='Human/Crab-eating Macaque')
plt.loglog(x4, y4, 'teal', marker='.', linewidth=2, linestyle='None', label='Human/Green Monkey')
plt.loglog(x5, y5, 'red', marker='.', linewidth=2, linestyle='None', label='Human/Proboscis Monkey')
plt.loglog(x6, y6, 'green', marker='.', linewidth=2, linestyle='None', label='Human/Sooty Mangabey')
plt.loglog(x7, y7, 'blue', marker='.', linewidth=2, linestyle='None', label='Human/Drill')
plt.loglog(x8, y8, 'darkorange', marker='.', linewidth=2, linestyle='None', label='Human/Pig-tailed Macaque')
plt.loglog(x9, y9, 'gold', marker='.', linewidth=2, linestyle='None', label='Human/Angolan Colobus')
plt.loglog(x10, y10, 'firebrick', marker='.', linewidth=2, linestyle='None', label='Human/Black Snub-nosed Monkey')
plt.loglog(x11, y11, 'black', marker='.', linewidth=2, linestyle='None', label='Human/Golden Snub-nosed Monkey')
plt.loglog(x12, y12, 'purple', marker='.', linewidth=2, linestyle='None', label='Human/Gorilla')
plt.loglog(x13, y13, 'darkgray', marker='.', linewidth=2, linestyle='None', label='Human/Mouse')


plt.title('Indel-Terminated Segment Length Distribution Comparison Among\nEqually Distant Primates (28MYA) vs Human/Gorilla and Human/Mouse', fontsize=15)
plt.xlabel('Length (log10)', fontsize=20)
plt.ylabel('Number of L-mers (log10)', fontsize=20)
plt.tick_params(axis='x', labelsize=15)
plt.tick_params(axis='y', labelsize=15)
plt.legend(fontsize=10, markerscale=3, loc='lower left')

plt.savefig('/bucket/.mabuya/MillerU/Abrar/PCS/Uppercase_analysis/indel_terminated/indel_terminated_comparison_same_28MYA_primates.png')























