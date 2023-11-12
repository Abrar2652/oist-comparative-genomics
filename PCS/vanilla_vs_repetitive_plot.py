# Importing necessary libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Preparing the data for the plot
genome_list = ['calJac4','danRer10','equCab3','felCat9','galGal6','ornAna2','oryCun2','oviAri4','petMar3','susScr11','tarSyr2','thaSir1','triMan1','xenTro10']
names_list = ['Marmoset','Tarsier','Rabbit','Horse','Cat','Sheep','Pig','Manatee','Platypus','Chicken','Garter Snake','X. tropicalis','Zebrafish','Lamprey']
for genome,name in zip(genome_list, names_list):

    df1 = pd.read_csv("/bucket/MillerU/Abrar/PCS/hg38-{0}/PCS_hg38{0}_small_to_upper.csv".format(genome)) #repetitive+vanilla PCS
    df2 = pd.read_csv("/bucket/MillerU/Abrar/PCS/hg38-{0}/PCS_hg38{0}.csv".format(genome)) #vanilla PCS
    
    # Preparing the data for the plot
    x1 = [i for i in df1.index][1:]
    y1 = df1['Number of PCS of Length L'][1:]
    x2 = [i for i in df2.index][1:]
    y2 = df2['Number of PCS of Length L'][1:]
    
    # Resizing the figure
    plt.figure(figsize=[10, 7])
    
    # Plotting the graph with Log ticks at x and y axis using loglog
    plt.loglog(x1, y1, '.r', linewidth=2, label='Human/{0} Vanilla + Repetitive L-mers'.format(name))
    plt.loglog(x2, y2, '.b', linewidth=2, label='Human/{0} Vanilla L-mers'.format(name))
    plt.title('PCS Length Distribution Comparison (Vanilla vs Vanilla + Repetitive)\n of hg38/{0} genome alignment (log-log) for L>=1'.format(genome), fontsize=14)
    plt.xlabel('Length (log10)', fontsize=20)
    plt.ylabel('Number of L-mers (log10)', fontsize=20)
    plt.tick_params(axis='x', labelsize=15)
    plt.tick_params(axis='y', labelsize=15)
    
    plt.legend(fontsize=13, markerscale=3)
    plt.savefig('/bucket/MillerU/Abrar/PCS/hg38-{0}/PCS_hg38{0}_repetitive_vs_vanilla_loglog.png'.format(genome))








