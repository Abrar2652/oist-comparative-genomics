# Importing necessary libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import csv
import scipy
import lmfit

# Preparing the data for the plot
# General Functions
def linear(x, m, c):
    return c / (x ** (m))


genome_list = ['mm39']
names_list = ['Mouse']
for genome,name in zip(genome_list, names_list):

    df1 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/PCS/PCS_common_intersection_primate_L20.csv") #repetitive+vanilla PCS
    df2 = pd.read_csv("/bucket/MillerU/Abrar/PCS/hg38-{0}/PCS_hg38{0}.csv".format(genome)) #vanilla PCS
    
    # Preparing the data for the plot
    x1 = [i for i in df1.index][20:]
    x1_fake = [i for i in df1.index][20:190]
    x1 = np.array(x1)
    x1_fake = np.array(x1_fake)
    y1 = df1['Number of PCS of Length L'][20:]
    y1 = np.array(y1)
    y1_fake = df1['Number of PCS of Length L'][20:190]
    x2 = [i for i in df2.index][20:]
    y2 = df2['Number of PCS of Length L'][20:]
    
    # Resizing the figure 
    plt.figure(figsize=[10, 7])
    
    y1_fit = linear(x1_fake, 4.5, 10**10.2)
    
    absError = np.log10(y1_fit) - np.log10(y1_fake)
    SE = np.square(absError) # squared errors
    MSE = np.mean(SE) # mean squared errors
    RMSE = np.sqrt(MSE) # Root Mean Squared Error, RMSE
    Rsquared = 1.0 - (np.var(absError) / np.var(np.log10(y1_fake)))

    # Plotting the graph with Log ticks at x and y axis using loglog
    plt.loglog(x1, y1, '.r', linewidth=2, label='Exactly Matched L-mers Among Primates')
    plt.loglog(x1_fake, y1_fit, 'b--', label="Fitted Curve:\nSlope = -4.5,\nLower L cutoff = 20")
    plt.loglog(x2, y2, '.g', linewidth=2, label='Human/{0} Vanilla L-mers'.format(name))
    plt.title('PCS Length Distribution Comparison Between Exactly Matched PCS Among Primates vs\n hg38/{0} Vanilla PCS (log-log) for L>=20'.format(genome), fontsize=14)
    plt.xlabel('Length (log10)', fontsize=20)
    plt.ylabel('Number of L-mers (log10)', fontsize=20)
    plt.tick_params(axis='x', labelsize=15)
    plt.tick_params(axis='y', labelsize=15)
    
    plt.legend(fontsize=13, markerscale=3)
    plt.savefig('/bucket/MillerU/Abrar/PCS/common_primates_pcs_vs_hg38{0}_loglog.png'.format(genome))








