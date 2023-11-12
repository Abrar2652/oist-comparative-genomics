# Importing necessary libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Preparing the data for the plot
genome_list1 = ['hg38']
names_list1 = ['Human']
genome_list2 = ['mm39']
names_list2 = ['Mouse']

for genome1,name1,genome2,name2 in zip(genome_list1,names_list1,genome_list2,names_list2):
    #df = pd.read_csv("/bucket/MillerU/Abrar/PCS/hg38-{0}/PCS_hg38{0}.csv".format(genome2))
    df = pd.read_csv("/bucket/MillerU/Abrar/PCS/{1}-{0}/PCS_{1}{0}.csv".format(genome2,genome1))
    #df1 = pd.read_csv("/bucket/MillerU/Abrar/PCS/hg19-mm10/PCS_hg19mm10.csv")
    
    '''
    x = [i for i in df.index][1:]
    y = df['Number of PCS of Length L'][1:]
    
    # Resizing the figure
    plt.figure(figsize=[10, 7])
    
    # Plotting the graph with Log ticks at x and y axis using loglog
    plt.plot(x, y, '.r', linewidth=2, 
               label='Human/{0} whole genome L-mers'.format(name2))
    plt.title('Length distributions of Perfectly Conserved\n\L-mers from hg38/{0} genome alignment for L>=1'.format(genome2), fontsize=14)
    plt.xlabel('Length', fontsize=20)
    plt.ylabel('Number of L-mers', fontsize=20)
    plt.tick_params(axis='x', labelsize=15)
    plt.tick_params(axis='y', labelsize=15)
    
    plt.legend(fontsize=15, markerscale=3)
    plt.savefig('/bucket/MillerU/Abrar/PCS/hg38-{0}/PCS_hg38{0}_linear.png'.format(genome2))
    '''
    
    
    # Importing necessary libraries
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    
    # Preparing the data for the plot
    
    x = [i for i in df.index][1:]
    y = df['Number of PCS of Length L'][1:]
    #x1 = [i for i in df1.index][1:]
    #y1 = df1['Number of PCS of Length L'][1:]
        
    # Resizing the figure
    plt.figure(figsize=[10, 7])
    
    # Plotting the graph with Log ticks at x and y axis using loglog
    #plt.loglog(x, y, '.r', linewidth=2, label='{1}/{0} whole genome L-mers'.format(name2,name1))
    plt.loglog(x, y, '.r', linewidth=2, label='{1}/{0} chrY genome L-mers'.format(name2,name1))
    plt.title('Length Distributions of Perfectly Conserved\nL-mers From {1}/{0} chrY Genome Alignment (log-log) For L>=1'.format(genome2,genome1), fontsize=14)
    #plt.title('PCS Length Distribution Comparison of hg19/mm10 And\ngorGor6/mm10 (log-log) For L>=1'.format(genome), fontsize=14)
    plt.xlabel('Length (log10)', fontsize=20)
    plt.ylabel('Number of L-mers (log10)', fontsize=20)
    plt.tick_params(axis='x', labelsize=15)
    plt.tick_params(axis='y', labelsize=15)
    
    plt.legend(fontsize=15, markerscale=3)
    plt.savefig('/bucket/MillerU/Abrar/PCS/{1}-{0}/chrY_PCS_{1}{0}_loglog.png'.format(genome2,genome1))
    #plt.savefig('/bucket/MillerU/Abrar/PCS/gorGor6-{0}/PCS_hgmm_ggmm_loglog.png'.format(genome2))
    
    '''
    
    # Preparing the data for the plot
    x = [i for i in df.index][1:]
    y = df['Number of PCS of Length L'][1:]
    
    plt.figure(figsize=[10, 7])
    
    plt.semilogy(x , y, '*m', linewidth=2, base=10, label='Human/{0} whole genome L-mers'.format(name2))
    plt.xlabel("Length (L)", fontsize=20)
    plt.ylabel('Number of L-mers (log10)', fontsize=20)
    plt.tick_params(axis='x', labelsize=15)
    plt.tick_params(axis='y', labelsize=15)
    
    
    plt.title("Length distributions of Perfectly Conserved\n\L-mers from hg38/{0} genome alignment (semi-log) for L>=1".format(genome2), fontsize=14)
    
    plt.legend(fontsize=15, markerscale=3)
    
    plt.savefig('/bucket/MillerU/Abrar/PCS/hg38-{0}/PCS_hg38{0}_semilog.png'.format(genome2))
    
    '''




