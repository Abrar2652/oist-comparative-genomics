# Importing necessary libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Preparing the data for the plot
genome_list1 = ['hg38']
names_list1 = ['Human']
genome_list2 = ['gorGor6']
names_list2 = ['Gorilla']

for genome1,name1,genome2,name2 in zip(genome_list1,names_list1,genome_list2,names_list2):
    df1 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/PCS/AR/indel_terminated/vanilla_indel_terminated_PCS_{0}{1}.csv".format(genome1,genome2)) # indel-terminated

    df2 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/PCS/AR/both_mismatch_indel_terminated/GA_CT/GA_CT_both_mismatch_indel_PCS_{0}{1}.csv".format(genome1,genome2)) # both-mismatch-indel-terminated
    df3 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/PCS/AR/both_mismatch_indel_terminated/GC_AT/GC_AT_both_mismatch_indel_PCS_{0}{1}.csv".format(genome1,genome2)) # both-mismatch-indel-terminated
    df4 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/PCS/AR/both_mismatch_indel_terminated/GT_AC/GT_AC_both_mismatch_indel_PCS_{0}{1}.csv".format(genome1,genome2)) # both-mismatch-indel-terminated
    df5 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/PCS/{0}-{1}/PCS_{0}{1}.csv".format(genome1,genome2)) # both-mismatch-indel-terminated            
    
    df6 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/PCS/AR/mismatch_terminated/GA_CT/GA_CT_mismatch_terminated_PCS_{0}{1}.csv".format(genome1,genome2)) # mismatch_terminated
    df7 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/PCS/AR/mismatch_terminated/GC_AT/GC_AT_mismatch_terminated_PCS_{0}{1}.csv".format(genome1,genome2)) # mismatch_terminated
    df8 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/PCS/AR/mismatch_terminated/GT_AC/GT_AC_mismatch_terminated_PCS_{0}{1}.csv".format(genome1,genome2)) # mismatch_terminated
    df9 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/PCS/AR/mismatch_terminated/vanilla/vanilla_mismatch_terminated_PCS_{0}{1}.csv".format(genome1,genome2)) # mismatch_terminated
        
    
    # Importing necessary libraries
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    
    # Preparing the data for the plot
    
    x1 = [i for i in df1.index][1:]    
    x2 = [i for i in df2.index][1:]
    x3 = [i for i in df3.index][1:]
    x4 = [i for i in df4.index][1:]
    x5 = [i for i in df5.index][1:]
    x6 = [i for i in df6.index][1:]
    x7 = [i for i in df7.index][1:]
    x8 = [i for i in df8.index][1:]
    x9 = [i for i in df9.index][1:]
    y1 = df1['Number of PCS of Length L'][1:]
    y2 = df2['Number of PCS of Length L'][1:]
    y3 = df3['Number of PCS of Length L'][1:]
    y4 = df4['Number of PCS of Length L'][1:]
    y5 = df5['Number of PCS of Length L'][1:]
    y6 = df6['Number of PCS of Length L'][1:]
    y7 = df7['Number of PCS of Length L'][1:]
    y8 = df8['Number of PCS of Length L'][1:]
    y9 = df9['Number of PCS of Length L'][1:]
        
    
    # Resizing the figure
    plt.figure(figsize=[10, 7])
    
    # Plotting the graph with Log ticks at x and y axis using loglog
    plt.loglog(x1, y1, 'red', marker='.', linewidth=2, linestyle='None', label='Indel terminated PCS')
    plt.loglog(x2, y2, 'blue', marker='.', linewidth=2, linestyle='None', label='G->A, C->T both mismatch and indel terminated PCS')
    plt.loglog(x3, y3, 'green', marker='.', linewidth=2, linestyle='None', label='G->C, A->T both mismatch and indel terminated PCS')
    plt.loglog(x4, y4, 'gold', marker='.', linewidth=2, linestyle='None', label='G->T, A->C both mismatch and indel terminated PCS')
    plt.loglog(x5, y5, 'olive', marker='.', linewidth=2, linestyle='None', label='Vanilla both mismatch and indel terminated PCS')
    plt.loglog(x6, y6, 'firebrick', marker='.', linewidth=2, linestyle='None', label='G->A, C->T mismatch terminated PCS')
    plt.loglog(x7, y7, 'darkorange', marker='.', linewidth=2, linestyle='None', label='G->C, A->T mismatch terminated PCS')
    plt.loglog(x8, y8, 'crimson', marker='.', linewidth=2, linestyle='None', label='G->T, A->C mismatch terminated PCS')
    plt.loglog(x9, y9, 'darkslategray', marker='.', linewidth=2, linestyle='None', label='Vanilla mismatch terminated PCS')
    
    
    #plt.title('All Indel Terminated PCS Length Distribution Combinations\nFrom {1}/{0} vs hg38/mm39 Whole Genome Alignments (log-log)'.format(genome2,genome1), fontsize=14)
    plt.title('All Possible Combinations of Vanilla PCS Length Distribution\nFrom {1}/{0} Whole Genome Alignment (log-log)'.format(genome2,genome1), fontsize=14)
    #plt.title('Both Mismatch And Indel Terminated PCS Length Distribution\nCombinations From {1}/{0} Whole Genome Alignment (log-log)'.format(genome2,genome1), fontsize=14)

    plt.xlabel('Length (log10)', fontsize=20)
    plt.ylabel('Number of L-mers (log10)', fontsize=20)
    plt.tick_params(axis='x', labelsize=15)
    plt.tick_params(axis='y', labelsize=15)
    
    #plt.legend(fontsize=11, markerscale=3, loc='upper right')
    plt.legend(fontsize=10, markerscale=2, loc='upper right')
    plt.savefig('/bucket/.mabuya/MillerU/Abrar/PCS/AR/all_possible_combinations_PCS_{1}{0}_loglog.png'.format(genome2,genome1))

    


    




