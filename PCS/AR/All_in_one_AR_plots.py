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
######## indel terminated ########
    df1 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/PCS/AR/indel_terminated/AR_both_lowercase_match/indel_terminated_AR_{0}{1}_both_lowercase_must.csv".format(genome1,genome2)) # indel-terminated
    df2 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/PCS/AR/indel_terminated/AR_only_hg_lowercase_match_LPH/indel_terminated_AR_{0}{1}_hg_lowercase_only.csv".format(genome1,genome2)) # indel-terminated
    df3 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/PCS/AR/indel_terminated/AR_only_hg_lowercase_ignore_{1}/indel_terminated_AR_{0}{1}_hg_lowercase_ignore_{1}.csv".format(genome1,genome2)) # indel-terminated
    
    ######## mismatch terminated ######## 
    df4 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/PCS/AR/mismatch_terminated/AR_exact_match/AR_mismatch_terminated_exact_match_PCS_{0}{1}.csv".format(genome1,genome2)) 
    df5 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/PCS/AR/mismatch_terminated/AR_ga_ct/AR_ga_ct_mismatch_terminated_{0}{1}.csv".format(genome1,genome2)) 
    df6 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/PCS/AR/mismatch_terminated/AR_gc_at/AR_gc_at_mismatch_terminated_{0}{1}.csv".format(genome1,genome2)) 
    df7 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/PCS/AR/mismatch_terminated/AR_gt_ac/AR_gt_ac_mismatch_terminated_{0}{1}.csv".format(genome1,genome2))          
    
    ######## both mismatch indel terminated ########
    df8 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/PCS/AR/both_mismatch_indel_terminated/AR_exact_match/AR_exact_match_both_mismatch_indel_PCS_{0}{1}.csv".format(genome1,genome2)) 
    df9 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/PCS/AR/both_mismatch_indel_terminated/AR_ga_ct/AR_ga_ct_both_mismatch_indel_PCS_{0}{1}.csv".format(genome1,genome2))
    df10 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/PCS/AR/both_mismatch_indel_terminated/AR_gc_at/AR_gc_at_both_mismatch_indel_PCS_{0}{1}.csv".format(genome1,genome2)) 
    df11 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/PCS/AR/both_mismatch_indel_terminated/AR_gt_ac/AR_gt_ac_both_mismatch_indel_PCS_{0}{1}.csv".format(genome1,genome2))

    ##### CMMS ########          
    df12 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/PCS/AR/CMMS_match_terminated/AR_CMMS_match_terminated_PCS_{0}{1}.csv".format(genome1,genome2)) # match terminated
    df13 = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/PCS/AR/CMMS_both_match_indel_terminated/AR_CMMS_both_match_indel_terminated_PCS_{0}{1}.csv".format(genome1,genome2))  # both match indel terminated
    
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
    x10 = [i for i in df10.index][1:]
    x11 = [i for i in df11.index][1:]
    x12 = [i for i in df12.index][1:]
    x13 = [i for i in df13.index][1:]
    y1 = df1['Number of PCS of Length L'][1:]
    y2 = df2['Number of PCS of Length L'][1:]
    y3 = df3['Number of PCS of Length L'][1:]
    y4 = df4['Number of PCS of Length L'][1:]
    y5 = df5['Number of PCS of Length L'][1:]
    y6 = df6['Number of PCS of Length L'][1:]
    y7 = df7['Number of PCS of Length L'][1:]
    y8 = df8['Number of PCS of Length L'][1:]
    y9 = df9['Number of PCS of Length L'][1:]
    y10 = df10['Number of PCS of Length L'][1:]
    y11 = df11['Number of PCS of Length L'][1:]
    y12 = df12['Number of PCS of Length L'][1:]
    y13 = df13['Number of PCS of Length L'][1:]
        
    
    # Resizing the figure
    plt.figure(figsize=[10, 7])
    
    # Plotting the graph with Log ticks at x and y axis using loglog
    plt.loglog(x1, y1, 'red', marker='.', linewidth=2, linestyle='None', label='Indel terminated AR (both genomes lowercase)')
    plt.loglog(x2, y2, 'purple', marker='.', linewidth=2, linestyle='None', label='Indel terminated AR (only hg38 lowercase)')
    plt.loglog(x3, y3, 'black', marker='.', linewidth=2, linestyle='None', label='Only hg38 genome indel terminated AR')

    plt.loglog(x4, y4, 'darkslategray', marker='.', linewidth=2, linestyle='None', label='Exactly matched (mismatch terminated) AR')
    plt.loglog(x5, y5, 'firebrick', marker='.', linewidth=2, linestyle='None', label='g->a, c->t (mismatch terminated) AR')
    plt.loglog(x6, y6, 'darkorange', marker='.', linewidth=2, linestyle='None', label='g->c, a->t (mismatch terminated) AR')
    plt.loglog(x7, y7, 'crimson', marker='.', linewidth=2, linestyle='None', label='g->t, a->c (mismatch terminated) AR')
        
    plt.loglog(x8, y8, 'olive', marker='.', linewidth=2, linestyle='None', label='Exactly matched (both mismatch+indel\nterminated) AR')
    plt.loglog(x9, y9, 'blue', marker='.', linewidth=2, linestyle='None', label='g->a, c->t (both mismatch+indel terminated) AR')
    plt.loglog(x10, y10, 'green', marker='.', linewidth=2, linestyle='None', label='g->c, a->t (both mismatch+indel terminated) AR')
    plt.loglog(x11, y11, 'gold', marker='.', linewidth=2, linestyle='None', label='g->t, a->c (both mismatch+indel terminated) AR')
    
    plt.loglog(x12, y12, 'darkgray', marker='.', linewidth=2, linestyle='None', label='CMMS (match terminated) AR')
    plt.loglog(x13, y13, 'teal', marker='.', linewidth=2, linestyle='None', label='CMMS (both match+indel terminated) AR')
    
    
    
    #plt.title('All Indel Terminated AR Length Distribution Combinations\nFrom {1}/{0} vs hg38/mm39 Whole Genome Alignments (log-log)'.format(genome2,genome1), fontsize=14)
    plt.title('All Possible Combinations of AR Length Distributions\nFrom {1}/{0} Whole Genome Alignment (log-log)'.format(genome2,genome1), fontsize=14)
    #plt.title('Both Mismatch And Indel Terminated AR Length Distribution\nCombinations From {1}/{0} Whole Genome Alignment (log-log)'.format(genome2,genome1), fontsize=14)

    plt.xlabel('Length (log10)', fontsize=20)
    plt.ylabel('Number of AR L-mers (log10)', fontsize=20)
    plt.tick_params(axis='x', labelsize=15)
    plt.tick_params(axis='y', labelsize=15)
    
    #plt.legend(fontsize=11, markerscale=3, loc='upper right')
    plt.legend(fontsize=8, markerscale=2, loc='upper right')
    plt.savefig('/bucket/.mabuya/MillerU/Abrar/PCS/AR/all_possible_combinations_AR_{1}{0}_loglog.png'.format(genome2,genome1))

    


    




