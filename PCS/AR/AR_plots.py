# Importing necessary libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Preparing the data for the plot
genome_list1 = ['hg38']
names_list1 = ['Human']
genome_list2 = ['mm39']
names_list2 = ['Mouse']

'''
## coversion to csv file from existing BED file having coordinates

df1 = pd.read_csv("/bucket/.mabuya/MillerU/Zifcakova/genome_coordinates2/AR/hg19_mm10_AR_positions_in_human.bed", header=None, sep='\t')
df1 = pd.read_csv("/bucket/.mabuya/MillerU/Zifcakova/genome_coordinates2/Lunter/Lunter.bed", header=None, sep='\t')

list1 = df1[1].values.tolist()
list2 = df1[2].values.tolist()

print(list1[:10])
print(list2[:10])

seq_count = list()

for i in range(len(list1)):
  item = list2[i] - list1[i] + 1
  seq_count.append(item)

serial_count=[]
for i in range(0, int(max(seq_count)+5)):
    serial_count.append(seq_count.count(i))
   
lmer = [x for x in range(int(max(seq_count)+5))]

serial_count=[]
for i in range(0, 500):
    serial_count.append(seq_count.count(i))
   
lmer = [x for x in range(500)]
# dictionary of lists  
dict = {'PCS length (L)': lmer, 'Number of PCS of Length L': serial_count}  
       
df = pd.DataFrame(dict)  
df.to_csv(r'/bucket/MillerU/Abrar/PCS/AR/Lunter_IGS_length_counts.csv', index = False) 
'''

for genome1,name1,genome2,name2 in zip(genome_list1,names_list1,genome_list2,names_list2):
    sub = 'Contiguosly Mis-Matched'
    sub_name = 'CMMS_both_match_indel'
    path = "/bucket/.mabuya/MillerU/Abrar/PCS/AR/CMMS_both_match_indel_terminated/"
    #df = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/PCS/AR/indel_terminated/indel_terminated_PCS_{0}{1}.csv".format(genome1,genome2)) # indel-terminated
    df = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/PCS/AR/indel_terminated/AR_both_lowercase_match/indel_terminated_AR_{0}{1}_both_lowercase_must.csv".format(genome1,genome2,sub_name)) # both-mismatch-indel-terminated
    #df = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/PCS/AR/mismatch_terminated/GC_AT/GC_AT_mismatch_terminated_PCS_{0}{1}.csv".format(genome1,genome2)) # mismatch_terminated
    #df = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/PCS/AR/indel_terminated/AR_only_hg_lowercase_match/indel_terminated_AR_{0}{1}_hg_lowercase_only.csv".format(genome1,genome2))
    #df = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/PCS/AR/mismatch_terminated/AR_exact_match/AR_mismatch_terminated_exact_match_PCS_hg38gorGor6.csv")
    
  
    x = [i for i in df.index][1:]
    y = df['Number of PCS of Length L'][1:]
    
    # Resizing the figure
    plt.figure(figsize=[10, 7])
    
    # Plotting the graph with Log ticks at x and y axis using loglog
    plt.plot(x, y, '.r', linewidth=2, label='Human/{0} {1}\n AR length counts'.format(name2,sub))
    #plt.title('Indel Terminated AR Length Distribution of {0}/{1} (linear)'.format(genome1,genome2), fontsize=14)
    plt.title("{2} (both match + indel terminated) AR Length\n Distributions of {0}/{1} (linear)".format(genome1,genome2,sub), fontsize=14)
    plt.xlabel('Length', fontsize=20)
    plt.ylabel('Number of AR L-mers', fontsize=20)
    plt.tick_params(axis='x', labelsize=15)
    plt.tick_params(axis='y', labelsize=15)
    
    plt.legend(fontsize=13, markerscale=3)
    plt.savefig(path+'/AR_{2}_terminated_{0}{1}_linear.png'.format(genome1,genome2,sub_name))

    
    
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
    plt.loglog(x, y, '.r', linewidth=2, label='{1}/{0} {2}\n AR length counts'.format(name2,name1,sub))
    #plt.title('Indel Terminated AR Length Distribution of {0}/{1} (log-log)'.format(genome1,genome2), fontsize=14)
    plt.title("{2} (both match + indel terminated) AR Length\n Distributions of {0}/{1} (loglog)".format(genome1,genome2,sub), fontsize=14)
    plt.xlabel('Length (log10)', fontsize=20)
    plt.ylabel('Number of AR L-mers (log10)', fontsize=20)
    plt.tick_params(axis='x', labelsize=15)
    plt.tick_params(axis='y', labelsize=15)
    
    plt.legend(fontsize=13, markerscale=3)
    plt.savefig(path+'/AR_{2}_terminated_{0}{1}_loglog.png'.format(genome1,genome2,sub_name))


    # Preparing the data for the plot
    x = [i for i in df.index][1:]
    y = df['Number of PCS of Length L'][1:]
    
    plt.figure(figsize=[10, 7])    
    plt.semilogy(x , y, '*m', linewidth=2, base=10, label='Human/{0} {1}\n AR length counts'.format(name2,sub))
    plt.xlabel("Length (L)", fontsize=20)
    plt.ylabel('Number of AR L-mers (log10)', fontsize=20)
    plt.tick_params(axis='x', labelsize=15)
    plt.tick_params(axis='y', labelsize=15)        
    plt.title("{2} (both match + indel terminated) AR Length\n Distributions of {0}/{1} (semilog)".format(genome1,genome2,sub), fontsize=14)
    plt.legend(fontsize=13, markerscale=3)
    
    plt.savefig(path+'/AR_{2}_terminated_{0}{1}_semilog.png'.format(genome1,genome2,sub_name))
    

    




