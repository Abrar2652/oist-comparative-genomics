

import pandas as pd
import numpy as np

seq_count = pd.read_csv("/bucket/MillerU/Abrar/hg19-mm10_outputs/filtered_pcs_hg19-mm10_all/all/hg19_mm10_all_identical_seqs.tsv", sep='\t')['first.width'].values.tolist()
print(len(seq_count))
serial_count=[]
for i in range(0, int(max(seq_count)+5)):
    serial_count.append(seq_count.count(i))


lmer = [x for x in range(int(max(seq_count)+5))]
# dictionary of lists  
dict = {'PCS length (L)': lmer, 'Number of PCS of Length L': serial_count}  
       
df = pd.DataFrame(dict) 
    
# saving the dataframe 
df.to_csv(r'/bucket/MillerU/Abrar/PCS/hg19-mm10/Filtered_PCS_hg19mm10.csv', index = False) 



# Importing necessary libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Preparing the data for the plot

df = pd.read_csv("/bucket/MillerU/Abrar/PCS/hg19-mm10/Filtered_PCS_hg19mm10.csv")
x = [i for i in df.index][1:]
y = df['Number of PCS of Length L'][1:]

# Resizing the figure
plt.figure(figsize=[10, 7])

# Plotting the graph with Log ticks at x and y axis using loglog
plt.plot(x, y, '.r', linewidth=2, label='Human/Mouse whole genome L-mers')
plt.title('Length distributions of Perfectly Conserved\n\
 L-mers (filtered exons) from hg19/mm10 genome alignment for L>=1', fontsize=14)
plt.xlabel('Length', fontsize=20)
plt.ylabel('Number of L-mers', fontsize=20)
plt.tick_params(axis='x', labelsize=15)
plt.tick_params(axis='y', labelsize=15)

plt.legend(fontsize=15, markerscale=3)
plt.savefig('/bucket/MillerU/Abrar/PCS/hg19-mm10/Filtered_PCS_hg19mm10_linear_L1.png')



# Importing necessary libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Preparing the data for the plot

x = [i for i in df.index][1:]
y = df['Number of PCS of Length L'][1:]

# Resizing the figure
plt.figure(figsize=[10, 7])

# Plotting the graph with Log ticks at x and y axis using loglog
plt.loglog(x, y, '.r', linewidth=2, label='Human/Mouse whole genome L-mers')
plt.title('Length distributions of Perfectly Conserved\n\
 L-mers (filtered exons) from hg19/mm10 genome alignment (log-log) for L>=1', fontsize=14)
plt.xlabel('Length (log10)', fontsize=20)
plt.ylabel('Number of L-mers (log10)', fontsize=20)
plt.tick_params(axis='x', labelsize=15)
plt.tick_params(axis='y', labelsize=15)

plt.legend(fontsize=15, markerscale=3)
plt.savefig('/bucket/MillerU/Abrar/PCS/hg19-mm10/Filtered_PCS_hg19mm10_loglog_L1.png')




# Preparing the data for the plot
x = [i for i in df.index][1:]
y = df['Number of PCS of Length L'][1:]

plt.figure(figsize=[10, 7])

plt.semilogy(x , y, '*m', linewidth=2, base=10, label='Human/Mouse whole genome L-mers')
plt.xlabel("Length (L)", fontsize=20)
plt.ylabel('Number of L-mers (log10)', fontsize=20)
plt.tick_params(axis='x', labelsize=15)
plt.tick_params(axis='y', labelsize=15)

plt.title("Length distributions of Perfectly Conserved\n\
 L-mers (filtered exons) from hg19/mm10 genome alignment (semi-log) for L>=1",fontsize=14)

plt.legend(fontsize=15, markerscale=3)
plt.savefig('/bucket/MillerU/Abrar/PCS/hg19-mm10/Filtered_PCS_hg19mm10_semilog_L1.png')













