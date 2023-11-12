# Importing necessary libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Preparing the data for the plot

df = pd.read_csv("/bucket/MillerU/Abrar/PCS/hg19-gorGor3/PCS_hg19gorGor3_small_to_upper.csv")
print(df.head())
x = [i for i in df.index][1:]
y = df['Number of PCS of Length L'][1:]

# Resizing the figure
plt.figure(figsize=[10, 7])

# Plotting the graph with Log ticks at x and y axis using loglog
plt.plot(x, y, '.r', linewidth=2, 
           label='Human/Gorilla whole genome L-mers')
plt.title('Length distributions of Perfectly Conserved\n\
 L-mers (with repetitive sequences) from hg19/gorGor3 genome alignment for L>=1', 
          fontsize=14)
plt.xlabel('Length', fontsize=20)
plt.ylabel('Number of L-mers', fontsize=20)
plt.tick_params(axis='x', labelsize=15)
plt.tick_params(axis='y', labelsize=15)

plt.legend(fontsize=13, markerscale=3)
plt.savefig('/bucket/MillerU/Abrar/PCS/hg19-gorGor3/PCS_hg19gorGor3_small_to_upper_linear.png')



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
plt.loglog(x, y, '.r', linewidth=2, 
           label='Human/Gorilla whole genome L-mers')
plt.title('Length distributions of Perfectly Conserved\n\
 L-mers (with repetitive sequences) from hg19/gorGor3 genome alignment (log-log) for L>=1', 
          fontsize=14)
plt.xlabel('Length (log10)', fontsize=20)
plt.ylabel('Number of L-mers (log10)', fontsize=20)
plt.tick_params(axis='x', labelsize=15)
plt.tick_params(axis='y', labelsize=15)

plt.legend(fontsize=13, markerscale=3)
plt.savefig('/bucket/MillerU/Abrar/PCS/hg19-gorGor3/PCS_hg19gorGor3_small_to_upper_loglog_L1.png')



# Preparing the data for the plot
x = [i for i in df.index][1:]
y = df['Number of PCS of Length L'][1:]

plt.figure(figsize=[10, 7])

plt.semilogy(x , y, '*m', linewidth=2, base=10, label='Human/Gorilla whole genome L-mers')
plt.xlabel("Length (L)", fontsize=20)
plt.ylabel('Number of L-mers (log10)', fontsize=20)
plt.tick_params(axis='x', labelsize=15)
plt.tick_params(axis='y', labelsize=15)

plt.title("Length distributions of Perfectly Conserved \n\
L-mers (with repetitive sequences) from hg19/gorGor3 genome alignment (semi-log) for L>=1",fontsize=14)

plt.legend(fontsize=13, markerscale=3)
plt.savefig('/bucket/MillerU/Abrar/PCS/hg19-gorGor3/PCS_hg19gorGor3_small_to_upper_semilog.png')






