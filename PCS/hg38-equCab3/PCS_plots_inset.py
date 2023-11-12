import numpy as np
import matplotlib.pyplot as plt
import csv
import pandas as pd
import scipy
import lmfit

# General Functions
def linear(x, m, c):
    return c / (x ** (m))

df = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/PCS/hg38-equCab3/PCS_hg38equCab3.csv") 
df = df[df['Number of PCS of Length L']!=0]
# Preparing the data for the plot
x = [i for i in df.index][:]
x_fake = [i for i in df.index][11:500]

print(x_fake[0], x_fake[-1])
y = df['Number of PCS of Length L'][:]
y_fake = df['Number of PCS of Length L'][11:500]
x=np.array(x)
x_fake=np.array(x_fake)
y=np.array(y)

# Resizing the figure
plt.figure(figsize=[10, 7])

y_fit = linear(x_fake, 4, 10**10.8)

absError = np.log10(y_fit) - np.log10(y_fake)
SE = np.square(absError) # squared errors
MSE = np.mean(SE) # mean squared errors
RMSE = np.sqrt(MSE) # Root Mean Squared Error, RMSE
Rsquared = 1.0 - (np.var(absError) / np.var(np.log10(y_fake)))

# Plotting the graph with Log ticks at x and y axis using loglog
fig,ax=plt.subplots(figsize=[10,7])
axIns = ax.inset_axes([0.05, 0.05, 0.47, 0.47]) #.49, 0.48, 0.5, 0.5
axIns.semilogy(x , y, '*y', linewidth=2, label='Semi-Log L-mers')
axIns.semilogy(x_fake, y_fit, 'b--', label="Fitted Curve:\nSlope = -4,\nLower L cutoff = 12")
axIns.legend(fontsize=13, markerscale=3, loc='upper right')

ax.loglog(x, y, '.r', linewidth=2, label='Log-Log L-mers')
ax.loglog(x_fake, y_fit, 'b--', label="Fitted Curve:\nSlope = -4,\nLower L cutoff = 12")

ax.set_title('PCS length distribution of hg38/equCab3 genome alignment for L>=1', fontsize=15)
ax.annotate("MSE: {:.3f}\nRMSE: {:.3f}\nR-squared: {:.3f}".format(MSE,RMSE,Rsquared),(2,10**4.5), fontsize=15, xytext=(2,10**4.5))
ax.set_xlabel('Length (log10)', fontsize=20)
ax.set_ylabel('Number of L-mers (log10)', fontsize=20)
ax.tick_params(axis='x', labelsize=15)
ax.tick_params(axis='y', labelsize=15)

ax.legend(fontsize=15, markerscale=3, loc='upper right')

plt.savefig("/bucket/.mabuya/MillerU/Abrar/PCS/hg38-equCab3/hgequCab3.png")
#plt.show()
