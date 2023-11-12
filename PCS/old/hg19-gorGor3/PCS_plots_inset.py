import numpy as np
import matplotlib.pyplot as plt
import csv
import pandas as pd
import scipy
#pip3 install lmfit
import lmfit

df = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/PCS/hg19-gorGor3/PCS_hg19gorGor3.csv") 
df = df[df['Number of PCS of Length L']!=0][6:]
# Preparing the data for the plot
x = [i for i in df.index][:]
y = df['Number of PCS of Length L'][:]
x_fake = [i for i in df.index][:1010]
print(x_fake[0],x_fake[-1])
y_fake = df['Number of PCS of Length L'][:1010]
y_fake_2 = np.log10(y_fake)

# Resizing the figure
plt.figure(figsize=[10, 7])

logx, logy = np.log(x_fake), np.log(y_fake)

p = np.polyfit(logx, logy, 3)
y_fit = np.exp(np.polyval(p, logx)) #plyval->Compute polynomial values.

absError = np.log10(y_fit) - np.log10(y_fake)
SE = np.square(absError) # squared errors
MSE = np.mean(SE) # mean squared errors
RMSE = np.sqrt(MSE) # Root Mean Squared Error, RMSE
Rsquared = 1.0 - (np.var(absError) / np.var(np.log10(y_fake)))

# Plotting the graph with Log ticks at x and y axis using loglog
fig,ax=plt.subplots(figsize=[10,7])
axIns = ax.inset_axes([0.07, 0.06, 0.5, 0.55])
axIns.semilogy(x , y, '*y', linewidth=2, label='Semi-Log L-mers')
axIns.semilogy(x_fake, y_fit, 'b--', label="Fitted Curve:\nDegree = 3,\nLower L cutoff = 7")
axIns.legend(fontsize=13, markerscale=2, loc='upper right')

ax.loglog(x, y, '.r', linewidth=2, label='Log-Log L-mers')
ax.loglog(x_fake, y_fit, 'b--', label="Fitted Curve:\nDegree = 3,\nLower L cutoff = 7")

ax.set_title('PCS length distribution of hg19/gorGor3 genome alignment for L>=7', fontsize=15)
ax.annotate("MSE: {:.3f}\nRMSE: {:.3f}\nR-squared: {:.3f}".format(MSE,RMSE,Rsquared),(1.2*10**1,10**4), fontsize=15)
ax.set_xlabel('Length (log10)', fontsize=20)
ax.set_ylabel('Number of L-mers (log10)', fontsize=20)
ax.tick_params(axis='x', labelsize=15)
ax.tick_params(axis='y', labelsize=15)
ax.legend(fontsize=15, markerscale=3, loc='upper right')

plt.savefig("/bucket/MillerU/Abrar/PCS/hg19-gorGor3/hggg.png")
#ax.show()