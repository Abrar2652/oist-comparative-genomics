

# importing necessary libraries 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats
import statsmodels.api as sm

sns.set_theme()
sns.set_palette(palette = "rainbow")
%matplotlib inline


df1 = pd.read_csv("/content/gorilla.csv").drop(['index'],axis=1).rename(columns={'0': 'Human/Gorilla'})
df2 = pd.read_csv("/content/galvar.csv").drop(['index'],axis=1).rename(columns={'0': 'Human/Malayan Flying Lemur'}) 
df3 = pd.read_csv("/content/micmur.csv").drop(['index'],axis=1).rename(columns={'0': 'Human/Mouse Lemur'}) 
df4 = pd.read_csv("/content/mouse.csv").drop(['index'],axis=1).rename(columns={'0': 'Human/Mouse'}) 
df5 = pd.read_csv("/content/orangutan.csv").drop(['index'],axis=1).rename(columns={'0': 'Human/Orangutan'}) 
df6 = pd.read_csv("/content/nomleu.csv").drop(['index'],axis=1).rename(columns={'0': 'Human/Gibbon'}) 
df7 = pd.read_csv("/content/rhesus.csv").drop(['index'],axis=1).rename(columns={'0': 'Human/Rhesus'}) 
df8 = pd.read_csv("/content/saibol.csv").drop(['index'],axis=1).rename(columns={'0': 'Human/Squirrel Monkey'})


df = pd.concat([df1,df2,df3,df4,df5,df6,df7,df8], axis=1)

df_table_1 = df.describe().transpose()


print(df_table_1)


# Skew function of Pandas
old_skew = df.skew().sort_values(ascending=False)
print(old_skew)

##### Q-Q PLOTS ##########

plt.figure(figsize=(17,9))
for i in list(enumerate(df.columns)):
    plt.subplot(2, 4,i[0]+1)
    stats.probplot(df[i[1]], dist="norm", plot=plt)   # QQ Plot
    plt.title("{}".format(i[1]), fontsize = 15, weight='bold')
    plt.xticks(fontsize = 12, weight='bold')
    plt.yticks(fontsize = 12, weight='bold')
    plt.xlabel('Theoretical Quantities', fontsize = 12, weight='bold')
    plt.ylabel('Sample Quantities', fontsize = 12, weight='bold')
plt.tight_layout()  
plt.savefig('Q-Q_plot_skewness.png')

####### KDE PLOTS ########

plt.figure(figsize=(25,9))
for i in list(enumerate(df.columns)):
    plt.subplot(2, 4, i[0]+1)
    sns.histplot(data = df[i[1]], kde=True)  # Histogram with KDE line
  
#plt.tight_layout()  
#plt.show()
plt.savefig('Histogram_KDE_plot_skewness.png')






















