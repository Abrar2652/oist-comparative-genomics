import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd

df = pd.read_excel('/content/displacement_vs_divergence.xlsx')
plt.figure(figsize=[10, 7])

sp = df['Query species'].tolist()
# Plotting the graph 
plt.xlim([0,100])
plt.ylim([0,350])
plt.plot(df['Divergence time (MYA)'], df['Displacement'], color='red', marker='.', markersize=20, linewidth=2, linestyle='-')

for i, txt in enumerate(sp):
  plt.annotate(txt, (df['Divergence time (MYA)'].tolist()[i], df['Displacement'].tolist()[i]), ha='center', va='bottom', textcoords='data', size=13,color="Green")
plt.title('PCS Power-law Curve Displacement from x=0 line along x-axis\n Vs Divergence Time', fontsize=15)
plt.ylabel('Displacement from x=0 line', fontsize=20)
plt.xlabel('Divergence Time (MYA)', fontsize=20)
plt.tick_params(axis='x', labelsize=15)
plt.tick_params(axis='y', labelsize=15)

green_patch = mpatches.Patch(color='green', label='Query Species')
plt.legend(handles=[green_patch], fontsize=15, markerscale=3)
plt.savefig('PCS_displacement_vs_divergence.png')
