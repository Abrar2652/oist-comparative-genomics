
### Count the total PCS bases out-of-bound AXT 
import pandas as pd
df2 = pd.read_csv("/bucket/MillerU/Zifcakova/genome_coordinates2/hg19_gorGor3_new/pcs_axt_hg19_gorGor3_intersect_loj_lenght", header=None)


import numpy as np
df = pd.DataFrame({'L' : []})
df['L_counts'] = df2[0].value_counts()
del(df['L'])

df['L'] = df.index
df['product'] = df['L'] * df['L_counts']
num = df['product'].sum()

### Count the total number of PCS bases
df3 = pd.read_csv("/bucket/MillerU/Zifcakova/genome_coordinates2/hg19_gorGor3_new/hg19_gorGor3_all_identical_seqs.bed", sep='\t', header=None)
del(df3[0])
del(df3[1])
del(df3[2])
del(df3[4])
del(df3[5])
del(df3[6])
del(df3[7])
del(df3[8])
del(df3[9])
df3 = df3.astype({3:'int64'})

df4 = pd.DataFrame({'L' : []})
df4['L_counts'] = df3[3].value_counts()
del(df4['L'])
df4['L'] = df4.index
df4['product'] = df4['L'] * df4['L_counts']
deno = df4['product'].sum()

print(num)
print(deno)
print((num/deno)*100,'%')