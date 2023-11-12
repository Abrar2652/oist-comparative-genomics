

import pandas as pd

###### Input is the chrY PCS coordinates of that window following the power law ########

'''
df=pd.read_csv("/content/chrY_idx0_10000kb_coordinates_hg38gorGor6.csv")

def par_df(df):
    if (df['start']>=10001) and (df['end']<=2781479):
        return df['end']-df['start']+1
    else:
        return 0

def nonpar_df(df):
    if (df['start']>=10001) and (df['end']<=2781479):
        return 0
    else:
        return df['end']-df['start']+1

df['PAR'] = df.apply(par_df, axis = 1)
df['nonPAR'] = df.apply(nonpar_df, axis = 1)
'''

######## Plots #######3

for idx in range(21,29):  ##### idx are the power law windows: (window-1) 
        import pandas as pd
        df=pd.read_csv("/content/chrY_idx{0}_100kb_coordinates_hg38gorGor6.csv".format(idx))
        def par_df(df):
            if (df['start']>=10001) and (df['end']<=2781479):
                return df['end']-df['start']+1
            else:
                return 0

        def nonpar_df(df):
            if (df['start']>=10001) and (df['end']<=2781479):
                return 0
            else:
                return df['end']-df['start']+1

        df['PAR'] = df.apply(par_df, axis = 1)
        df['nonPAR'] = df.apply(nonpar_df, axis = 1)


        #!pip install lmfit
        import numpy as np
        import matplotlib.pyplot as plt
        import csv
        import pandas as pd
        import scipy
        genome = 'gorGor6'
        name = 'Gorilla'
        window_size = int(100000/1000)

        ##### PAR #####
        PAR=df['PAR'].values.tolist()
        serial_count = []
        for j in range(0, int(max(PAR)+5)):
            serial_count.append(PAR.count(j))
        s1 = [x for x in range(int(max(PAR)+5))]
        t1 = serial_count
        t1 = [i for i in t1]
        ##### nonPAR #####
        nonPAR=df['nonPAR'].values.tolist()
        serial_count = []
        for j in range(0, int(max(nonPAR)+5)):
            serial_count.append(nonPAR.count(j))
        s2 = [x for x in range(int(max(nonPAR)+5))]
        t2 = serial_count
        t2 = [i for i in t2]
        ##### whole window #####
        win=(df['end']-df['start']+1).values.tolist()
        serial_count = []
        for j in range(0, int(max(win)+5)):
            serial_count.append(win.count(j))
        s3 = [x for x in range(int(max(win)+5))]
        t3 = serial_count
        t3 = [i for i in t3]
        ######### PLOT ###########
        import matplotlib.pyplot as plt

        df2 = pd.read_csv("/content/PCS_hg38{0}.csv".format(genome))
        x2 = [i for i in df2.index]
        y2 = df2['Number of PCS of Length L']
        y2 = [i for i in y2]

        df3 = pd.read_csv("/content/PCS_hg38mm39.csv")            
        x3 = [i for i in df3.index]
        y3 = df3['Number of PCS of Length L']
        y3 = [i for i in y3]

        total1 = 0    
        total2 = 0      
        total3 = 0  
        total4 = 0
        total5 = 0

        for ele in range(0, len(t1)):
            total1 = total1 + t1[ele]
        for ele in range(0, len(t2)):
            total4 = total4 + t2[ele]
        for ele in range(0, len(t3)):
            total5 = total5 + t3[ele]
        for ele in range(0, len(y2)):
            total2 = total2 + y2[ele]
        for ele in range(0, len(y3)):
            total3 = total3 + y3[ele]  

        t1 = [i/total1 for i in t1]
        t2 = [i/total4 for i in t2]
        t3 = [i/total5 for i in t3]
        y2 = [i/total2 for i in y2]  
        y3 = [i/total3 for i in y3]     

        # Resizing the figure
        plt.figure(figsize=[10, 7])
        plt.loglog(x3, y3, '.b', linewidth=2, label='Human/Mouse Vanilla L-mers'.format(name))    
        plt.loglog(s3, t3, '.r', linewidth=2, label='Human/{0} {1}kb windowed\n L-mers for chrY, window number = {2}'.format(name,window_size,idx+1))
        plt.loglog(s1, t1, 'gold', marker='.', linewidth=2, linestyle='None', label='Human/{0} {1}kb windowed\n PAR1 L-mers for chrY, window number = {2}'.format(name,window_size,idx+1))
        plt.loglog(s2, t2, 'teal', marker='.', linewidth=2, linestyle='None', label='Human/{0} {1}kb windowed\n non-PAR L-mers for chrY, window number = {2}'.format(name,window_size,idx+1))
        plt.loglog(x2, y2, '.g', linewidth=2, label='Human/{0} Vanilla L-mers'.format(name))
                    

        plt.title('All, Par1, and non-PAR Vanilla PCS Length Normalized Distribution\nComparison of Window = {2} for {1}kb Window Size in hg38/{0} chrY'.format(genome,window_size,idx+1), fontsize=14)

        plt.xlabel('PCS Length', fontsize=20)
        plt.ylabel('Counts', fontsize=20)
        plt.tick_params(axis='x', labelsize=15)
        plt.tick_params(axis='y', labelsize=15)

        plt.legend(fontsize=11, markerscale=3)


        plt.savefig('hggg_chrY_{1}kb_allParnonPar_idx{0}.png'.format(idx,window_size))
        plt.show()
        
        

######## download files from colab to local ##########

!zip -r /content/file.zip  /content/
from google.colab import files
files.download("/content/file.zip")
























