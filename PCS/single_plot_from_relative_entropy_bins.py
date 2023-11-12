


import os
import glob
import pandas as pd
#genome_list = ['rheMac10','papAnu4','macFas5','chlSab2','nasLar1','cerAty1','manLeu1','macNem1','colAng1','rhiBie1','rhiRox1']
#names = ['Rhesus','Baboon','Crab-eating Macaque','Green Monkey','Proboscis Monkey','Sooty Mangabey','Drill','Pig-tailed Macaque','Angolan Colobus','Black Snub-nosed Monkey','Golden Snub-nosed Monkey']
genome_list = ['gorGor6']
names = ['Gorilla']
dct1 = {}
dct2 = {}
#window_size = int(10000000/1000) #(mb -> kb) [10000000, 5000000, 2500000, 1200000, 15000, 50000, 100000]
for genome,name in zip(genome_list,names):
    for window_size in [10000]:  
        selected_windows = [
        [2,3,10,15,22,23],#1
        [0,1,3,4,7,8,10,11,12,13,23],#2
        [0,6,7,18],#3
        [1,2,3,4,5,6,16,18],#4
        [0,1,2,11,15,16,17],#5
        [0,2,3,6,15,16],#6
        [1,7,10,14],#7
        [2,4,5,8,12,13],#8
        [0,2,3,8,10,11,12],#9 
        [0,1,2,3,4,10,12],#10
        [0,2,4,6,7,8,9,12],#11
        [0,1,3,10,11,12],#12
        [2,8,10],#13
        [2,9],#14
        [3,5,7,8,9],#15
        [2,4,5,7],#16
        [1,3,6,7],#17
        [6,7],#18
        [1,3],#19
        [3,4,5],#20
        [1,2,3],#21
        [3],#22
        [4,7,8,11],#X
        []#Y
        ]
        path = '/bucket/MillerU/Abrar/PCS/hg38-{0}/relative_entropy_pcs_comparisons_with_hgmm/{1}kb/'.format(genome,window_size)
        #os.makedirs(path)
        for i,m in zip([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X','Y'], selected_windows):
            lines=[]

            with open("/bucket/MillerU/Abrar/hg38-{1}_outputs/hg38-{1}_tiled_pcs_{2}kb/chr{0}_{2}kb_min_window_tiled_pcs.txt".format(i, genome,window_size)) as file:
                while True:
                    # Get next line from file
                    line = file.readline().rstrip('\n')
                    # if line is empty
                    # end of file is reached
                    if not line:
                        break
                    elif line==("\n"):
                        continue
                    #print(line)
                    lines.append([line])
            lines_frame = pd.DataFrame(lines, columns=['PCS'])
            
            ########### Single window analysis ############
            
            #backtrack_index1 = [0,1,2,3,4] #one index at a time will generate pcs distribution for a single window
            #colors = ['blue','gold','red','orange','black']

            for k in m:
                print(i,k)
                pcs=[]
                pcs.append(lines_frame.iloc[k].values.tolist())
                flat_pcs = [item for sublist in pcs for item in sublist]
                final_pcs = [int(k) for k in ','.join(flat_pcs).split(',')]        
                serial_count = []
                for j in range(0, int(max(final_pcs)+5)):
                    serial_count.append(final_pcs.count(j))
                
                dct1['s{0}_{1}'.format(i,k)] = [x for x in range(int(max(final_pcs)+5))]
                dct2['t{0}_{1}'.format(i,k)] = serial_count
                   
        import matplotlib.pyplot as plt    
        from matplotlib.colors import hsv_to_rgb
        from cycler import cycler    
        # Resizing the figure
        plt.figure(figsize=[10, 7])
        
        # 1000 distinct colors:
        colors = [hsv_to_rgb([(i * 0.618033988749895) % 1.0, 1, 1])
                  for i in range(1000)]
        plt.rc('axes', prop_cycle=(cycler('color', colors)))
        print((dct1))
        for k,t in zip(dct1,dct2):           
            plt.loglog(dct1[k], dct2[t], marker='.', alpha=0.5, linewidth=2, linestyle='None')     
            #plt.loglog(dct1['s{0}_{1}'.format(i,k)], dct2['t{0}_{1}'.format(i,k)], marker='.', alpha=0.5, linewidth=2, linestyle='None')
            
        df2 = pd.read_csv("/bucket/MillerU/Abrar/PCS/hg38-{0}/PCS_hg38{0}.csv".format(genome))          
        x2 = [i for i in df2.index]
        y2 = df2['Number of PCS of Length L']
        df3 = pd.read_csv("/bucket/MillerU/Abrar/PCS/hg38-mm39/PCS_hg38mm39.csv")            
        x3 = [i for i in df3.index]
        y3 = df3['Number of PCS of Length L']
        
        plt.loglog(x2, y2, '.g', linewidth=2, label='Human/{0} Vanilla L-mers'.format(name))
        plt.loglog(x3, y3, '.b', linewidth=2, label='Human/Mouse Vanilla L-mers'.format(name))
        
        ### chosen idx from chrY references: [1:10MB,2:5000kb,3:2500kb,12:1200kb,14:500kb]
        plt.title('Vanilla PCS Length Distribution Comparison of{1}kb Window Sized Samples\n For Relative Entropy Between (0.7, 0.9] vs chrY 2nd Window in hg38/{0}'.format(genome,window_size), fontsize=14)
        
        plt.xlabel('PCS Length', fontsize=20)
        plt.ylabel('Counts', fontsize=20)
        plt.tick_params(axis='x', labelsize=15)
        plt.tick_params(axis='y', labelsize=15)
        
        plt.legend(title="{0} Number of Windows\nEach of Size {1}kb\nFor Relative Entropy Between (0.7, 0.9]".format(len(dct1),window_size), title_fontsize='large', fontsize=13, markerscale=3)
        
        plt.savefig('/bucket/MillerU/Abrar/PCS/hg38-{0}/relative_entropy_pcs_comparisons_with_hgmm/{1}kb/hg38{0}_RE_3rd_bin_comparisons.png'.format(genome,window_size))
            
            
            
            
            
            
            
            
            
            
            
            
            