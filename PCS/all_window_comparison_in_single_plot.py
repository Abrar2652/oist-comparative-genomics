                                                                  


import os
import glob
import pandas as pd

genome_list = ['mm39']
names = ['Mouse']
#window_size = int(500000/1000) #(mb -> kb) [10000000, 5000000, 2500000, 1200000, 500000, 15000, 50000, 100000]
for genome,name in zip(genome_list,names):
    #for window_size in [10000, 5000, 2500, 1200, 15, 50, 100, 30]: 
    for window_size in [500]: 
        dct1 = {}
        dct2 = {}
        k_s = []
        for i in [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X','Y']:
            lines=[]
            path = '/bucket/MillerU/Abrar/PCS/hg38-{0}/all_window_comparisons/'.format(genome,window_size,i)
            #os.makedirs(path)
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
            print(lines_frame.shape)
            ########### Single window analysis ############
    
            #backtrack_index1 = [0,1,2,3,4] #one index at a time will generate pcs distribution for a single window
            k_s.append(len(lines_frame))
            print(k_s)
            #for k in backtrack_index1:
            for k in range(len(lines_frame)):
                pcs=[]
                pcs.append(lines_frame.iloc[k].values.tolist())
                flat_pcs = [item for sublist in pcs for item in sublist]
                final_pcs = [int(k) for k in ','.join(flat_pcs).split(',')]        
                serial_count = []
                for j in range(0, int(max(final_pcs)+5)):
                    serial_count.append(final_pcs.count(j))
                dct1['s{0}_{1}'.format(i,k)] = [x for x in range(int(max(final_pcs)+5))]
                dct2['t{0}_{1}'.format(i,k)] = serial_count
                  
        from matplotlib import pyplot as plt
        from matplotlib.colors import hsv_to_rgb
        from cycler import cycler      
        # Resizing the figure 
        plt.figure(figsize=[10, 7])
        # 1000 distinct colors:
        colors = [hsv_to_rgb([(i * 0.618033988749895) % 1.0, 1, 1])
                  for i in range(1000)]
        plt.rc('axes', prop_cycle=(cycler('color', colors)))
        print(len(dct1),len(k_s))
        for i,j in zip([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X','Y'], k_s):
            for k in range(j):
              plt.loglog(dct1['s{0}_{1}'.format(i,k)], dct2['t{0}_{1}'.format(i,k)], marker='.', alpha=0.5, linewidth=2, linestyle='None')
              #print(dct1['s{0}_{1}'.format(i,k)])
        
        df2 = pd.read_csv("/bucket/MillerU/Abrar/PCS/hg38-{0}/PCS_hg38{0}.csv".format(genome))
        x2 = [i for i in df2.index]
        y2 = df2['Number of PCS of Length L']
        plt.loglog(x2, y2, 'black', marker='.', linewidth=2, linestyle='None', label='Human/{0} Vanilla L-mers'.format(name))
        
        plt.title('Vanilla PCS Length Distribution Comparison of All\n {1}kb Window Sized Samples in hg38/{0} for All Chromosomes'.format(genome,window_size,i,k+1), fontsize=14)
        
        plt.xlabel('PCS Length', fontsize=20)
        plt.ylabel('Counts', fontsize=20)
        plt.tick_params(axis='x', labelsize=15)
        plt.tick_params(axis='y', labelsize=15)
        
        plt.legend(title="{0} Number of Windows\nEach of Size {1}kb".format(len(dct1),window_size), title_fontsize='large', fontsize=13, markerscale=3)
        plt.savefig('/bucket/MillerU/Abrar/PCS/hg38-{0}/all_window_comparisons/hg38{0}_{1}kb_all_window_comparisons.png'.format(genome,window_size,i))
            

        
        
    
    
