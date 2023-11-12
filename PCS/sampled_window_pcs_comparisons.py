


import os
import glob
import pandas as pd
#genome_list = ['rheMac10','papAnu4','macFas5','chlSab2','nasLar1','cerAty1','manLeu1','macNem1','colAng1','rhiBie1','rhiRox1']
#names = ['Rhesus','Baboon','Crab-eating Macaque','Green Monkey','Proboscis Monkey','Sooty Mangabey','Drill','Pig-tailed Macaque','Angolan Colobus','Black Snub-nosed Monkey','Golden Snub-nosed Monkey']
#genome_list = ['gorGor6']
#names = ['Gorilla']
genome_list = ['nasLar1','cerAty1','manLeu1','macNem1','colAng1','rhiBie1','rhiRox1']
names = ['Proboscis Monkey','Sooty Mangabey','Drill','Pig-tailed Macaque','Angolan Colobus','Black Snub-nosed Monkey','Golden Snub-nosed Monkey']
#window_size = int(10000000/1000) #(mb -> kb) [10000000, 5000000, 2500000, 1200000, 15000, 50000, 100000]
for genome,name in zip(genome_list,names):
    for window_size in [10000, 5000, 2500, 1200, 500, 100, 50]: 
        for i in [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X','Y']:
            lines=[]
            path = '/bucket/MillerU/Abrar/PCS/hg38-{0}/sampled_window_pcs_comparisons_with_hgmm/{1}kb/{1}kb_chr{2}/'.format(genome,window_size,i)
            os.makedirs(path)
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
            dct1 = {}
            dct2 = {}
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
                    
            import matplotlib.pyplot as plt    
            from matplotlib.colors import hsv_to_rgb
            from cycler import cycler    
            # Resizing the figure
            plt.figure(figsize=[10, 7])
            
            # 1000 distinct colors:
            colors = [hsv_to_rgb([(i * 0.618033988749895) % 1.0, 1, 1])
                      for i in range(1000)]
            plt.rc('axes', prop_cycle=(cycler('color', colors)))
            print(len(lines_frame), len(dct1))
            for k in range(len(lines_frame)):                
                plt.loglog(dct1['s{0}_{1}'.format(i,k)], dct2['t{0}_{1}'.format(i,k)], marker='.', alpha=0.5, linewidth=2, linestyle='None')
                  
            df2 = pd.read_csv("/bucket/MillerU/Abrar/PCS/hg38-{0}/PCS_hg38{0}.csv".format(genome))            
            x2 = [i for i in df2.index]
            y2 = df2['Number of PCS of Length L']
            df3 = pd.read_csv("/bucket/MillerU/Abrar/PCS/hg38-mm39/PCS_hg38mm39.csv")            
            x3 = [i for i in df3.index]
            y3 = df3['Number of PCS of Length L']
            
            plt.loglog(x2, y2, '.g', linewidth=2, label='Human/{0} Vanilla L-mers'.format(name))
            plt.loglog(x3, y3, '.b', linewidth=2, label='Human/Mouse Vanilla L-mers'.format(name))
            
            plt.title('Vanilla PCS Length Distribution Comparison of\n {1}kb Window Sized Samples in hg38/{0} for chr{2}'.format(genome,window_size,i,k+1), fontsize=14)
            
            plt.xlabel('PCS Length', fontsize=20)
            plt.ylabel('Counts', fontsize=20)
            plt.tick_params(axis='x', labelsize=15)
            plt.tick_params(axis='y', labelsize=15)
            
            plt.legend(title="{0} Number of Windows\nEach of Size {1}kb\nFor chr{2}".format(len(dct1),window_size,i), title_fontsize='large', fontsize=13, markerscale=3)
            
            plt.savefig('/bucket/MillerU/Abrar/PCS/hg38-{0}/sampled_window_pcs_comparisons_with_hgmm/{1}kb/{1}kb_chr{2}/hg38{0}_chr{2}_comparisons.png'.format(genome,window_size,i))
            
            
            
            
            
            
            
            
            
            
            
            
            