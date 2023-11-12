


import os
import glob
import pandas as pd
#genome_list = ['rheMac10','papAnu4','macFas5','chlSab2','nasLar1','cerAty1','manLeu1','macNem1','colAng1','rhiBie1','rhiRox1']
#names = ['Rhesus','Baboon','Crab-eating Macaque','Green Monkey','Proboscis Monkey','Sooty Mangabey','Drill','Pig-tailed Macaque','Angolan Colobus','Black Snub-nosed Monkey','Golden Snub-nosed Monkey']
genome_list = ['gorGor6']
names = ['Gorilla']
window_size = int(5000000/1000) #(mb -> kb) [10000000, 5000000, 2500000, 1200000, 500000, 100000, 50000, 30000, 15000]
for genome,name in zip(genome_list,names):
        for i in ['Y',1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X']:
            lines=[]
            path = '/bucket/MillerU/Abrar/PCS/hg38-{0}/{1}kb/{1}kb_chr{2}/'.format(genome,window_size,i)
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
         
            #backtrack_index1 = [0,1,2,3,4] 
            # range(len(lines_frame)) for all windows
            #one index at a time will generate pcs distribution for a single window
            for k in range(len(lines_frame)):
                pcs=[]
                pcs.append(lines_frame.iloc[k].values.tolist())
                flat_pcs = [item for sublist in pcs for item in sublist]
                final_pcs = [int(k) for k in ','.join(flat_pcs).split(',')]        
                serial_count = []
                for j in range(0, int(max(final_pcs)+5)):
                    serial_count.append(final_pcs.count(j))
                
                ######### PLOT ###########
                import matplotlib.pyplot as plt
                df2 = pd.read_csv("/bucket/MillerU/Abrar/PCS/hg38-{0}/PCS_hg38{0}.csv".format(genome))
                s = [x for x in range(int(max(final_pcs)+5))]
                t = serial_count
                x2 = [i for i in df2.index]
                y2 = df2['Number of PCS of Length L']
                df3 = pd.read_csv("/bucket/MillerU/Abrar/PCS/hg38-mm39/PCS_hg38mm39.csv")            
                x3 = [i for i in df3.index]
                y3 = df3['Number of PCS of Length L']
                                          
                # Resizing the figure
                plt.figure(figsize=[10, 7])
                plt.loglog(x3, y3, '.b', linewidth=2, label='Human/Mouse Vanilla L-mers'.format(name))    
                plt.loglog(s, t, '.r', linewidth=2, label='Human/{0} {1}kb windowed L-mers for\n chr{2}, window number = {3}'.format(name,window_size,i,k+1))
                plt.loglog(x2, y2, '.g', linewidth=2, label='Human/{0} Vanilla L-mers'.format(name))
                           
                
                plt.title('Vanilla PCS Length Distribution of\n Window Number = {3} for {1}kb Window Size in hg38/{0} chr{2}'.format(genome,window_size,i,k+1), fontsize=14)
            
                plt.xlabel('PCS Length', fontsize=20)
                plt.ylabel('Counts', fontsize=20)
                plt.tick_params(axis='x', labelsize=15)
                plt.tick_params(axis='y', labelsize=15)
            
                plt.legend(fontsize=13, markerscale=3)
                
                plt.savefig('/bucket/MillerU/Abrar/PCS/hg38-{0}/{1}kb/{1}kb_chr{2}/hg38{0}_chr{2}_idx{3}.png'.format(genome,window_size,i,k))
                plt.close('all')
        
        
        
    
    
    
