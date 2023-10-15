
import os
import pandas as pd
percent = '100' # 100, 98, 96.6, 96, 90
window_size_list = ['15','30'] #in kb  # 15, 120, 30, 60
species1 = 'hg38'
#species2 = 'micMur3' 
#column = '500' #For Gorilla
genome_list = ['calJac4','danRer10','equCab3','felCat9','galGal6','ornAna2','oryCun2','oviAri4','petMar3','susScr11','tarSyr2','thaSir1','triMan1','xenTro10']

for species2 in genome_list:
  for window_size in window_size_list:
    #####################   KURTOSIS (QUANTILE) #######################
    
    stat = 'quantile_kurtosis'
    #parent_dir = '/bucket/MillerU/Abrar/' + species1 + '-' + species2 + '_outputs/' + window_size + 'kb/binned_' + stat +'_genome_ave_' + percent + 'pc_' + column + 'col/' #for Gorilla
    parent_dir = '/bucket/MillerU/Abrar/' + species1 + '-' + species2 + '_outputs/' + window_size + 'kb/binned_' + stat +'_genome_ave_' + percent + 'pc/' #for others
    appended_data = []
    chrom = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X','Y']
    for i in chrom:
        data = pd.read_csv(parent_dir + species1 + '_' + species2 +'_chr{}_'.format(i)+ window_size+ 'kb_min_window.wig', header=None, index_col=False, skiprows = 1)
        # store DataFrame in list
        appended_data.append(data)
    # see pd.concat documentation for more info
    appended_data = pd.concat(appended_data)
    # write DataFrame to an excel sheet 
    appended_data = appended_data.reset_index()
    os.makedirs('/bucket/MillerU/Abrar/Moment_val_bins/' + species1 + species2 + '/All_bin_' + window_size + 'kb/')
    #appended_data.to_csv('/bucket/MillerU/Abrar/Moment_val_bins/' + species1 + species2 + '/All_bin_' + window_size + 'kb/' + stat + '_' + percent + 'pc_' + column + 'col' + '.csv', index=False) #for Gorilla
    appended_data.to_csv('/bucket/MillerU/Abrar/Moment_val_bins/' + species1 + species2 + '/All_bin_' + window_size + 'kb/' + stat + '_' + percent + 'pc.csv', index=False) #for others
    
    
    #####################   SKEWNESS (QUANTILE) #######################
    
    stat = 'quantile_skewness'
    #parent_dir = '/bucket/MillerU/Abrar/' + species1 + '-' + species2 + '_outputs/' + window_size + 'kb/binned_' + stat +'_genome_ave_' + percent + 'pc_' + column + 'col/' #for Gorilla
    parent_dir = '/bucket/MillerU/Abrar/' + species1 + '-' + species2 + '_outputs/' + window_size + 'kb/binned_' + stat +'_genome_ave_' + percent + 'pc/' #for others
    appended_data = []
    chrom = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X','Y']
    for i in chrom:
        data = pd.read_csv(parent_dir + species1 + '_' + species2 +'_chr{}_'.format(i)+ window_size+ 'kb_min_window.wig', header=None, index_col=False, skiprows = 1)
        # store DataFrame in list
        appended_data.append(data)
    # see pd.concat documentation for more info
    appended_data = pd.concat(appended_data)
    # write DataFrame to an excel sheet 
    appended_data = appended_data.reset_index()
    #appended_data.to_csv('/bucket/MillerU/Abrar/Moment_val_bins/' + species1 + species2 + '/All_bin_' + window_size + 'kb/' + stat + '_' + percent + 'pc_' + column + 'col' + '.csv', index=False) #for Gorilla
    appended_data.to_csv('/bucket/MillerU/Abrar/Moment_val_bins/' + species1 + species2 + '/All_bin_' + window_size + 'kb/' + stat + '_' + percent + 'pc.csv', index=False) #for others
    

'''
#####################   MEAN #######################

stat = 'mean'
#parent_dir = '/bucket/MillerU/Abrar/' + species1 + '-' + species2 + '_outputs/' + window_size + 'kb/binned_' + stat +'_genome_ave_' + percent + 'pc_' + column + 'col/' #for Gorilla
parent_dir = '/bucket/MillerU/Abrar/' + species1 + '-' + species2 + '_outputs/' + window_size + 'kb/binned_' + stat +'_genome_ave_' + percent + 'pc/'
appended_data = []
chrom = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X','Y']
for i in chrom:
    data = pd.read_csv(parent_dir + species1 + '_' + species2 +'_chr{}_'.format(i)+ window_size+ 'kb_min_window.wig', header=None, index_col=False, skiprows = 1)
    # store DataFrame in list
    appended_data.append(data)
# see pd.concat documentation for more info
appended_data = pd.concat(appended_data)
# write DataFrame to an excel sheet 
appended_data = appended_data.reset_index()
#appended_data.to_csv('/bucket/MillerU/Abrar/Moment_val_bins/' + species1 + species2 + '/All_bin_' + window_size + 'kb/' + stat + '_' + percent + 'pc_' + column + 'col' + '.csv', index=False) #for Gorilla
appended_data.to_csv('/bucket/MillerU/Abrar/Moment_val_bins/' + species1 + species2 + '/All_bin_' + window_size + 'kb/' + stat + '_' + percent + 'pc.csv', index=False) #for others



#####################   VARIANCE #######################

stat = 'variance'
#parent_dir = '/bucket/MillerU/Abrar/' + species1 + '-' + species2 + '_outputs/' + window_size + 'kb/binned_' + stat +'_genome_ave_' + percent + 'pc_' + column + 'col/' #for Gorilla
parent_dir = '/bucket/MillerU/Abrar/' + species1 + '-' + species2 + '_outputs/' + window_size + 'kb/binned_' + stat +'_genome_ave_' + percent + 'pc/'
appended_data = []
chrom = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X','Y']
for i in chrom:
    data = pd.read_csv(parent_dir + species1 + '_' + species2 +'_chr{}_'.format(i)+ window_size+ 'kb_min_window.wig', header=None, index_col=False, skiprows = 1)
    # store DataFrame in list
    appended_data.append(data)
# see pd.concat documentation for more info
appended_data = pd.concat(appended_data)
# write DataFrame to an excel sheet 
appended_data = appended_data.reset_index()
#appended_data.to_csv('/bucket/MillerU/Abrar/Moment_val_bins/' + species1 + species2 + '/All_bin_' + window_size + 'kb/' + stat + '_' + percent + 'pc_' + column + 'col' + '.csv', index=False) #for Gorilla
appended_data.to_csv('/bucket/MillerU/Abrar/Moment_val_bins/' + species1 + species2 + '/All_bin_' + window_size + 'kb/' + stat + '_' + percent + 'pc.csv', index=False) #for others



#####################   SKEWNESS #######################

stat = 'skewness'
#parent_dir = '/bucket/MillerU/Abrar/' + species1 + '-' + species2 + '_outputs/' + window_size + 'kb/binned_' + stat +'_genome_ave_' + percent + 'pc_' + column + 'col/' #for Gorilla
parent_dir = '/bucket/MillerU/Abrar/' + species1 + '-' + species2 + '_outputs/' + window_size + 'kb/binned_' + stat +'_genome_ave_' + percent + 'pc/'
appended_data = []
chrom = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X','Y']
for i in chrom:
    data = pd.read_csv(parent_dir + species1 + '_' + species2 +'_chr{}_'.format(i)+ window_size+ 'kb_min_window.wig', header=None, index_col=False, skiprows = 1)
    # store DataFrame in list
    appended_data.append(data)
# see pd.concat documentation for more info
appended_data = pd.concat(appended_data)
# write DataFrame to an excel sheet 
appended_data = appended_data.reset_index()
#appended_data.to_csv('/bucket/MillerU/Abrar/Moment_val_bins/' + species1 + species2 + '/All_bin_' + window_size + 'kb/' + stat + '_' + percent + 'pc_' + column + 'col' + '.csv', index=False) #for Gorilla
appended_data.to_csv('/bucket/MillerU/Abrar/Moment_val_bins/' + species1 + species2 + '/All_bin_' + window_size + 'kb/' + stat + '_' + percent + 'pc.csv', index=False) #for others



#####################  EXCESS KURTOSIS (FISHER) #######################

stat = 'excess_kurtosis'
#parent_dir = '/bucket/MillerU/Abrar/' + species1 + '-' + species2 + '_outputs/' + window_size + 'kb/binned_' + stat +'_genome_ave_' + percent + 'pc_' + column + 'col/' #for Gorilla
parent_dir = '/bucket/MillerU/Abrar/' + species1 + '-' + species2 + '_outputs/' + window_size + 'kb/binned_' + stat +'_genome_ave_' + percent + 'pc/'
appended_data = []
chrom = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X','Y']
for i in chrom:
    data = pd.read_csv(parent_dir + species1 + '_' + species2 +'_chr{}_'.format(i)+ window_size+ 'kb_min_window.wig', header=None, index_col=False, skiprows = 1)
    # store DataFrame in list
    appended_data.append(data)
# see pd.concat documentation for more info
appended_data = pd.concat(appended_data)
# write DataFrame to an excel sheet 
appended_data = appended_data.reset_index()
#appended_data.to_csv('/bucket/MillerU/Abrar/Moment_val_bins/' + species1 + species2 + '/All_bin_' + window_size + 'kb/' + stat + '_' + percent + 'pc_' + column + 'col' + '.csv', index=False) #for Gorilla
appended_data.to_csv('/bucket/MillerU/Abrar/Moment_val_bins/' + species1 + species2 + '/All_bin_' + window_size + 'kb/' + stat + '_' + percent + 'pc.csv', index=False) #for others



#####################   6TH MOMENT #######################


stat = 'sixth_moment'
#parent_dir = '/bucket/MillerU/Abrar/' + species1 + '-' + species2 + '_outputs/' + window_size + 'kb/binned_' + stat +'_genome_ave_' + percent + 'pc_' + column + 'col/' #for Gorilla
parent_dir = '/bucket/MillerU/Abrar/' + species1 + '-' + species2 + '_outputs/' + window_size + 'kb/binned_' + stat +'_genome_ave_' + percent + 'pc/'
appended_data = []
chrom = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X','Y']
for i in chrom:
    data = pd.read_csv(parent_dir + species1 + '_' + species2 +'_chr{}_'.format(i)+ window_size+ 'kb_min_window.wig', header=None, index_col=False, skiprows = 1)
    # store DataFrame in list
    appended_data.append(data)
# see pd.concat documentation for more info
appended_data = pd.concat(appended_data)
# write DataFrame to an excel sheet 
appended_data = appended_data.reset_index()
#appended_data.to_csv('/bucket/MillerU/Abrar/Moment_val_bins/' + species1 + species2 + '/All_bin_' + window_size + 'kb/' + stat + '_' + percent + 'pc_' + column + 'col' + '.csv', index=False) #for Gorilla
appended_data.to_csv('/bucket/MillerU/Abrar/Moment_val_bins/' + species1 + species2 + '/All_bin_' + window_size + 'kb/' + stat + '_' + percent + 'pc.csv', index=False) #for others


'''










