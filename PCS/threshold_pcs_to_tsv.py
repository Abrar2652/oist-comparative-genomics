

import pandas as pd
import numpy as np
import os

genome_list = ['calJac4','danRer10','equCab3','felCat9','galGal6','ornAna2','oryCun2','oviAri4','petMar3','susScr11','tarSyr2','thaSir1','triMan1','xenTro10']
cutoff_list = [30,8,12,10,8,7,10,11,9,10,14,8,10,10]
for genome,cutoff in zip(genome_list,cutoff_list):
    for i in ['all','chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY']:
         
        # Vanilla threshold pcs
        df = pd.read_csv("/bucket/MillerU/Abrar/hg38-{1}_outputs/get_identical_seq_loc_hg38-{1}/{0}/hg38_{1}_{0}_identical_seqs.tsv".format(i,genome), sep='\t')
        cutoff = cutoff 
        df = df[df['first.width']>=cutoff] #PCS threshold L>=20'
        path = '/bucket/MillerU/Abrar/hg38-{2}_outputs/get_identical_seq_loc_hg38-{2}_vanilla_L{0}/{1}/'.format(cutoff, i, genome)
        os.makedirs(path)
        df.to_csv(r'/bucket/MillerU/Abrar/hg38-{2}_outputs/get_identical_seq_loc_hg38-{2}_vanilla_L{0}/{1}/hg38_{2}_{1}_identical_seqs.tsv'.format(cutoff, i, genome), sep="\t", index=False) 
         
        '''
        # Filtered threshold pcs
        df = pd.read_csv("/bucket/MillerU/Abrar/hg38-micMur3_outputs/filtered_pcs_hg38-micMur3_all/{0}/hg38_micMur3_{0}_identical_seqs.tsv".format(i), sep='\t')
        df = df[df['first.width']>=20] #PCS threshold L>=20'
        path = '/bucket/MillerU/Abrar/hg38-micMur3_outputs/filtered_pcs_hg38-micMur3/{0}'.format(i)
        os.makedirs(path)
        df.to_csv(r'/bucket/MillerU/Abrar/hg38-micMur3_outputs/filtered_pcs_hg38-micMur3/{0}/hg38_micMur3_{0}_identical_seqs.tsv'.format(i), sep="\t", index=False) 
        
        '''
    
    
    