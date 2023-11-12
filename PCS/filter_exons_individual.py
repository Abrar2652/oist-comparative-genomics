
import os
import pandas as pd
import numpy as np

def range_subset(range1, range2):
    """Whether range1 is a subset of range2."""
    if not range1:
        return True  # empty range is subset of anything
    if not range2:
        return False  # non-empty range can't be subset of empty range
    if len(range1) > 1 and range1.step % range2.step:
        return False  # must have a single value or integer multiple step
    return range1.start in range2 and range1[-1] in range2
    
# import HCNE bed files    
percent = str(100) # CNE identification threshold
bed = pd.read_csv(("/bucket/MillerU/Abrar/HCNE_hg19_mm10/" + percent + "/HCNE_hg19_mm10_" + percent + "pc_30col.bed"), sep='\t', skiprows=1, header=None) # Mouse and others
#bed = pd.read_csv(("/bucket/MillerU/Abrar/HCNE_hg19_mm10/" + percent + "_500col/hg19_mm10_" + percent + "pc_500col.bed"), sep='\t', header=None) # Gorilla
bed = bed.rename(columns={0: 'first.seqnames', 1: 'first.start', 2:'first.end'})

for i in [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X','Y']: #All_files_merged for whole genome PCS calculation
    pcs = pd.read_csv('/bucket/MillerU/Abrar/hg19-mm10_outputs/get_identical_seq_loc_hg19-mm10_all/chr{0}/hg19_mm10_chr{0}_identical_seqs.tsv'.format(i), sep='\t')
    path = '/bucket/MillerU/Abrar/hg19-mm10_outputs/filtered_pcs_hg19-mm10_all/chr{0}'.format(i)
    os.mkdir(path) 
    chrom = 'chr{}'.format(i)
    print(chrom)
    temp_bed = bed[bed['first.seqnames']==chrom]
    
    # for array broadcasting
    m = temp_bed['first.start'].to_numpy()[:, None]
    n = temp_bed['first.end'].to_numpy()[:, None]
    
    # A chunk_size that is too small or too big will lower performance.
    # Experiment to find a sweet spot
    chunk_size = 500_00
    offset = 0
    mask = []
    mask1 = []
    mask2 = []
    mask3 = []
    
    while offset < len(pcs):
        x = pcs['first.start'].to_numpy()[offset:offset+chunk_size] #main first_start
        y = pcs['first.end'].to_numpy()[offset:offset+chunk_size] #main first_end

        # Converting pcs lengths to their 50% because we'll filter the exons having atleast 50% overlaps 
        pcs_length = pcs['first.width'].to_numpy()[offset:offset+chunk_size]
        pcs_percent = [(.5 * i) for i in pcs_length]

        mask.append(((m <= x) & (n >= y) ).any(axis=0))
        mask1.append((((m <= x) & (x <= n)).any(axis=0)) & ((n-x+1) > pcs_percent).any(axis=0))
        mask2.append((((m >= x) & (y >= n)).any(axis=0)) & ((n-m+1) > pcs_percent).any(axis=0))
        mask3.append((((m >= x) & (y <= n) & (y >= m)).any(axis=0)) & ((y-m+1) > pcs_percent).any(axis=0))
        
        offset += chunk_size
    

    mask = np.hstack(mask)
    mask1 = np.hstack(mask1)
    mask2 = np.hstack(mask2)
    mask3 = np.hstack(mask3)
    mask_all = (mask | mask1 | mask2 | mask3)
    
    print(pcs[mask_all])
    pcs[~mask_all].to_csv(r'/bucket/MillerU/Abrar/hg19-mm10_outputs/filtered_pcs_hg19-mm10_all/chr{0}/hg19_mm10_chr{0}_identical_seqs.tsv'.format(i), sep="\t", index=False)


   
   
   
   
   
   
   
   
   
    
