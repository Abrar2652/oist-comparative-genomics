#!/usr/bin/env python
# coding: utf-8

# ## Md. Abrar Jahin
# ### Research Intern, Miller Unit, OIST

#get_ipython().system('ls')

import bx.align.axt
import os

genome_list = ['mm39','gorGor6']
for genome in genome_list:
    count = 0
    lines=[]
    header=[]
    for i in [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X','Y']: #All_files_merged for whole genome PCS calculation
        with open('/bucket/MillerU/Abrar/axtNet/hg38{1}/chr{0}.hg38.{1}.net.axt'.format(i, genome)) as file:
            while True: 
                # Get next line from file
                line = bx.align.axt.readline(file, skip_blank=False)
    
                # if line is empty
                # end of file is reached
                if not line:
                    break
                elif line==("\n"):
                    continue
                count += 1
                if (count%3-1)==0:
                    header.append(line)
                if (count%3-1)!=0:
                    #print("Line{}: {}".format(count, line.strip()))
                    lines.append(line)
            
    import pandas as pd
    lines_frame = pd.DataFrame(lines, columns=['Genome Sequence'])
    header_frame = pd.DataFrame(header, columns=['header'])
    header_frame2 = header_frame['header'].str.split('\n', 2, expand=True)
    
    header_frame = header_frame['header'].str.split(' ', 9, expand=True).rename(columns={0:'index', 1:'first.seqnames', 
                                                                         4:"second.seqnames"
                                                                         })
    
    header_frame.drop('index', axis=1, inplace=True)
    header_frame.drop(7, axis=1, inplace=True)
    header_frame.drop(8, axis=1, inplace=True)
    
    #Inserting a new column at desired location
    #df.insert(loc=idx, column='A', value=new_col)
    header_frame = header_frame.astype({'first.seqnames':'str', 'second.seqnames':'str'})

    header_frame['summary_line'] = header_frame2[0] # including summary line of the alignment block for identifying the PCS under the same sequence blocks
    header_frame.drop(2, axis=1, inplace=True)
    header_frame.drop(3, axis=1, inplace=True)
    header_frame.drop(5, axis=1, inplace=True)
    header_frame.drop(6, axis=1, inplace=True)
    print(header_frame.head())  
        
    import re
    from itertools import starmap
    def compare_strings(df, column, j):
        # First, I removed the split... it is already an array
        str1 = df["Genome Sequence"][j] #even rows
        str2 = df["Genome Sequence"][j+1] #odd rows
    
        #then creating a new variable to store the result after  
        #comparing the strings. You note that I added result2 because 
        #if string 2 is longer than string 1 then you have extra characters 
        #in result 2, if string 1 is  longer then the result you want to take 
        #a look at is result 2
    
        result1 = ''
        result2 = ''
    
        #handle the case where one string is longer than the other
        maxlen=len(str2) if len(str1)<len(str2) else len(str1)
        
        #loop through the characters
        for i in range(maxlen):
          #use a slice rather than index in case one string longer than other
            letter1=str1[i:i+1]
            letter2=str2[i:i+1]
            #create string with differences
            ## Contiguous mis-matched sequences (CMMS), runs of mismatches that are terminated not by mismatches, but by a match.
            if (letter1.isupper() and letter2.isupper() and letter1==letter2):
                result1+='\n'
                result2+='\n'
            if ((letter1 not in ['a','t','c','g']) or (letter2 not in ['a','t','c','g'])) and letter1!=letter2: # smallercase terminate
                result1+='\n'
                result2+='\n'
            if ((letter1!=letter2)  and (letter1!='-') and (letter2!='-') and 
                (letter1 in ['a','t','c','g']) and (letter2 in ['a','t','c','g'])):
                result1+=letter1
                result2+=letter2
            if ((letter1 == letter2) or (letter1=='-') or (letter2=='-')): # match and indel terminate
                result1+='\n'
                result2+='\n'

        word=re.sub(r'\n+', '\n', result1).strip()
        return word
         
    
    # seq_count stores the PCS lengths
    seq=[]    
    seq_count = []
    pcs_per_seq = []
    ultimate = []
    first_start = []
    first_end = []
    second_start = []
    second_end = []
    for i in range(int(len(lines_frame))-1):  #Iterates through the whole genome
        seq_count_replica = []
        if i%2==0: # only even rows are input, but inside the function both even and odd lines are compared
            #store start and end indices: relative ranges in each sequence
            for x in (compare_strings(lines_frame, "Genome Sequence", i).split()):
                seq.append(x)
                seq_count.append(len(x))
                seq_count_replica.append(len(x))
            ultimate.append(seq_count_replica)#handle null pcs
    for i in ultimate:
        pcs_per_seq.append(len(i))
    
    
    header_frame['pcs_per_seq'] = pcs_per_seq
    #repeat rows in header PCS_per_seq number of times
    header_frame = header_frame.loc[header_frame.index.repeat(header_frame.pcs_per_seq)].reset_index(drop=True)

    del(seq_count_replica)
    del(ultimate)
    
    header_frame.drop('pcs_per_seq', inplace=True, axis=1)
    
    header_frame.insert(loc=2, column = "first.sequence", value = seq)
    header_frame.insert(loc=4, column = "second.sequence", value = seq)
    header_frame.insert(loc=5, column = "pcs.width", value = seq_count)
        
    chrom = 'all'

    #header_frame.to_csv(r'/bucket/.mabuya/MillerU/Abrar/PCS/AR/mismatch_terminated/AR_exact_match/hg38_{1}_{0}_AR_mismatch_terminated_exact_match_seqs.tsv'.format(chrom, genome), sep="\t", index=False)
    header_frame.to_csv(r'hg38_{1}_{0}_AR_CMMS_both_match_indel_terminated_seqs.tsv'.format(chrom, genome), sep="\t", index=False)
    
    print(len(seq_count))
    
    # RUN THE FOLLOWING CODE ONLY FOR THE ALL_MERGED FILE
    
    serial_count=[]
    for i in range(0, int(max(seq_count)+5)):
        serial_count.append(seq_count.count(i))
       
    lmer = [x for x in range(int(max(seq_count)+5))]
    # dictionary of lists  
    dict = {'PCS length (L)': lmer, 'Number of PCS of Length L': serial_count}  
           
    df = pd.DataFrame(dict) 
  
    # saving the dataframe 

    #df.to_csv(r'/bucket/.mabuya/MillerU/Abrar/PCS/AR/mismatch_terminated/AR_exact_match/AR_mismatch_terminated_exact_match_PCS_hg38{0}.csv'.format(genome), index = False) 
    df.to_csv(r'AR_CMMS_both_match_indel_terminated_PCS_hg38{0}.csv'.format(genome), index = False)
    
    
    
    
    
    
    




