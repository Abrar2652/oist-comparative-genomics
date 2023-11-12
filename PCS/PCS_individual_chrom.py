#!/usr/bin/env python
# coding: utf-8

# ## Md. Abrar Jahin
# ### Research Intern, Miller Unit, OIST

#get_ipython().system('ls')


import bx.align.axt
import os
genome_list = ['panPan3','panTro6']
chrom_list = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X','Y'] # Iterate through each chromosome axt files individually
for genome in genome_list:
    for chrom in chrom_list:
        count = 0
        lines=[]
        header=[]    
        with open('/bucket/MillerU/Abrar/axtNet/hg38{1}/chr{0}.hg38.{1}.net.axt'.format(chrom, genome)) as file:
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
        header_frame = header_frame['header'].str.split(' ', 9, expand=True).rename(columns={0:'index', 1:'first.seqnames', 
                                                                             2:"first.start.old", 3:"first.end.old",
                                                                             4:"second.seqnames", 5:"second.start.old",
                                                                             6:"second.end.old"
                                                                             })
        
        header_frame.drop('index', axis=1, inplace=True)
        header_frame.drop(7, axis=1, inplace=True)
        header_frame.drop(8, axis=1, inplace=True)
        
        #Inserting a new column at desired location
        #df.insert(loc=idx, column='A', value=new_col)
        header_frame.insert(loc=3, column="first.strand", value="*")
        header_frame.insert(loc=7, column="second.strand", value="*")
        header_frame = header_frame.astype({'first.seqnames':'str', 'first.start.old':'int64', 
                           'first.end.old':'int64', 'first.strand':'str',
                           'second.seqnames':'str', 'second.start.old':'int64', 
                           'second.end.old':'int64', 'second.strand':'str'})
        
        
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
                if (letter1.islower() and letter2.islower() and letter1==letter2):
                    result1+='\n'
                    result2+='\n'
                if ((letter1 not in ['A','T','C','G']) or (letter2 not in ['A','T','C','G'])) and letter1==letter2:
                    result1+='\n'
                    result2+='\n'
                if ((letter1 == letter2) and letter1.isupper() and letter2.isupper() and (letter1!='-') and (letter2!='-') and 
                    letter1 in ['A','T','C','G'] and letter2 in ['A','T','C','G']):
                    result1+=letter1
                    result2+=letter2
                if ((letter1 != letter2) or (letter1=='-') or (letter2=='-')):
                    result1+='\n'
                    result2+='\n'
        
            word=re.sub(r'\n+', '\n', result1).strip()
            return word
        
        
        def matched_indices(df, column, j):
            str1 = lines_frame["Genome Sequence"][j] #even rows
            str2 = lines_frame["Genome Sequence"][j+1] #odd rows
            
            exclusive1 = [] #new indices skipping indels for 1st sequence
            exclusive2 = [] #new indices skipping indels for 2nd sequence
            count=0
            for idx, val in enumerate(str1):
                if (val not in ['-', '\n', '+', '-', '.']) and (idx!=len(str1)-1) and (str1[idx+1] not in ['-', '\n', '+', '-', '.']): #str is not '-' and not the last idx and next is not '-'
                    exclusive1.append(count)
                    count+=1
                elif (val not in ['-', '\n', '+', '-', '.']) and (idx==len(str1)-1): #str is not '-' and the last idx
                    exclusive1.append(count)
                elif (val in ['-', '\n', '+', '-', '.']) and (idx!=len(str1)-1) and (str1[idx+1] not in ['-', '\n', '+', '-', '.']): #str is '-' and not the last idx and next is not '-'
                    exclusive1.append(count)
                    count+=1        
                else: #str is '-'
                    exclusive1.append(count)
            count=0
            for idx, val in enumerate(str2):
                if (val not in ['-', '\n', '+', '-', '.']) and (idx!=len(str2)-1) and (str2[idx+1] not in ['-', '\n', '+', '-', '.']):
                    exclusive2.append(count)
                    count+=1
                elif (val not in ['-', '\n', '+', '-', '.']) and (idx==len(str2)-1):
                    exclusive2.append(count)
                elif (val in ['-', '\n', '+', '-', '.']) and (idx!=len(str2)-1) and (str2[idx+1] not in ['-', '\n', '+', '-', '.']):
                    exclusive2.append(count)
                    count+=1        
                else:
                    exclusive2.append(count)
        
            matches = []
            for i,(letter1, letter2) in enumerate(zip(str1,str2)):#i=index, letter1=str1, letter2=str2
                if ((letter1 == letter2) and (letter1 in ['A','T','C','G']) and (letter2 in ['A','T','C','G'])):            
                    if not matches or matches[-1][1] != i-1:
                        matches.append([i,i])
                    else:
                        matches[-1][1] += 1
            
            start = [k[0] for k in matches]
            end = [k[1] for k in matches]
            
            first_start = [exclusive1[i] for i in start]
            second_end = [exclusive2[i] for i in end]
            first_end = [exclusive1[i] for i in end]
            second_start = [exclusive2[i] for i in start]
        
            return first_start, first_end, second_start, second_end
        
            
        
        # seq_count stores the PCS lengths
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
                for _ in starmap(list.extend, zip([first_start, first_end, second_start, second_end], matched_indices(lines_frame, "Genome Sequence", i))):
                    pass
                for x in (compare_strings(lines_frame, "Genome Sequence", i).split()):
                    seq_count.append(len(x))
                    seq_count_replica.append(len(x))
                ultimate.append(seq_count_replica)#handle null pcs
        for i in ultimate:
            pcs_per_seq.append(len(i))
        
        
        header_frame['pcs_per_seq'] = pcs_per_seq
        #repeat rows in header PCS_per_seq number of times
        header_frame = header_frame.loc[header_frame.index.repeat(header_frame.pcs_per_seq)].reset_index(drop=True)
        header_frame['first_start'] = first_start
        header_frame['first_end'] = first_end
        header_frame['second_start'] = second_start
        header_frame['second_end'] = second_end
        header_frame['first.start.r'] = header_frame['first.start.old'] + header_frame['first_start'] 
        header_frame['first.end.r'] = header_frame['first.start.old'] + header_frame['first_end'] 
        header_frame['second.start.r'] = header_frame['second.start.old'] + header_frame['second_start'] 
        header_frame['second.end.r'] = header_frame['second.start.old'] + header_frame['second_end']
        del(first_start)
        del(first_end)
        del(second_start)
        del(second_end)
        del(seq_count_replica)
        del(ultimate)
        header_frame.drop('first.start.old', inplace=True, axis=1)
        header_frame.drop('first.end.old', inplace=True, axis=1)
        header_frame.drop('second.start.old', inplace=True, axis=1)
        header_frame.drop('second.end.old', inplace=True, axis=1)
        header_frame.drop('first_start', inplace=True, axis=1)
        header_frame.drop('first_end', inplace=True, axis=1)
        header_frame.drop('second_start', inplace=True, axis=1)
        header_frame.drop('second_end', inplace=True, axis=1)
        header_frame.drop('pcs_per_seq', inplace=True, axis=1)
        header_frame.insert(loc=1, column = "first.start", value = header_frame["first.start.r"])
        header_frame.insert(loc=2, column = "first.end", value = header_frame["first.end.r"])
        header_frame.insert(loc=3, column = "first.width", value = seq_count)
        header_frame.insert(loc=6, column = "second.start", value = header_frame["second.start.r"])
        header_frame.insert(loc=7, column = "second.end", value = header_frame["second.end.r"])
        header_frame.insert(loc=8, column = "second.width", value = seq_count)
        header_frame.drop('first.start.r', inplace=True, axis=1)
        header_frame.drop('first.end.r', inplace=True, axis=1)
        header_frame.drop('second.start.r', inplace=True, axis=1)
        header_frame.drop('second.end.r', inplace=True, axis=1)
        
        path = '/bucket/MillerU/Abrar/hg38-{1}_outputs/get_identical_seq_loc_hg38-{1}/chr{0}'.format(chrom, genome)
        os.makedirs(path) 
        header_frame.to_csv(r'/bucket/MillerU/Abrar/hg38-{1}_outputs/get_identical_seq_loc_hg38-{1}/chr{0}/hg38_{1}_chr{0}_identical_seqs.tsv'.format(chrom, genome), sep="\t", index=False)
    
        print(len(seq_count))
        
    






