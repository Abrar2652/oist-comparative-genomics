#!/usr/bin/env python
# coding: utf-8

# ## Md. Abrar Jahin
# ### Research Intern, Miller Unit, OIST
# 
# This script contains alignments of the following assemblies:
# 
#   - target/reference: Human (hg19, Feb. 2009 (GRCh37/hg19), GRCh37 Genome Reference Consortium Human Reference 37 (GCA_000001405.1))
# 
#   - query: Mouse (mm10, Dec. 2011 (GRCm38/mm10), Genome Reference Consortium Mouse Build 38 (GCA_000001635.2))
# 

# In[1]:


#get_ipython().system('ls')


# In[3]:

import numpy as np
import bx.align.axt
import os
genome_list = ['calJac4','danRer10','equCab3','felCat9','galGal6','ornAna2','oryCun2','oviAri4','petMar3','susScr11','tarSyr2','thaSir1','triMan1','xenTro10']
for genome in genome_list:
    count = 0
    lines=[]
    for i in [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X','Y']:
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
                if (count%3-1)!=0:
                    #print("Line{}: {}".format(count, line.strip()))
                    lines.append(line)
    

    print(len(lines))
    
    import pandas as pd
    lines_frame = pd.DataFrame(lines, columns=['Genome Sequence'])
    lines_frame
    

    import re
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
            if ((letter1.upper() == letter2.upper()) and letter1.upper() in ['A','T','C','G'] and letter2.upper() in ['A','T','C','G']):
                result1+=letter1
                result2+=letter2
            else:
                result1+='\n'
                result2+='\n'
    
        word=re.sub(r'\n+', '\n', result1).strip()
        return word
    
    
    # seq_count stores the PCS lengths
    seq_count = []
    for i in range(int(len(lines_frame))-1): #Iterates throught the whole genome
        if i%2==0: # only even rows are input, but inside the function both even and odd lines are compared
            for x in (compare_strings(lines_frame, "Genome Sequence", i).split()):
                seq_count.append(len(x))
    

    # serial_count counts the total number of a individual PCS lengths
    serial_count=[]
    for i in range(0, int(max(seq_count)+5)):
        serial_count.append(seq_count.count(i))
    
    
    lmer = [x for x in range(int(max(seq_count)+5))]
    # dictionary of lists  
    dict = {'PCS length (L)': lmer, 'Number of PCS of Length L': serial_count}  
    df = pd.DataFrame(dict) 
    # saving the dataframe 
    df.to_csv(r'/bucket/MillerU/Abrar/PCS/hg38-{0}/PCS_hg38{0}_small_to_upper.csv'.format(genome), index=False) 



