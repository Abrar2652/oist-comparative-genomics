#!/usr/bin/env python
# coding: utf-8

# ## Md. Abrar Jahin
# ### Research Intern, Miller Unit, OIST


#get_ipython().system('ls')


# In[2]:


import bx.align.axt
import os
count = 0
lines=[]
for i in [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X','Y']:
    with open('/bucket/MillerU/Abrar/axtNet/hg38nomLeu3/chr{0}.hg38.nomLeu3.net.axt'.format(i)) as file:
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


# In[3]:


print(len(lines))


# In[4]:


import pandas as pd
lines_frame = pd.DataFrame(lines, columns=['Genome Sequence'])
lines_frame


# In[5]:


len(lines_frame)/2


# In[6]:


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

# seq_count stores the PCS lengths
seq_count = []
for i in range(int(len(lines_frame))-1): #Iterates throught the whole genome
    if i%2==0: # only even rows are input, but inside the function both even and odd lines are compared
        for x in (compare_strings(lines_frame, "Genome Sequence", i).split()):
            seq_count.append(len(x))


# In[7]:


print(len(seq_count))


# In[8]:
# serial_count counts the total number of a individual PCS lengths
# For mouse, range=1000. dog:1050, gorilla:2500 (their maximum PCS are close to these values)
serial_count=[]
for i in range(0, int(max(seq_count)+5)):
    #print("count of {}: {}".format(i,seq_count.count(i)))
    serial_count.append(seq_count.count(i))


lmer = [x for x in range(int(max(seq_count)+5))]
# dictionary of lists  
dict = {'PCS length (L)': lmer, 'Number of PCS of Length L': serial_count}  
       
df = pd.DataFrame(dict) 
    
# saving the dataframe 
df.to_csv(r'/bucket/MillerU/Abrar/PCS/hg38-nomLeu3/PCS_hg38nomLeu3.csv', index = False) 








