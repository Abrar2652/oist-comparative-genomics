

import bx.align.axt
import os
count = 0
file = open("/bucket/MillerU/Abrar/axtNet/hg38galVar1/hg38.galVar1.net.axt")
lines=[]
header=[]
while True: 
    # Get next line from file
    line = bx.align.axt.readline(file, skip_blank=False)

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
    #if(count==50):
    #    break
    
    
import pandas as pd
header_frame1 = pd.DataFrame(header, columns=['header info'])
header_frame2=header_frame1.copy()
print(header_frame1.head())

lines_frame = pd.DataFrame(lines, columns=['Genome Sequence'])
print(lines_frame.head())

header_frame2 = header_frame2['header info'].str.split(' ', 9, expand=True).rename(columns={0:'index', 1:'first.seqnames', 
                                                                     2:"first.start.old", 3:"first.end.old",
                                                                     4:"second.seqnames", 5:"second.start.old",
                                                                     6:"second.end.old"
                                                                     })

header_frame2.drop('index', axis=1, inplace=True)
header_frame2.drop(7, axis=1, inplace=True)
header_frame2.drop(8, axis=1, inplace=True)

#Inserting a new column at desired location
#df.insert(loc=idx, column='A', value=new_col)
header_frame2.insert(loc=3, column="first.strand", value="*")
header_frame2.insert(loc=7, column="second.strand", value="*")
header_frame2 = header_frame2.astype({'first.seqnames':'str', 'first.start.old':'int64', 
                   'first.end.old':'int64', 'first.strand':'str',
                   'second.seqnames':'str', 'second.start.old':'int64', 
                   'second.end.old':'int64', 'second.strand':'str'})
                   

header_idx=[i for i in lines_frame.index if i%2==0]

header_frame1['index']=pd.DataFrame(header_idx)
header_frame1=header_frame1.set_index(header_frame1['index'])
header_frame2['index']=pd.DataFrame(header_idx)
header_frame2=header_frame2.set_index(header_frame2['index'])

import csv
for i in [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X','Y']:
    k='chr{}'.format(i)
    #idx->0,1,2,3,... for header
    #even->0,2,4,8,...for header
    final=[]
    for idx, (even, j) in enumerate(zip(header_frame2.index,header_frame2['first.seqnames'])):
        print(j)
        if(k==j):
            final.append(header[idx])
            final.append(lines[even])
            final.append(lines[even+1])
    pd.DataFrame(final).to_csv("/bucket/MillerU/Abrar/axtNet/hg38galVar1/chr{}.hg38.galVar1.net.axt".format(i),index=False,header=False, quoting=csv.QUOTE_NONE, escapechar = ' ')    
    































                   