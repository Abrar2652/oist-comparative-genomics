
import pandas as pd
## w0 = 60 # 3117730805 
### w1=100 # 71039592

other = 'felCat9'  # other distant genome that I want to merge with mm39
mm = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/hg38-mm39_outputs/pcs_sequences/all/hg38_mm39_all_pcs_seqs.tsv", sep='\t')
temp = pd.read_csv("/bucket/.mabuya/MillerU/Abrar/hg38-{0}_outputs/pcs_sequences/all/hg38_{0}_all_pcs_seqs.tsv".format(other), sep='\t')
cutoff = 60 ## taking only pcs in the power-law regime
mm = mm[mm['pcs.width']>=cutoff]
temp = temp[temp['pcs.width']>=cutoff]
mm = mm['first.sequence'].tolist() 
temp = temp['first.sequence'].tolist() 

############### merging overlapped sequence functions #####################
def concat(*args):
    result = ''
    for arg in args:
        result = concat2(result, arg)
    return result

def _concat(a, b):
    la = len(a)
    lb = len(b)
    for i in range(la):
        j = i
        k = 0
        while j < la and k < lb and a[j] == b[k]:
            j += 1
            k += 1
        if j == la:
            n = k
            break
    else:
        n = 0
    return a + b[n:]

def concat2(a,b):
    la = len(a)
    lb = len(b)
    for i in range(lb):
        j = i
        k = 0
        while j < lb and k < la and b[j] == a[k]:
            j += 1
            k += 1
        if j == lb:
            n = k
            break
    else:
        n = 0
    return b + a[n:]
######################################################################

from itertools import product
merged=[]
c=0
for i, j in product(mm, temp):
  if((len(concat(i,j)) < len(i+j)) and (len(concat(i,j)) < len(concat2(j,i)))):
    merged.append(concat(i,j))
    c+=1
    print(c)
    #print(concat(i,j))
  elif((len(concat2(j,i)) < len(i+j)) and (len(concat2(j,i)) < len(concat(i,j)))):
    merged.append(concat2(j,i))
    c+=1
    print(c)
    #print(concat2(j,i))
  elif((len(concat2(j,i)) < len(i+j)) and (len(concat(i,j)) < len(i+j)) and (len(concat(i,j)) == len(concat2(j,i)))):
    merged.append(concat(i,j))
    c+=1
    print(c)
    #print(concat(i,j))
    

final = list(dict.fromkeys(merged))  ### removing duplicates from the list 
pcs = [len(i) for i in final]  ## pcs length calculation
merged_dict = {'Merged Sequence': final, 'PCS Length': pcs}  
        
merged_df = pd.DataFrame(merged_dict)
merged_df.to_csv("hg38_mm39_{0}_pcs_merged_extended_L{1}.tsv".format(other,cutoff), sep='\t', index=False)

##################################################### 

serial_count=[]
for i in range(0, int(max(pcs)+5)):
    serial_count.append(pcs.count(i))
    
lmer = [x for x in range(int(max(pcs)+5))]
# dictionary of lists  
dict = {'PCS length (L)': lmer, 'Number of PCS of Length L': serial_count}  
        
df3 = pd.DataFrame(dict)

df3.to_csv("PCS_hg38-mm39-{0}_merged_extended_L{1}.csv".format(other,cutoff), index = False) 

























