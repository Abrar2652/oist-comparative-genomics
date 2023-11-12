

import pandas as pd

df1=pd.read_csv('/bucket/MillerU/Abrar/hg19-mm10_outputs/filtered_pcs_hg19-mm10_all/chr1/hg19_mm10_chr1_identical_seqs.tsv', sep='\t')
df2=pd.read_csv('/bucket/MillerU/Abrar/hg19-mm10_outputs/filtered_pcs_hg19-mm10_all/chr2/hg19_mm10_chr2_identical_seqs.tsv', sep='\t')
df3=pd.read_csv('/bucket/MillerU/Abrar/hg19-mm10_outputs/filtered_pcs_hg19-mm10_all/chr3/hg19_mm10_chr3_identical_seqs.tsv', sep='\t')
df4=pd.read_csv('/bucket/MillerU/Abrar/hg19-mm10_outputs/filtered_pcs_hg19-mm10_all/chr4/hg19_mm10_chr4_identical_seqs.tsv', sep='\t')
df5=pd.read_csv('/bucket/MillerU/Abrar/hg19-mm10_outputs/filtered_pcs_hg19-mm10_all/chr5/hg19_mm10_chr5_identical_seqs.tsv', sep='\t')
df6=pd.read_csv('/bucket/MillerU/Abrar/hg19-mm10_outputs/filtered_pcs_hg19-mm10_all/chr6/hg19_mm10_chr6_identical_seqs.tsv', sep='\t')
df7=pd.read_csv('/bucket/MillerU/Abrar/hg19-mm10_outputs/filtered_pcs_hg19-mm10_all/chr7/hg19_mm10_chr7_identical_seqs.tsv', sep='\t')
df8=pd.read_csv('/bucket/MillerU/Abrar/hg19-mm10_outputs/filtered_pcs_hg19-mm10_all/chr8/hg19_mm10_chr8_identical_seqs.tsv', sep='\t')
df9=pd.read_csv('/bucket/MillerU/Abrar/hg19-mm10_outputs/filtered_pcs_hg19-mm10_all/chr9/hg19_mm10_chr9_identical_seqs.tsv', sep='\t')
df10=pd.read_csv('/bucket/MillerU/Abrar/hg19-mm10_outputs/filtered_pcs_hg19-mm10_all/chr10/hg19_mm10_chr10_identical_seqs.tsv', sep='\t')
df11=pd.read_csv('/bucket/MillerU/Abrar/hg19-mm10_outputs/filtered_pcs_hg19-mm10_all/chr11/hg19_mm10_chr11_identical_seqs.tsv', sep='\t')
df12=pd.read_csv('/bucket/MillerU/Abrar/hg19-mm10_outputs/filtered_pcs_hg19-mm10_all/chr12/hg19_mm10_chr12_identical_seqs.tsv', sep='\t')
df13=pd.read_csv('/bucket/MillerU/Abrar/hg19-mm10_outputs/filtered_pcs_hg19-mm10_all/chr13/hg19_mm10_chr13_identical_seqs.tsv', sep='\t')
df14=pd.read_csv('/bucket/MillerU/Abrar/hg19-mm10_outputs/filtered_pcs_hg19-mm10_all/chr14/hg19_mm10_chr14_identical_seqs.tsv', sep='\t')
df15=pd.read_csv('/bucket/MillerU/Abrar/hg19-mm10_outputs/filtered_pcs_hg19-mm10_all/chr15/hg19_mm10_chr15_identical_seqs.tsv', sep='\t')
df16=pd.read_csv('/bucket/MillerU/Abrar/hg19-mm10_outputs/filtered_pcs_hg19-mm10_all/chr16/hg19_mm10_chr16_identical_seqs.tsv', sep='\t')
df17=pd.read_csv('/bucket/MillerU/Abrar/hg19-mm10_outputs/filtered_pcs_hg19-mm10_all/chr17/hg19_mm10_chr17_identical_seqs.tsv', sep='\t')
df18=pd.read_csv('/bucket/MillerU/Abrar/hg19-mm10_outputs/filtered_pcs_hg19-mm10_all/chr18/hg19_mm10_chr18_identical_seqs.tsv', sep='\t')
df19=pd.read_csv('/bucket/MillerU/Abrar/hg19-mm10_outputs/filtered_pcs_hg19-mm10_all/chr19/hg19_mm10_chr19_identical_seqs.tsv', sep='\t')
df20=pd.read_csv('/bucket/MillerU/Abrar/hg19-mm10_outputs/filtered_pcs_hg19-mm10_all/chr20/hg19_mm10_chr20_identical_seqs.tsv', sep='\t')
df21=pd.read_csv('/bucket/MillerU/Abrar/hg19-mm10_outputs/filtered_pcs_hg19-mm10_all/chr21/hg19_mm10_chr21_identical_seqs.tsv', sep='\t')
df22=pd.read_csv('/bucket/MillerU/Abrar/hg19-mm10_outputs/filtered_pcs_hg19-mm10_all/chr22/hg19_mm10_chr22_identical_seqs.tsv', sep='\t')
dfX=pd.read_csv('/bucket/MillerU/Abrar/hg19-mm10_outputs/filtered_pcs_hg19-mm10_all/chrX/hg19_mm10_chrX_identical_seqs.tsv', sep='\t')
dfY=pd.read_csv('/bucket/MillerU/Abrar/hg19-mm10_outputs/filtered_pcs_hg19-mm10_all/chrY/hg19_mm10_chrY_identical_seqs.tsv', sep='\t')


i='all'
pd.concat([df1,df2,df3,df4,df5,df6,df7,df8,df9,df10,df11,df12,df13,df14,df15,df16,df17,df18,df19,df20,df21,df22,dfX,dfY],ignore_index=True).to_csv(r'/bucket/MillerU/Abrar/hg19-mm10_outputs/filtered_pcs_hg19-mm10_all/{0}/hg19_mm10_{0}_identical_seqs.tsv'.format(i), sep="\t", index=False) 





