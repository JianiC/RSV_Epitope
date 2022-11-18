#!/usr/bin/env python
# coding: utf-8

# In[17]:


import pandas as pd
import numpy as np
import sys
import os, re
from progressbar import ProgressBar
import argparse
#import multiprocessing as mp


# In[18]:


def core_tcrf(core):
    core_tcrf = f'___{core[3]}{core[4]}{core[5]}{core[6]}{core[7]}_'
    return core_tcrf


# In[19]:


## non-cross conserved  class I eptiopes determined by TCRf 4,5,6,7,8
def epitope_distance(core1,core2,score1,score2):
    if core1[3]==core2[3] and core1[4]==core2[4] and core1[5]==core2[5] and core1[6]==core2[6] and core1[7]==core2[7]:
        return 0
    else:
        return score1+score2


# In[20]:


## example to calculate sequence distance with seq[1] and seq[2]
## i and are j are the index of seq list

def seq_dist(df_epitope,i,j):
    seq = df_epitope.ID.unique()
    seq1_epitope = df_epitope.loc[df_epitope['ID']==seq[i]]
    seq2_epitope = df_epitope.loc[df_epitope['ID']==seq[j]]
    epitope_compare = pd.merge(seq1_epitope, seq2_epitope, how="outer", on=["Pos","core_tcrf"])
    epitope_compare['EL-score_x']=epitope_compare['EL-score_x'].fillna(0)
    epitope_compare['EL-score_y']=epitope_compare['EL-score_y'].fillna(0)
    epitope_compare['core_x']=epitope_compare['core_x'].fillna('---------')
    epitope_compare['core_y']=epitope_compare['core_y'].fillna('---------')
    epitope_compare['epi_dist'] = epitope_compare.apply(lambda x: epitope_distance(str(x['core_x']),str(x['core_y']),x['EL-score_x'],x['EL-score_y']), axis =1)
    seq_distance = epitope_compare['epi_dist'].sum()
    return(seq_distance)


# In[21]:


## loop through all seqs in the eitope prediction dataframe to get the no-conserved epitope distance

def non_conserved_d(df_epitope):
    pbar = ProgressBar()
    seq = df_epitope.ID.unique()
    seq_distance=[]
    for i in pbar(range(len(seq))):
        for j in range(len(seq)):
            if i<j:
                distance = seq_dist(df_epitope,i,j)
                seq_distance.append([seq[i], seq[j], distance])                
                df_distance = pd.DataFrame(seq_distance, columns=["seq1", "seq2","epitope_distance"])           
    return(df_distance)


# In[22]:


def get_allele(allel_xls):
    with open(allel_xls, "r") as file:
        first_line = file.readline()
        allele=re.sub(r"\t|\n", "", first_line)
    return allele    


# In[24]:


def main():
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required = True)  ## input is the netMPHpan xlsx for multiple sequence but single allel
    parser.add_argument('-o', '--output', required=True) ## output is the distance matrisx for each allele, required output csv name
    args = parser.parse_args()
    print(args.input, 'is input, ', args.output, 'output.')
    
    ## input is the netMPHpan xlsx for multiple sequence but single allel
    ## for single allel, load the xls file and filter by cutoff
    data = pd.read_csv(args.input,sep="\t",skiprows=1)
    data_epitope = data.loc[data['EL_Rank'] <2] 
    data_epitope['core_tcrf'] = data_epitope.apply(lambda x:core_tcrf(x['core']),axis=1)
    # . df_epitope.to_csv(output,index=False)
    ## calculate distance
    seq_epidist=non_conserved_d(data_epitope)
    allele=get_allele(args.input)
    seq_epidist['allele']=allele
    seq_epidist.to_csv(args.output,index=False)

main()


# In[ ]:




