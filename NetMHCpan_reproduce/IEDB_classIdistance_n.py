#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
import sys
from functools import reduce


# In[2]:


def epitope_distance(x4,y4,x5,y5,x6,y6,x7,y7,x8,y8,score_x,score_y):
    if x4==y4 and x5==y5 and x6==y6 and x7==y7 and x8==y8:
        return 0
    else:
        return score_x+score_y


# In[3]:


def sequence_distance(df,i,j):
    seq1=df[df["seq_num"] ==i]
    seq2=df[df["seq_num"] ==j]
    
    seq1_pep=seq1['peptide'].apply(lambda x: pd.Series(list(x)))
    seq1_merge = seq1.merge(seq1_pep, left_index=True, right_index=True)
    seq2_pep=seq2['peptide'].apply(lambda x: pd.Series(list(x)))
    seq2_merge = seq2.merge(seq2_pep, left_index=True, right_index=True)
    merged_df = pd.merge(seq1_merge, seq2_merge, left_on="start", right_on="start")
    merged_df['distance'] = merged_df.apply(lambda row : epitope_distance(row['4_x'],
                     row['4_y'], row['5_x'],row['5_y'],row['6_x'],row['6_y'],row['7_x'],row['7_y'],
                                                                     row['8_x'],row['8_y'],row['score_x'],row['score_y']), axis = 1)
    score=merged_df["distance"].sum()
    
    return(score)


# In[4]:


## wrap the sequence comparasion for a single allel to function, 
## df is the IEDB output for signle class I allele,n is the number of sequence
def seq_compare(df,n):
    data = []
   
    ## loop thorough n*n comparasion
    for i in range (1,n):
        for j in range (1,n):
            distance=sequence_distance(df,i,j)
            data.append([i, j, distance])
            output = pd.DataFrame(data, columns=["seq1", "seq2","distance"])
            
          
    return(output)       


# In[17]:


## wrap different allel summarise into a function
## df is output for MHC class I binding evaluaion from IEDB, wihich are build with multiople allels
## n is the number you want to compare, from seq num 1 to n
def classIepi_dist(df,n):
    i=0
    result = pd.DataFrame(columns=['seq1','seq2']) 
    for allele in df['allele'].unique():
        i=i+1
        df_allele = df[(df.allele == allele)]
        print("start working on allele" + allele)
        allele_out=seq_compare(df_allele,n) ## see seq_compare function, write output 1,2,3...
        print("finish the calulation for alle "+ str(i))
        
        result=pd.merge(result,allele_out,on=['seq1','seq2'], how='outer') 
    

    print(i+1)
    result['classI_distance']= result.iloc[:,2:i+1].sum(axis=1) ## need to adjust by the number of allel
    result = result.loc[:,~result.columns.duplicated()] ## deal with the error caused by duplicate column name with multiple allele
    result = result.filter(['seq1','seq2','classI_distance'], axis=1)
    ## convert to n*n matrix
    dis_matrix=result.pivot(index='seq1', columns='seq2', values='classI_distance')
    return(dis_matrix)
    


# In[6]:


def main():
    classI_data = pd.read_table(sys.argv[1]) ## input IEDB MHC class I result
    n=int(sys.argv[2]) ## number of sequence
    result=classIepi_dist(classI_data,n)
    
    result.to_csv(sys.argv[3])## result csv file name
main() 


# In[ ]:





# In[ ]:




