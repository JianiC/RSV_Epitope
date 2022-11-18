#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pandas as pd
import numpy as np
from sklearn.manifold import MDS
from matplotlib import pyplot as plt
import sklearn.datasets as dt
import seaborn as sns         
from sklearn.metrics.pairwise import manhattan_distances, euclidean_distances
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
## k means to find the best number of group
from sklearn.preprocessing import MinMaxScaler
from sklearn.cluster import KMeans
import argparse


# In[ ]:


## from distance dataframe to generate a full matrix
def fullmatrix(df1):
    df2 = df1 
    data=pd.concat([df1,df2.rename(columns={'seq1':'seq2','seq2':'seq1'})], ignore_index=True) ## concate but switch seq1 and seq2
    #data= data.filter(['seq1','seq2','epitope_distance'])
    data= data.pivot_table(index='seq1',columns='seq2',values='epitope_distance')
    data=data.fillna(0)
    return(data)


# In[ ]:


## print MDS use user define group
## print strain and coordinate
def mdsplot_usergroup(strain, coords,meta):
    df=pd.DataFrame({'strain':strain, 'mds1':coords[:, 0], 'mds2': coords[:, 1]})
    df_meta = pd.merge(df, meta, on='strain')
    ## plot
    fig, ax = plt.subplots(figsize=(10, 10))
    ax = sns.scatterplot(data = df_meta, x="mds1", y="mds2", hue="group", cmap='rainbow')
    fig.savefig('mds_usergroup.png')   # save the figure to file
    plt.close(fig)    # close the figure window
    


# In[ ]:


def main():
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required = True)  ## input is the netMPHpan xlsx for multiple sequence but single allel
    parser.add_argument("-m", "--meta",required= False,
                    help= "user defined csv file")
    #parser.add_argument('-o', '--output', required=True) ## output is the mds coordiant
    args = parser.parse_args()
    print(args.input, 'is input, ', args.output, 'output.', args.meta, "user defined groups")
    ## get data that generated a full matrix from distance dataframe
    data = fullmatrix(data1)
    ## perform MDS
    mds = MDS(n_components=2, dissimilarity="precomputed", random_state=0)## create mds module
    results = mds.fit(data)
    coords = results.embedding_
    strain = data.columns
    meta = pd.read_csv(args.meta)
    if args.meta:
        mdsplot_usergroup(strain,coords,meta)
    


main()

