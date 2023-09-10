import pandas as pd
import numpy as np
import sys
import matplotlib.pyplot as plt

def orthogroup_filter_list_generate():
    threshold_list = []
    deleted_clusters = []
    threshold = 50
    while threshold<150:
        df = pd.read_csv('Gene_counts_all.txt',sep="\t")
        df = df.set_index('Orthogroup')
   #     del df["Total"]
        total_cluster = df.shape[0]
        new_df = main_filter(df, threshold)
        new_cluster_count = new_df.shape[0]
        diff =  total_cluster - new_cluster_count
        threshold_list.append(threshold)
        deleted_clusters.append(diff)
        new_df.to_csv('filtered_og_' + str(threshold) + '.txt', sep = '\t', index=True)             
        threshold = threshold + 20
    plotting(threshold_list, deleted_clusters)  
        
        
def main_filter(dataf, limiter):
    df= dataf
    for index, row in df.iterrows():
        for i in row:
            if i > limiter:
                df.drop(index, inplace=True)
                break
    return(df)
                                 
def plotting(list1, list2):
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(list1, 
            list2)
    ax.set(title = "Assessment of filtering",
           xlabel = "Threshold",
           ylabel = "deleted_clusters")

    fig.savefig('assessment.png')     
                    
orthogroup_filter_list_generate()

