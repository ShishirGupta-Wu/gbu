import pandas as pd
import numpy as np
import sys
#import matplotlib.pyplot as plt
#import scipy.stats as stats
#from astropy.stats import median_absolute_deviation

### interquartile range (IQR) filtering script

df = pd.read_csv('bad_clusters_50.txt',sep="\t",names=['Orthogroup','tconor','trange','tcruCLb','tbrubru','tvivax','tcruYc','tcruBc','tbruequ','tcruDm','tthei','tgrayi','tcrumar','tcongol','tbrugam','tcruSx'])

df = df.set_index('Orthogroup')
threshold = 3
for index, row in df.iterrows():
    a  =row.to_numpy()
    q1_a = np.percentile(a, 25)
    q3_a = np.percentile(a, 75)
    iqr = q3_a - q1_a
    threshold = (1.5*iqr)+q3_a
    for i in row:
        if i >= threshold:
            df.drop(index, inplace=True)
            break
			
df1 = pd.read_csv('good_clusters_50.txt',sep="\t",names=['Orthogroup','tconor','trange','tcruCLb','tbrubru','tvivax','tcruYc','tcruBc','tbruequ','tcruDm','tthei','tgrayi','tcrumar','tcongol','tbrugam','tcruSx'])
df1 = df1.set_index('Orthogroup')

print(df1.shape[0])

frames =[df, df1]
result = pd.concat(frames)
result.to_csv('Final_filtered_gene_count.tsv', sep = '\t', index=True)
