import pandas as pd
import numpy as np
import sys
import matplotlib.pyplot as plt
import scipy.stats as stats
#from astropy.stats import median_absolute_deviation


##ADJUST names block and read_csv

df = pd.read_csv('bad_clusters_20.txt',sep="\t",names=['Orthogroup','Cgri','Csab','Ecab','Ggal','Hsap','Mjav','Mmul','Mput','Nvis','Rfer','Sscr'])

df = df.set_index('Orthogroup')
threshold = 3
for index, row in df.iterrows():
    std= stats.zscore(row)
    for i in std:
        a = abs(i)
        if a >= threshold:
            df.drop(index, inplace=True)
            break

##ADJUST names block and read_csv
			
df1 = pd.read_csv('good_clusters_20.txt',sep="\t",names=['Orthogroup','Cgri','Csab','Ecab','Ggal','Hsap','Mjav','Mmul','Mput','Nvis','Rfer','Sscr'])
df1 = df1.set_index('Orthogroup')

print(df1.shape[0])

frames =[df, df1]
result = pd.concat(frames)
result = result.iloc[1:]
result.to_csv('Final_filtering_with20.tsv', sep = '\t', index=True)
