from __future__ import division
import pandas as pd



plates = ['GSF2056-run1-plates1-2-demultiplexing-summary.csv',
            'GSF2056-run2-plates3-4-demultiplexing-summary.csv',
            'SampleSheet-GSF2124-run3-plates1-2.csv',
            'GSF2124 Lennon Run 3 Plates 3-4 Run Summary Sorted.csv',
            'GSF2124-run5-plates5-6-demultiplexing-summay.csv']


df = pd.read_csv('/Users/WRShoemaker/GitHub/Phylo_Evol_Timeseries/data/library_metadata/final_to_seq.txt', sep = '\t')


sample_dict = {}

for plate in plates:
    df_plate = pd.read_csv('/Users/WRShoemaker/GitHub/Phylo_Evol_Timeseries/data/library_metadata/'+plate, sep = ',')
    print(df_plate)
