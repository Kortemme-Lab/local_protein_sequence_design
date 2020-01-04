#!/usr/bin/env python3
'''Count the sizes of designed regions.
Usage:
   ./count_designed_region_sizes.py summary_table.tsv 
'''
import sys
import pandas as pd
import json
import numpy as np

if __name__ == '__main__':
   
    summary_table = sys.argv[1]

    df = pd.read_csv(summary_table, sep='\t')

    lengths = []

    for dr_s in df['designable_residues']:
        dr = json.loads(dr_s)

        lengths.append(len(dr))

    print('max length = {0}'.format(np.max(lengths)))
    print('min length = {0}'.format(np.min(lengths)))
    print('median length = {0}'.format(np.median(lengths)))
