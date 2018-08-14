#!/usr/bin/env python3
'''Print the statistics of designs
Usage:
    ./print_design_statistics.py summary_table_file
'''
import sys
import json

import pandas as pd


def print_aa_composition(data_frame, only_count_designable_residues=True):
    '''Print amino acid composition'''
    # Cound amino acids

    aa_counts = {}

    for i in data_frame.index:
        sequence = data_frame['sequence'][i]
        
        if only_count_designable_residues:
            designable_residues = json.loads(data_frame['designable_residues'][i])
            sequence = [sequence[j - 1] for j in designable_residues ]
        
        for aa in sequence:
            if aa in aa_counts.keys():
                aa_counts[aa] += 1
            else:
                aa_counts[aa] = 1

    total_aas = sum(aa_counts[aa] for aa in aa_counts.keys())

    # Print the results

    if only_count_designable_residues:
        print("AA composition for designable residues:")
    else:
        print("AA composition:")

    keys = sorted(aa_counts.keys(), key=lambda x : aa_counts[x], reverse=True)

    for k in keys:
        print('{0}: {1:.2f}%'.format(k, 100 * aa_counts[k] / total_aas))


if __name__ == '__main__':
    summary_table_file = sys.argv[1]

    df = pd.read_csv(summary_table_file, sep='\t')

    print_aa_composition(df, only_count_designable_residues=False)
