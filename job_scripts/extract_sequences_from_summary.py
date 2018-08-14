#!/usr/bin/env python3
'''Extract sequences of a summary table and save the sequences into a fasta file
Usage:
    ./extract_sequences_from_summary.py summary_table_file output_fasta_file 
'''

import sys

import pandas as pd


def sequence_to_fast_string(title, sequence):
    '''Convert a sequence to a fasta string.'''
    string = ''
    
    string += '>' + str(title) + '\n'
        
    start = 0
    end = min(80, len(sequence))
    string += sequence[start:end] + '\n'
    start += 80

    while start < len(sequence):
        end = min(start + 80, len(sequence))
        string += sequence[start:end] + '\n'
        start += 80

    return string

if __name__ == '__main__':
    summary_table_file = sys.argv[1]
    output_fasta_file = sys.argv[2]

    df = pd.read_csv(summary_table_file, sep='\t')

    with open(output_fasta_file, 'w') as f:
        for i in df.index:
            f.write(sequence_to_fast_string(df['design_id'][i], df['sequence'][i]))
