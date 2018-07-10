#!/usr/bin/env python3
'''Aggregate the information for designs and filters into a summary table.
Usage:
    ./aggregate_data_into_a_table.py path_to_the_dataset
'''

import sys
import os
import json

import pandas as pd


def get_data_frame_for_design(design_path, design_id):
    '''Get the data frame for a design.
    Return None if there is missing information.
    '''
    design_info_file = os.path.join(design_path, 'design_info.json')
    filter_info_file = os.path.join(design_path, 'filter_info.json')

    # Load the information of the design

    if (not os.path.exists(design_info_file)) or (not os.path.exists(filter_info_file)):
        return None

    with open(design_info_file, 'r') as f:
        design_info = json.load(f)
    
    with open(filter_info_file, 'r') as f:
        filter_info = json.load(f)

    # Merge the information

    for k in filter_info.keys():
        design_info[k] = filter_info[k]

    # Make a pandas data frame

    keys = list(design_info.keys())
    values = [design_info[k] for k in keys]

    return pd.DataFrame([[design_id] + values], columns=['design_id'] + keys)

def generate_summary_table_for_dataset(path_to_the_dataset):
    '''Generate a summary table for a dataset'''
    df = None
    
    for design in os.listdir(path_to_the_dataset):
        df_for_design = get_data_frame_for_design(os.path.join(path_to_the_dataset, design), design)

        if not (df_for_design is None):
            if df is None:
                df = df_for_design

            else:
                df = df.append(df_for_design, ignore_index=True)

    df.to_csv(os.path.join(path_to_the_dataset, 'summary_table.tsv'),
        sep='\t', index=False)
        

if __name__ == '__main__':
    path_to_the_dataset = sys.argv[1]

    generate_summary_table_for_dataset(path_to_the_dataset)
