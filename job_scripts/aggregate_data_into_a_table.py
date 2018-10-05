#!/usr/bin/env python3
'''Aggregate the information for designs and filters into a summary table.
Usage:
    ./aggregate_data_into_a_table.py path_to_the_dataset
'''

import sys
import os
import json
import time

import pandas as pd


def get_dictionary_for_design(design_path, design_id):
    '''Get the information dictionary for a design.
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

    design_info['design_id'] = design_id

    return design_info


    # Make a pandas data frame

    keys = list(design_info.keys())
    values = [design_info[k] for k in keys]

    return pd.DataFrame([[design_id] + values], columns=['design_id'] + keys)

def generate_summary_table_for_dataset(path_to_the_dataset):
    '''Generate a summary table for a dataset'''
    df = None
    designs = os.listdir(path_to_the_dataset)
    start_time = time.time()
   
    design_dictionarys = []

    for i, design in enumerate(designs): 
        d_for_d = get_dictionary_for_design(os.path.join(path_to_the_dataset, design), design)

        if not (d_for_d is None):
            design_dictionarys.append(d_for_d)

        if i % 100 == 0:
            current_time = time.time()
            print('\r finish loading {0}/{1} designs in {2} seconds.'.format(i + 1, len(designs), current_time - start_time), end='')

    print('')

    df = pd.DataFrame(design_dictionarys)

    df.to_csv(os.path.join(path_to_the_dataset, 'summary_table.tsv'),
        sep='\t', index=False)
        

if __name__ == '__main__':
    path_to_the_dataset = sys.argv[1]

    generate_summary_table_for_dataset(path_to_the_dataset)
