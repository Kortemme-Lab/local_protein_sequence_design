#!/usr/bin/env python3

import itertools
import pandas as pd
import sys

if __name__ == '__main__':
"""quick script to select designs from tsv file that are in a list file -
list is intended to be ids of starting models that have loop 
geometries that comply with the rules in Koga publication
note: sensitive to relationship bewteen design id and starting model (line 17)
"""    
    selected_design_tsv = sys.argv[1]
    loop_geometry_compliance_file = sys.argv[2]

    selected_df = pd.read_csv(selected_design_tsv, sep='\t', header=0)
    selected_df['starting_model'] = selected_df['design_id'] // 10

    loop_geometry_compliance_list = []
    with open(loop_geometry_compliance_file, 'r') as f:
        for line in f:
            starting_model = line.split('model_')[-1].split('.pdb')[0]
            loop_geometry_compliance_list.append(starting_model)

    model_len_dict = {}
    model_helix_dict = {}
    with open('classify_loop_geometry_starting_units.log', 'r') as f:
        for model_path, sequence, dssp_str, abego_str in itertools.zip_longest(*[f] * 4, fillvalue=None):
            try:
                starting_model = int(model_path.split('model_')[-1].split('.pdb')[0])
            except ValueError:
                continue
            model_len_dict[starting_model] = len(sequence.rstrip())
            model_helix_dict[starting_model] = dssp_str.count('H')

    output_df = selected_df[selected_df['starting_model'].isin(loop_geometry_compliance_list)]
    output_df['length'] = output_df['starting_model'].apply(lambda x: model_len_dict[x])
    output_df['helix_length'] = output_df['starting_model'].map(lambda x: model_helix_dict[x])
    output_df.to_csv('{0}_loop_length.tsv'.format(selected_design_tsv.split('.tsv')[0]), sep='\t', index=False)
