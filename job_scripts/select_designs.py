#!/usr/bin/env python3
'''Select designs from a summary table. Save the selected
designs into a table file called selected_designs.tsv
Usage:
    ./select_designs.py summary_file
'''

import sys
import os

import pandas as pd

def row_path_selection(row):
    '''Return if the given row passes the selection criteria'''
    # Filter by the fragment quality analysis 
    if row['bb_remodeled_worst_fragment_crmsd'] > 1:
        return False

    # Filter by holes
    if row['movable_local_holes_score'] > 0 or row['movable_local_holes_score'] > 0:
        return False

    # Filter by helix complimentarity
    if row['bb_remodeled_residues_helix_complementarity'] < 0.6:
        return False

    # Filter by the number of buried unsatisfied hbonds
    if row['buried_unsat_for_movable_residues'] > 0:
        return False

    # Filter by the number of over saturated hbond acceptors
    if row['num_over_saturated_hbond_acceptors_for_movable_residues'] > 0:
        return False

    # Filter by the hydrophobic SASA
    if row['relative_hydrophobic_sasa_for_movable_residues'] > 0.58:
        return False

    return True


def select_designs(summary_file):
    '''Select designs from a summary table. Save the selected
    designs into a table file called selected_designs.tsv
    '''
    df_all = pd.read_csv(summary_file, sep='\t')
    df_selected = pd.DataFrame(columns=df_all.columns)

    for i, row in df_all.iterrows():
        if row_path_selection(row):
            df_selected = df_selected.append(row)

    # Save the selected designs

    dirname = os.path.dirname(summary_file)

    df_selected.to_csv(os.path.join(dirname, 'selected_designs.tsv'),
        sep='\t', index=False)


if __name__ == '__main__':
    summary_file = sys.argv[1]

    select_designs(summary_file)
