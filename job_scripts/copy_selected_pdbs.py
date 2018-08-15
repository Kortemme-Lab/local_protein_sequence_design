#!/usr/bin/env python3
'''Copy the seleceted designs into a directory called selected_pdbs
Usage:
    ./copy_selected_pdbs.py selected_designs_file source_dir
'''

import sys
import os
import subprocess

import pandas as pd


if __name__ == '__main__':
    selected_designs_file = sys.argv[1]
    source_dir = sys.argv[2]

    # Get the selected designs

    df = pd.read_csv(selected_designs_file, sep='\t')
    
    selected_designs = []
    for i, row in df.iterrows():
        selected_designs.append(row['design_id'])

    # Copy the designs

    target_dir = 'selected_pdbs'
    os.makedirs(target_dir, exist_ok=True)

    for d_id in selected_designs:
        subprocess.call(['scp', os.path.join(source_dir, str(d_id), 'design.pdb.gz'), 
            os.path.join(target_dir, str(d_id) + '.pdb.gz')])
