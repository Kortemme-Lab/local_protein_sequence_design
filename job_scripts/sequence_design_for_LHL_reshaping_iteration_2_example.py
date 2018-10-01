#!/usr/bin/env python3

import os
import sys
import time
import json

import pandas as pd

import pyrosetta
from pyrosetta import rosetta

import local_protein_sequence_design as  LPSD


def print_pymol_selection_for_residues(pose, residues):
    '''Print the pymol selection command for the residues.'''
    res_commands = ['(c. {0} and res {1})'.format(pose.pdb_info().chain(i), pose.pdb_info().number(i))
                    for i in residues]

    print('sele ' + ' or '.join(res_commands))

def design(input_dir, data_path, selected_designs_file, num_jobs, job_id, num_seq_per_model=1, do_ex_rot_run=True):
    pyrosetta.init(options='-mute all')

    # Get all the tasks

    df_selected = pd.read_csv(selected_designs_file, sep='\t')
    tasks = []

    for i, row in df_selected.iterrows():
        for j in range(num_seq_per_model):
            tasks.append(row['design_id'])

    # Make designs

    for i, t in enumerate(tasks): 
        if i % num_jobs == job_id:
            output_path = os.path.join(data_path, str(i))
            os.makedirs(output_path, exist_ok=True)
          
            input_pdb = os.path.join(input_dir, str(t), 'design.pdb.gz')
            
            with open(os.path.join(input_dir, str(t), 'design_info.json'), 'r') as f:
                input_design_info = json.load(f)

            if not ('sequence_symmetry_map' in input_design_info.keys()):
                input_design_info['sequence_symmetry_map'] = None

            LPSD.sequence_design.make_one_design(output_path, input_pdb, input_design_info['bb_remodeled_residues'],
                    designable_residues=input_design_info['designable_residues'],repackable_residues=input_design_info['repackable_residues'],
                    do_ex_rot_run=do_ex_rot_run, sequence_symmetry_map=input_design_info['sequence_symmetry_map'])

if __name__ == '__main__':
    data_path = sys.argv[1]
    
    num_jobs = 1
    job_id = 0
    
    if len(sys.argv) > 3:
        num_jobs = int(sys.argv[2])
        job_id = int(sys.argv[3]) - 1

    
    start_time = time.time()
   
    input_dir = 'data/sequence_design_for_LHL_reshaping_example'
    selected_designs_file = 'data/sequence_design_for_LHL_reshaping_example/selected_designs_for_next_iteration.tsv'

    design(input_dir, data_path, selected_designs_file, num_jobs, job_id, num_seq_per_model=1, do_ex_rot_run=True)

    end_time = time.time()
    print('Finish job in {0} seconds.'.format(int(end_time - start_time)))
