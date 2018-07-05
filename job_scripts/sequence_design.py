#!/usr/bin/env python3

import os
import sys
import time

import pyrosetta
from pyrosetta import rosetta

import local_protein_sequence_design as  LPSD


def print_pymol_selection_for_residues(pose, residues):
    '''Print the pymol selection command for the residues.'''
    res_commands = ['(c. {0} and res {1})'.format(pose.pdb_info().chain(i), pose.pdb_info().number(i))
                    for i in residues]

    print('sele ' + ' or '.join(res_commands))


def design(input_dir, data_path, num_jobs, job_id):
    pyrosetta.init(options='-mute all')

    # Get all the tasks

    tasks = []
    for s in sorted(os.listdir(input_dir)):
        if not s.endswith('.pdb.gz'):
            continue
        
        s_s = s.split('.')[0].split('_')
        bb_remodeled_residues = list(range(int(s_s[4]), int(s_s[5]) + 1))
        
        # Make 5 designs for each of the structure

        for i in range(5):
            tasks.append( (os.path.join(input_dir, s), bb_remodeled_residues) )

    # Make designs

    for i, t in enumerate(tasks): 
        if i % num_jobs == job_id:
            output_path = os.path.join(data_path, str(i))
            os.makedirs(output_path, exist_ok=True)

            LPSD.sequence_design.make_one_design(output_path, t[0], t[1])

if __name__ == '__main__':
    data_path = sys.argv[1]
    
    num_jobs = 1
    job_id = 0
    
    if len(sys.argv) > 3:
        num_jobs = int(sys.argv[2])
        job_id = int(sys.argv[3]) - 1

    
    start_time = time.time()
   
    #input_dir = 'test_inputs/screen_lhl_units_2lta_small'
    #input_dir = '/netapp/home/xingjiepan/Softwares/loop_helix_loop_reshaping/data/screen_lhl_units_2lta'
    input_dir = '/netapp/home/xingjiepan/Softwares/loop_helix_loop_reshaping/data/screen_lhl_units_2lv8'
    design(input_dir, data_path, num_jobs, job_id)

    end_time = time.time()
    print('Finish job in {0} seconds.'.format(int(end_time - start_time)))
