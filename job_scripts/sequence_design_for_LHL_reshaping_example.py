#!/usr/bin/env python3

import os
import sys
import time
import json

import pyrosetta
from pyrosetta import rosetta

import local_protein_sequence_design as  LPSD


def print_pymol_selection_for_residues(pose, residues):
    '''Print the pymol selection command for the residues.'''
    res_commands = ['(c. {0} and res {1})'.format(pose.pdb_info().chain(i), pose.pdb_info().number(i))
                    for i in residues]

    print('sele ' + ' or '.join(res_commands))

def get_bb_remodeled_residues_for_LHL_designs(file_for_insertion_points):
    '''Get the backbone remodeled residues for models made by
    LHL reshaping. The two achoring residues on the flanking
    secondary structures are included.
    '''
    with open(file_for_insertion_points, 'r') as f:
        insertion_points = json.load(f)

    bb_remodeled_residues = []

    for ip in insertion_points:
        start = ip['start'] + 1 if ip['start_ss'] == 'N' else ip['start']
        stop = ip['stop'] - 1 if ip['stop_ss'] == 'C' else ip['stop']

        for i in range(start, stop + 1):
            bb_remodeled_residues.append(i)

    return bb_remodeled_residues
def design(input_dir, data_path, pre_moved_bb_pdb, file_for_pre_moved_bb_insertion_points, num_jobs, job_id, num_seq_per_model=1, do_ex_rot_run=True,
        sequence_symmetry_map_generator=None):
    pyrosetta.init(options='-mute all')

    # Get all the tasks

    tasks = []
    for s in sorted(os.listdir(input_dir)):
        if not (s.endswith('.pdb.gz') or s.endswith('.pdb')):
            continue
       
        model_id = s.split('.')[0].split('_')[-1]
        file_for_insertion_points = os.path.join(input_dir, 'insertion_points_{0}.json'.format(model_id))

        # Get the bb remodeled residues

        bb_remodeled_residues = get_bb_remodeled_residues_for_LHL_designs(file_for_insertion_points)

        # Get the sequence symmetry map

        if sequence_symmetry_map_generator:
            sequence_symmetry_map = sequence_symmetry_map_generator(file_for_insertion_points)
        else:
            sequence_symmetry_map = None

        # Make multiple designs for each of the structure

        for i in range(num_seq_per_model):
            tasks.append( (os.path.join(input_dir, s), bb_remodeled_residues, sequence_symmetry_map) )

    # Get the pose for the pre-moved structure for the remodeled residues. 

    pre_moved_pose_whole = rosetta.core.import_pose.pose_from_file(pre_moved_bb_pdb)
    pre_moved_bb_insertion_points = get_bb_remodeled_residues_for_LHL_designs(file_for_pre_moved_bb_insertion_points)
    
    v_pre_moved_bb_insertion_points = rosetta.utility.vector1_unsigned_long()
    for p in pre_moved_bb_insertion_points:
        v_pre_moved_bb_insertion_points.append(p)
   
    pre_moved_bb_pose = rosetta.core.pose.Pose()
    rosetta.core.pose.pdbslice(pre_moved_bb_pose, pre_moved_pose_whole, v_pre_moved_bb_insertion_points)

    # Make designs

    for i, t in enumerate(tasks): 
        if i % num_jobs == job_id:
            output_path = os.path.join(data_path, str(i))
            os.makedirs(output_path, exist_ok=True)
            
            LPSD.sequence_design.make_one_design(output_path, t[0], t[1], pre_moved_bb_pose=pre_moved_bb_pose, do_ex_rot_run=do_ex_rot_run,
                    sequence_symmetry_map=t[2])

if __name__ == '__main__':
    data_path = sys.argv[1]
    
    num_jobs = 1
    job_id = 0
    
    if len(sys.argv) > 3:
        num_jobs = int(sys.argv[2])
        job_id = int(sys.argv[3]) - 1

    
    start_time = time.time()
   
    input_dir = 'test_inputs/two_lhl_units_2lv8'
    pre_moved_bb_pdb = 'test_inputs/2lv8_inputs/2lv8_cleaned.pdb'
    file_for_pre_moved_bb_insertion_points = 'test_inputs/2lv8_inputs/2lv8_insertion_points.json'
   
    design(input_dir, data_path, pre_moved_bb_pdb, file_for_pre_moved_bb_insertion_points, num_jobs, job_id, num_seq_per_model=1, do_ex_rot_run=False)

    end_time = time.time()
    print('Finish job in {0} seconds.'.format(int(end_time - start_time)))
