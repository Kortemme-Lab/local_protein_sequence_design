#!/usr/bin/env python3

import os
import sys
import time

import pyrosetta
from pyrosetta import rosetta

import binding_site_sequence_design as  BSSD


def print_pymol_selection_for_residues(pose, residues):
    '''Print the pymol selection command for the residues.'''
    res_commands = ['(c. {0} and res {1})'.format(pose.pdb_info().chain(i), pose.pdb_info().number(i))
                    for i in residues]

    print('sele ' + ' or '.join(res_commands))

def design_for_test_inputs(data_path, num_jobs, job_id):
    pyrosetta.init(options='-extra_res_fa test_inputs/PCB.params')

    input_pdb = 'test_inputs/2fpc_pcb_matched_site.pdb'
    motif_residues = [120, 136, 245, 292]
    ligand_residue = -1
    bb_remodeled_residues = []
    num_outputs = 2
   
    for i in range(num_outputs):
        if i % num_jobs == job_id:
        
            output_path = os.path.join(data_path, str(i))
            os.makedirs(output_path, exist_ok=True)

            BSSD.sequence_design.make_one_design(output_path, input_pdb, motif_residues, bb_remodeled_residues, ligand_residue, ligand_jump_id=1)

def design_for_spinach_fluorogen_and_designed_lhl_units(data_path, num_jobs, job_id):
    pyrosetta.init(options='-extra_res_fa /netapp/home/xingjiepan/Projects/fuzz_ball_match_to_lhl_units/inputs/38E.params -mute all')

    match_dir = '/netapp/home/xingjiepan/Projects/fuzz_ball_match_to_lhl_units/outputs'
   
    # Get all the tasks

    tasks = []
    for d1 in sorted(os.listdir(match_dir)):
        d1_s = d1.split('_')
        bb_remodeled_residues = list(range(int(d1_s[4]), int(d1_s[5]) + 1))
        
        for d2 in sorted(os.listdir(os.path.join(match_dir, d1))):
            
            for target_pdb in os.listdir(os.path.join(match_dir, d1, d2)):
                if target_pdb.startswith('target_pose_'):
                    break
            
            motifs_str = target_pdb.split('.')[0][13:-1] 

            motif_residues = [int(i) for i in motifs_str.split('_')]

            # Make 10 designs for each of the match

            for i in range(10):
                tasks.append( (os.path.join(match_dir, d1, d2, target_pdb), motif_residues, bb_remodeled_residues) )

    # Make designs

    for i, t in enumerate(tasks): 
        if i % num_jobs == job_id:
            output_path = os.path.join(data_path, str(i))
            os.makedirs(output_path, exist_ok=True)

            pose = rosetta.core.pose.Pose()
            rosetta.core.import_pose.pose_from_file(pose, t[0])

            BSSD.sequence_design.make_one_design(output_path, t[0], t[1], t[2], pose.size(), ligand_jump_id=1)

if __name__ == '__main__':
    data_path = sys.argv[1]
    
    num_jobs = 1
    job_id = 0
    
    if len(sys.argv) > 3:
        num_jobs = int(sys.argv[2])
        job_id = int(sys.argv[3]) - 1

    
    start_time = time.time()
   
    #design_for_test_inputs(data_path, num_jobs, job_id)
    design_for_spinach_fluorogen_and_designed_lhl_units(data_path, num_jobs, job_id)

    end_time = time.time()
    print('Finish job in {0} seconds.'.format(int(end_time - start_time)))
