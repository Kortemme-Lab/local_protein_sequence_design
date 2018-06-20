#!/usr/bin/env python3

import os
import sys
import time
import json

import pyrosetta
from pyrosetta import rosetta

import local_protein_sequence_design as  LPSD


def test_filter(data_path, num_jobs, job_id):
    cwd = os.getcwd()
    LPSD.site_settings.binding_site_sequence_design_home = cwd
    pyrosetta.init(options='-mute all')

    for i in os.listdir(data_path):
        if i.isdigit() and int(i) % num_jobs == job_id:
            design_path = os.path.join(data_path, i)
            os.chdir(design_path)

            pose = rosetta.core.pose.Pose()
            rosetta.core.import_pose.pose_from_file(pose, 'design.pdb.gz')
            filter_info_file = 'filter_info.json'
            design_info_file = 'design_info.json'

            with open(design_info_file, 'r') as f:
                design_info = json.load(f)

            LPSD.filter.generate_filter_scores(filter_info_file, pose,
                    design_info['designable_residues'], design_info['repackable_residues'], design_info['bb_remodeled_residues'])

            os.chdir(cwd)

if __name__ == '__main__':
    data_path = sys.argv[1]
    
    num_jobs = 1
    job_id = 0
    
    if len(sys.argv) > 3:
        num_jobs = int(sys.argv[2])
        job_id = int(sys.argv[3]) - 1

    
    start_time = time.time()
   
    test_filter(data_path, num_jobs, job_id)

    end_time = time.time()
    print('Finish job in {0} seconds.'.format(int(end_time - start_time)))
