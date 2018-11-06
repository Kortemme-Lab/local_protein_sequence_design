#!/usr/bin/env python3
'''Superimpose the fixed regions of designs 
Usage:
    ./superimpose_fixed_regions.py data_path 
'''

import os
import sys
import json

import numpy as np
import matplotlib.pyplot as plt

import pyrosetta
from pyrosetta import rosetta


def load_design(design_path):
    '''Load a design.
    Return:
        pose_design, bb_fixed_residues
    '''
    pose_design = rosetta.core.import_pose.pose_from_file(os.path.join(design_path, 'design.pdb.gz'))
    
    with open(os.path.join(design_path, 'design_info.json'), 'r') as f:
        bb_remodeled_residues_all = sorted(json.load(f)['bb_remodeled_residues'])

    bb_fixed_residues = [i for i in range(1, pose_design.size() + 1) if not (i in bb_remodeled_residues_all)]

    return pose_design, bb_fixed_residues

def xyzV_to_np_array(xyz):
    return np.array([xyz.x, xyz.y, xyz.z])

def np_array_to_xyzV(a):
    return rosetta.numeric.xyzVector_double_t(a[0], a[1], a[2])

def np_array_to_xyzM(a):
    return rosetta.numeric.xyzMatrix_double_t.rows(
            a[0][0], a[0][1], a[0][2],
            a[1][0], a[1][1], a[1][2],
            a[2][0], a[2][1], a[2][2])

def get_superimpose_transformation(P1, P2):
    '''Get the superimpose transformation that transfoms a list of
    points P1 to another list of points P2.'''
    if len(P1) != len(P2):
        raise Exception("Sets to be superimposed must have same number of points.")

    com1 = np.mean(P1, axis=0)
    com2 = np.mean(P2, axis=0)

    R = np.dot(np.transpose(np.array(P1) - com1), np.array(P2) - com2)
    V, S, W = np.linalg.svd(R)

    if (np.linalg.det(V) * np.linalg.det(W)) < 0.0:
        V[:, -1] = -V[:, -1]

    M = np.transpose(np.array(np.dot(V, W)))

    return M, com2 - np.dot(M, com1)

def get_backbone_points(pose, residues):
    '''Get backbone points for residues in a pose.'''
    points = []

    for res in residues:
        for atom in ['N', 'CA', 'C']:
            points.append(xyzV_to_np_array(pose.residue(res).xyz(atom)))

    return points

def superimpose_poses_by_residues(pose_source, residues_source, pose_target, residues_target):
    '''Superimpose residues in a source pose into residues in a target pose.
    Only backbone atoms are used for the superimposition.
    '''
    assert(len(residues_source) == len(residues_target))

    # Get the points to be superimposed

    points_source = get_backbone_points(pose_source, residues_source)
    points_target = get_backbone_points(pose_target, residues_target)

    # Get the rigid body transformation

    M, t = get_superimpose_transformation(points_source, points_target)

    # Transform the source pose

    pose_source.apply_transform_Rx_plus_v(np_array_to_xyzM(M), 
            np_array_to_xyzV(t))

def superimpose_fixed_regions(data_path):
    '''Superimpose all fixed regions of designs in a data path.'''

    designs = os.listdir(data_path)
    designs = [d for d in designs if d.isdigit() and os.path.exists(os.path.join(data_path, d, 'design.pdb.gz'))]

    ref_pose, ref_fixed_residues = load_design(os.path.join(data_path, designs[0])) 

    for i in range(1, len(designs)):
        pose, fixed_residues = load_design(os.path.join(data_path, designs[i]))

        superimpose_poses_by_residues(pose, fixed_residues, ref_pose, ref_fixed_residues)

        pose.dump_pdb(os.path.join(data_path, designs[i], 'design.pdb.gz'))

if __name__ == '__main__':
    data_path = sys.argv[1]

    pyrosetta.init()

    superimpose_fixed_regions(data_path)
