#!/usr/bin/env python3
'''Generate sequence aligment for a given set of designs. The fixed backbone region
of the designs are aligned. Then a structure based sequence aligment was performed.
Write the output to a fasta file sequence_alignment.fasta
Usage:
    ./generate_sequence_alignment_for_designs.py design_path1 design_path2 [design_path3 ...]
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
        pose_design, pose_lowest_energy, bb_remodeled_residues, bb_fixed_residues
    '''
    pose_design = rosetta.core.import_pose.pose_from_file(os.path.join(design_path, 'design.pdb.gz'))
    
    with open(os.path.join(design_path, 'design_info.json'), 'r') as f:
        bb_remodeled_residues = json.load(f)['bb_remodeled_residues']

    bb_fixed_residues = [i for i in range(1, pose_design.size() + 1) if not (i in bb_remodeled_residues)]

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

def get_target_to_source_residue_map(pose_source, pose_target):
    '''Get the residue map from the target pose to the source pose.
    Note that the length of the target design should be longer or equal
    to the source design.
    The two poses should be pre-aligned by the fixed residues.
    '''
    # Find all pairwise distances between all source residues and target residues

    s_t_distances = []

    for i in range(1, pose_source.size() + 1):
        for j in range(1, pose_target.size() + 1):
            s_t_distances.append((i, j, pose_source.residue(i).xyz('CA').distance(pose_target.residue(j).xyz('CA'))))

    s_t_distances_sorted = sorted(s_t_distances, key=lambda x : x[2])

    # Get the map from target residues to source residues

    res_map = {}

    for d in s_t_distances_sorted:
        if (not (d[1] in res_map.keys())) and (not (d[0] in res_map.values())):
            res_map[d[1]] = d[0]

        if len(res_map.keys()) == pose_source.size():
            break

    return res_map

def get_sequence_alignment_for_two_designs(design_path_source, design_path_target):
    '''Get the structure based sequence alignment for two designs.
    Note that the length of the target design should be longer or equal
    to the source design.
    Return:
        A string of aligned source sequence
    '''
    pose_source, bb_fixed_residues_source = load_design(design_path_source)
    pose_target, bb_fixed_residues_target = load_design(design_path_target)

    # Align the poses

    superimpose_poses_by_residues(pose_source, bb_fixed_residues_source, pose_target, bb_fixed_residues_target)

    # Generate the aligned sequence

    res_map = get_target_to_source_residue_map(pose_source, pose_target)
    aligned_seq = [] 

    for i in range(1, pose_target.size()):
        if i in res_map.keys():
            aligned_seq.append(pose_source.residue(res_map[i]).name1())
        else:
            aligned_seq.append('-')

    return ''.join(aligned_seq)

def generate_sequence_alignment_for_designs(design_paths):
    '''Generate structure based sequence alignment for designs.
    ''' 
    # Update the design_paths to exclude paths that have designed structures

    design_paths = [d for d in design_paths if os.path.exists(os.path.join(d, 'design.pdb.gz'))]
    design_paths = [d for d in design_paths if os.path.exists(os.path.join(d, 'design_info.json'))]
 
    # Find the design with the maximum sequence length

    max_len_design = 0
    max_len = 0
    for dp in design_paths:
        with open(os.path.join(dp, 'design_info.json'), 'r') as f:
            seq_len = len(json.load(f)['sequence'])
            
            if seq_len > max_len:
                max_len = seq_len
                max_len_design = dp

    # Get the aligned seque

    aligned_sequences = {}

    for dp in design_paths:
        aligned_sequences[dp] = get_sequence_alignment_for_two_designs(dp, max_len_design)

    # Write a fast file

    with open('sequence_alignment.fasta', 'w') as f:
        for dp in design_paths:
            f.write('>{0}\n'.format(dp))
            f.write(aligned_sequences[dp] + '\n')


if __name__ == '__main__':
    design_paths = sys.argv[1:]

    pyrosetta.init()

    generate_sequence_alignment_for_designs(design_paths)
