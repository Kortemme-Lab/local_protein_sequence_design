#!/usr/bin/env python3
'''Calculate the backbone RMSD of the remodeled region. 
Usage:
    ./calculate_bb_remodeled_region_rmsd.py design_path1 design_path2
'''

import os
import sys
import json

import numpy as np

import pyrosetta
from pyrosetta import rosetta


def load_design(design_path):
    '''Load a design.
    Return:
        pose, bb_remodeled_residues, bb_fixed_residues
    '''
    pose = rosetta.core.pose.Pose()
    rosetta.core.import_pose.pose_from_file(pose, os.path.join(design_path, 'design.pdb.gz'))
    
    with open(os.path.join(design_path, 'design_info.json'), 'r') as f:
        bb_remodeled_residues = json.load(f)['bb_remodeled_residues']

    bb_fixed_residues = [i for i in range(1, pose.size() + 1) if not (i in bb_remodeled_residues)]

    return pose, bb_remodeled_residues, bb_fixed_residues

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

def calc_backbone_RMSD(pose1, residues1, pose2, residues2):
    '''Calculate backbone RMSD between two poses for specific positions.'''
    assert(len(residues1) == len(residues2))

    def RMSD(points1, poinsts2):
        '''Calcualte RMSD between two lists of numpy points.'''
        diff = [points1[i] - poinsts2[i] for i in range(len(points1))]
        return np.sqrt(sum(np.dot(d, d) for d in diff) / len(diff))

    points1 = get_backbone_points(pose1, residues1)
    points2 = get_backbone_points(pose2, residues2)

    return RMSD(points1, points2)

def rmsd_between_segments(pose1, segment1, pose2, segment2):
    '''Calculate backbone RMSD between two segments.
    A segment is defined as (start, stop).
    '''
    # Swap the poses if the length of the second sse is longer than the first one

    if segment2[1] - segment2[0] < segment1[1] - segment1[0]:
        pose1, segment1, pose2, segment2 = pose2, segment2, pose1, segment1

    length1 = segment1[1] - segment1[0] + 1
    length_diff = (segment2[1] - segment2[0]) - (segment1[1] - segment1[0])

    # Calculate RMSDs for different windows

    rmsds = []
    residues1 = list(range(segment1[0], segment1[1] + 1))

    for shift in range(length_diff + 1):
        residues2 = list(range(segment2[0] + shift, segment2[0] + length1 + shift))
    
        rmsds.append(calc_backbone_RMSD(pose1, residues1, pose2, residues2))

    # Return the lowest RMSD

    return min(rmsds)

#def get_helix_segment_for_residues(pose, residues):
#    '''Return the helix segment within the given residues.'''
#    dssp = rosetta.core.scoring.dssp.Dssp(pose)
#    secstruct = dssp.get_dssp_secstruct()
#
#    helix_residues = [r for r in residues if secstruct[r - 1] == 'H']
#
#    return (min(helix_residues), max(helix_residues))

def calculate_bb_remodeled_region_rmsd(design_path1, design_path2):
    '''Calculate the backbone RMSD of the remodeled region.''' 
    pose1, bb_remodeled_residues1, bb_fixed_residues1 = load_design(design_path1)
    pose2, bb_remodeled_residues2, bb_fixed_residues2 = load_design(design_path2)
  
    # Superimpose design2 to design1

    superimpose_poses_by_residues(pose2, bb_fixed_residues2, pose1, bb_fixed_residues1)

    segment1 = (min(bb_remodeled_residues1), max(bb_remodeled_residues1))
    segment2 = (min(bb_remodeled_residues2), max(bb_remodeled_residues2))

    seg_rmsd = rmsd_between_segments(pose1, segment1, pose2, segment2)
    print('The bb remodeled region rmsd is', seg_rmsd) 


if __name__ == '__main__':
    design_path1 = sys.argv[1]
    design_path2 = sys.argv[2]

    pyrosetta.init()

    calculate_bb_remodeled_region_rmsd(design_path1, design_path2)

