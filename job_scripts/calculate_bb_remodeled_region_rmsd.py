#!/usr/bin/env python3
'''Calculate the backbone RMSD of the remodeled region. 
Usage:
    ./calculate_bb_remodeled_region_rmsd.py design_path1 design_path2 [design_path3 ...]
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
    
    if os.path.exists(os.path.join(design_path, 'lowest_energy_model.pdb')):
        pose_lowest_energy = rosetta.core.import_pose.pose_from_file(os.path.join(design_path, 'lowest_energy_model.pdb'))
    else:
        pose_lowest_energy = pose_design.clone()

    with open(os.path.join(design_path, 'design_info.json'), 'r') as f:
        bb_remodeled_residues_all = sorted(json.load(f)['bb_remodeled_residues'])

    # Separate remodeled residues into segments 

    separators = [i + 1 for i in range(len(bb_remodeled_residues_all) - 1) if bb_remodeled_residues_all[i] + 1 != bb_remodeled_residues_all[i + 1]] 
    separators = [0] + separators + [len(bb_remodeled_residues_all)]
   
    bb_remodeled_residues = []
    for i in range(len(separators) - 1):
        bb_remodeled_residues.append(list(range(bb_remodeled_residues_all[separators[i]], 
            bb_remodeled_residues_all[separators[i + 1] - 1] + 1)))
    
    bb_fixed_residues = [i for i in range(1, pose_design.size() + 1) if not (i in bb_remodeled_residues_all)]

    return pose_design, pose_lowest_energy, bb_remodeled_residues, bb_fixed_residues

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

def calculate_bb_remodeled_region_rmsds(design_paths):
    '''Calculate the backbone RMSDs of the remodeled region.
    Return a two dimensional list where list[i][j] element
    is the RMSD between design i and lowest scoring model j.
    ''' 
    # Update the design_paths to exclude paths that have designed structures

    design_paths = [d for d in design_paths if os.path.exists(os.path.join(d, 'design.pdb.gz'))]
   
    # Load the designs

    pose_designs = list(range(len(design_paths)))
    pose_lowest_energies = list(range(len(design_paths)))
    all_bb_remodeled_residues = list(range(len(design_paths)))
    all_bb_fixed_residues = list(range(len(design_paths)))

    for i in range(len(design_paths)):
        pose_designs[i], pose_lowest_energies[i], all_bb_remodeled_residues[i], all_bb_fixed_residues[i] = load_design(design_paths[i])

    # Loop through all pair of designs and calculate RMSDs
   
    RMSDs = []

    for i in range(len(design_paths)):
        RMSDs.append([])
        for j in range(len(design_paths)):

            # Superimpose design j to design i

            superimpose_poses_by_residues(pose_lowest_energies[j], all_bb_fixed_residues[j], pose_designs[i], all_bb_fixed_residues[i])
        
            # Find the maped residues between the two structures

            maped_residues_i = []
            maped_residues_j = []

            for k in range(len(all_bb_remodeled_residues[i])):
                seg_i = all_bb_remodeled_residues[i][k]
                seg_j = all_bb_remodeled_residues[j][k]

                min_len = min(len(seg_i), len(seg_j))
                maped_residues_i += seg_i[:min_len]
                maped_residues_j += seg_j[:min_len]

            # Calculate RMSD 

            rmsd = calc_backbone_RMSD(pose_lowest_energies[j], maped_residues_j, pose_lowest_energies[i], maped_residues_i)

            RMSDs[i].append(rmsd)
          
    return RMSDs

def plot_bb_remodeled_region_rmsds(design_paths):
    '''Calculate the backbone RMSDs of the remodeled region.'''
    # Update the design_paths to exclude paths that have designed structures

    design_paths = [d for d in design_paths if os.path.exists(os.path.join(d, 'design_info.json'))]

    rmsds = calculate_bb_remodeled_region_rmsds(design_paths)

    # Plot the heat map

    plt.close()
    
    fig, ax = plt.subplots()
    cax = ax.imshow(np.transpose(rmsds), origin='lower', interpolation='nearest', cmap='bwr', vmin=0, vmax=6)

    ## Plot RMSD values
    #for i in range(len(design_paths)):
    #    for j in range(len(design_paths)):
    #        ax.text(i - 0.35, j, "{:6.3f}".format(rmsds[i][j]), fontsize=10)

    plt.xticks(range(len(design_paths)), design_paths, fontsize=15, rotation='vertical')
    plt.yticks(range(len(design_paths)), design_paths, fontsize=15)
    
    plt.xlabel('Designs', fontsize=15)
    plt.ylabel('Lowest energy predictions', fontsize=15)

    plt.title('RMSD between designs and predicted structures', fontsize=15)

    cbar = fig.colorbar(cax)

    plt.tight_layout()  # Make sure there is room for the labels

    #plt.show()
    plt.savefig('bb_remodeled_RMSD_comparison.png')


if __name__ == '__main__':
    design_paths = sys.argv[1:]

    pyrosetta.init()

    plot_bb_remodeled_region_rmsds(design_paths)

