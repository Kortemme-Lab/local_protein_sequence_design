import os
import json
import time

import pyrosetta
from pyrosetta import rosetta

from local_protein_sequence_design.basic import *

def find_surrounding_seqposes_noGP(pose, central_residue_ids, cutoff_distance=10, pre_moved_bb_pose=None):
    '''Return the residue ids that surround a given list of
    central residues and the pose of unmoved backbone.
    The selected residues are within a cutoff CA distance and
    the CA-CB vectors are pointing to the central residues.
    GLYs and PROs are ignored
    '''
    rest_of_residues = [i for i in range(1, pose.size() + 1) if not i in central_residue_ids]
    surrounding_residues = set()
    
    central_residues = [pose.residue(i) for i in central_residue_ids]
    if not (pre_moved_bb_pose is None):
        for i in range(1, pre_moved_bb_pose.size() + 1):
            central_residues.append(pre_moved_bb_pose.residue(i))

    for res1_id in rest_of_residues:
        if pose.residue(res1_id).name3() in ['GLY', 'PRO']: continue 
        ca1 = pose.residue(res1_id).xyz('CA')
        cb1 = pose.residue(res1_id).xyz('CB')

        for res2 in central_residues:
            nbra2 = res2.nbr_atom_xyz()

            if ca1.distance(nbra2) > cutoff_distance:
                continue

            if (cb1 - ca1).normalize().dot((nbra2 - ca1).normalize()) > 0.5:
                surrounding_residues.add(res1_id)

    return list(surrounding_residues)

def select_designable_residues(pose, bb_remodeled_residues, ignore_GP=True, pre_moved_bb_pose=None):
    '''Select residues that should be set to designable.
    Return a list of residue ids. 
    '''
    raw_designable_residues = bb_remodeled_residues + find_surrounding_seqposes_noGP(pose, bb_remodeled_residues, cutoff_distance=10, pre_moved_bb_pose=pre_moved_bb_pose)
    
    designable_residues = []
    
    for r in raw_designable_residues:
        if ignore_GP and pose.residue(r).name3() in ['GLY', 'PRO']:
            continue

        designable_residues.append(r)

    return designable_residues

def get_move_map(bb_movable_residues, sc_movable_residues, movable_jumps):
    '''Get a move map given the backbone movable residues, side chain
    movable residues and movable jumps. 
    Note that the side chains of bb_movable_residues are also movable.
    '''
    mm = rosetta.core.kinematics.MoveMap()
    mm.set_bb(False)
    mm.set_chi(False)
    mm.set_jump(False)

    for i in bb_movable_residues:
        mm.set_bb(i, True)
        mm.set_chi(i, True)

    for i in sc_movable_residues:
        mm.set_chi(i, True)
    
    for i in movable_jumps:
        mm.set_jump(i, True)
    
    return mm

def fast_design(pose, bb_remodeled_residues, flex_bb=True, pre_moved_bb_pose=None, do_ex_rot_run=True, sequence_symmetry_map=None):
    '''Do fast design
    Return:
        designable_residues_all, repackable_residues
    '''
    
    rosetta.basic.options.set_boolean_option('relax:constrain_relax_to_start_coords', True)
    xmlobj = rosetta.protocols.rosetta_scripts.XmlObjects.create_from_string(
    '''
    <MOVERS>
        <FastDesign name="fastdes" clear_designable_residues="0" repeats="1" ramp_down_constraints="1"/>
        <RotamerTrialsMover name="rot_trial" />
    </MOVERS>
    ''')
    fast_design = xmlobj.get_mover('fastdes')
    rot_trial = xmlobj.get_mover('rot_trial')

    # Find designable and repackable residues
   
    designable_residues = select_designable_residues(pose, bb_remodeled_residues, pre_moved_bb_pose=pre_moved_bb_pose)
    repackable_residues = find_surrounding_seqposes_noGP(pose, designable_residues, cutoff_distance=8)

    # Set score function

    sfxn = rosetta.core.scoring.get_score_function()
    if sequence_symmetry_map:
        set_non_crystal_symmetry(pose, sfxn, sequence_symmetry_map)

    fast_design.set_scorefxn(sfxn)
    rot_trial.score_function(sfxn)

    # Design everything

    task_factory = get_task_factory(pose, designable_residues, repackable_residues, extra_rotamers=False)
    
    if flex_bb:
        move_map = get_move_map([i for i in range(1, pose.size() + 1)], [], [])
    else:
        move_map = get_move_map([], designable_residues + repackable_residues, [])
    
    fast_design.set_task_factory(task_factory)
    fast_design.set_movemap(move_map)
    fast_design.apply(pose)
    
    rot_trial.task_factory(task_factory)
    for i in range(3):
        rot_trial.apply(pose)

    # Do the final design with extra rotamers
   
    if do_ex_rot_run: 

        task_factory_ex_rot = get_task_factory(pose, designable_residues, repackable_residues, extra_rotamers=True)
        fast_design.set_task_factory(task_factory_ex_rot)
        fast_design.apply(pose)

        rot_trial.task_factory(task_factory_ex_rot)
        for i in range(3):
            rot_trial.apply(pose)

    return designable_residues, repackable_residues

def make_one_design(output_path, input_pdb, bb_remodeled_residues, pre_moved_bb_pose=None, do_ex_rot_run=True, sequence_symmetry_map=None):
    '''Make one design and dump the relative information.
    Args:
        output_path: path for the outputs
        input_pdb: path to the input pdb file
        bb_remodeled_residues: a list of sequence positions for backbone
            remodeled residues
    '''
    start_time = time.time()

    # Design
    
    pose = rosetta.core.pose.Pose()
    rosetta.core.import_pose.pose_from_file(pose, input_pdb)

    designable_residues, repackable_residues = fast_design(pose, bb_remodeled_residues, pre_moved_bb_pose=pre_moved_bb_pose, 
            do_ex_rot_run=do_ex_rot_run, sequence_symmetry_map=sequence_symmetry_map)

    # Dump information

    pose.dump_pdb(os.path.join(output_path, 'design.pdb.gz'))

    end_time = time.time()
    
    info_dict = {
            'sequence' : pose.sequence(),
            'bb_remodeled_residues' : bb_remodeled_residues,
            'designable_residues' : designable_residues,
            'repackable_residues' : repackable_residues,
            'score' :  pose.energies().total_energy(), 
            'design_time' : end_time - start_time,
            }
    
    if sequence_symmetry_map:
        info_dict['sequence_symmetry_map'] = sequence_symmetry_map

    with open(os.path.join(output_path, 'design_info.json'), 'w') as f:
        json.dump(info_dict, f)
