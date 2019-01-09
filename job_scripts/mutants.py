#!/usr/bin/env python3


# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   mutants.py
## @brief
## @author Evan H. Baugh, Johns Hopkins University

# mutate_residue is adapted from an original script by Sid Chaudhury


import pyrosetta
from pyrosetta import rosetta

#from local_protein_sequence_design.basic import get_task_factory

import sys


def mutate_residue(pose, mutant_position, mutant_aa, pack_radius=0.0):
    """
    Replaces the residue at  <mutant_position>  in  <pose>  with  <mutant_aa>
        and repack any residues within  <pack_radius>  Angstroms of the mutating
        residue's center (nbr_atom) using  <pack_scorefxn>
    note: <mutant_aa>  is the single letter name for the desired ResidueType

    example:
        mutate_residue(pose,30,A)
    See also:
        Pose
        PackRotamersMover
        MutateResidue
        pose_from_sequence
    """
    if not pose.is_fullatom():
        IOError('mutate_residue only works with fullatom poses')

    test_pose = rosetta.core.pose.Pose()
    test_pose.assign(pose)

    # create a standard scorefxn by default
    pack_scorefxn = rosetta.core.scoring.get_score_function()

    # task = get_task_factory(test_pose, designable_residues=[], repackable_residues=[], layered_design='')
    task = pyrosetta.standard_packer_task(test_pose)

    # the Vector1 of booleans (a specific object) is needed for specifying the
    #    mutation, this demonstrates another more direct method of setting
    #    PackerTask options for design
    aa_bool = rosetta.utility.vector1_bool()

    # PyRosetta uses several ways of tracking amino acids (ResidueTypes)
    # the numbers 1-20 correspond individually to the 20 proteogenic amino acids
    # aa_from_oneletter returns the integer representation of an amino acid
    #    from its one letter code
    # convert mutant_aa to its integer representation
    mutant_aa = rosetta.core.chemical.aa_from_oneletter_code(mutant_aa)

    # mutation is performed by using a PackerTask with only the mutant
    #    amino acid available during design
    # to do this, construct a Vector1 of booleans indicating which amino acid
    #    (by its numerical designation, see above) to allow
    for i in range(1, 21):
        # in Python, logical expression are evaluated with priority, thus the
        #    line below appends to aa_bool the truth (True or False) of the
        #    statement i == mutant_aa
        aa_bool.append(i == mutant_aa)

    # modify the mutating residue's assignment in the PackerTask using the
    #    Vector1 of booleans across the proteogenic amino acids
    task.nonconst_residue_task(mutant_position).restrict_absent_canonical_aas(aa_bool)

    # prevent residues from packing by setting the per-residue "options" of
    #    the PackerTask
    center = pose.residue(mutant_position).nbr_atom_xyz()
    for i in range(1, pose.total_residue() + 1):
        # only pack the mutating residue and any within the pack_radius
        if not i == mutant_position or center.distance_squared(
                test_pose.residue(i).nbr_atom_xyz()) > pack_radius**2:
            task.nonconst_residue_task(i).prevent_repacking()

    # apply the mutation and pack nearby residues
    packer = rosetta.protocols.minimization_packing.PackRotamersMover(pack_scorefxn, task)
    packer.apply(test_pose)

    return test_pose


if __name__ == '__main__':
    input_pdb = sys.argv[1]
    mutation_args = sys.argv[2:]
    pyrosetta.init()
    input_pose = rosetta.core.import_pose.pose_from_file(input_pdb)

    # mutations = [(19, 'V'), (23, 'I')]
    mutation_str = '_'.join([s.replace(',', '') for s in mutation_args])

    for mutation_arg in mutation_args:
        position, amino_acid = mutation_arg.split(',')
        input_pose = mutate_residue(input_pose, int(position), amino_acid, pack_radius=6)

    input_pose.dump_pdb('{0}_{1}.pdb.gz'.format(input_pdb.split('/')[-1].split('.pdb')[0], mutation_str))
