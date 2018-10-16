import pyrosetta
from pyrosetta import rosetta


def get_atom_sasa(pose):
    '''Get atom SASA for a pose.
    Return a rosetta.core.id.AtomID_Map_double_t of SASA
    '''
    atom_sasa = rosetta.core.id.AtomID_Map_double_t() 
    rsd_sasa = rosetta.utility.vector1_double()

    rosetta.core.scoring.calc_per_atom_sasa(pose, atom_sasa, rsd_sasa, 1.4) 

    return atom_sasa

def get_buhs_for_each_atom(pose):
    '''Return a list of numbers of buried unsatisfied
    hbonds for each atom.
    Return:
        buried_unsat_acceptors, buried_unsat_donors each is a list of pairs (seqpos, atom_number).
    '''
    # Score the pose before finding the hydrogen bonds. Otherwise Hbond won't be found correctly
    
    sfxn = rosetta.core.scoring.get_score_function()
    sfxn(pose)

    # Get all hbonds of the structure

    hbset = rosetta.core.scoring.hbonds.HBondSet(pose, bb_only=False)
    hb_acceptors = [(hbset.hbond(i).acc_res(), hbset.hbond(i).acc_atm()) for i in range(1, hbset.nhbonds() + 1)] 
    hb_donors = [(hbset.hbond(i).don_res(), hbset.hbond(i).don_hatm()) for i in range(1, hbset.nhbonds() + 1)] 

    # Get the atom SASAs

    atom_sasa = get_atom_sasa(pose)

    # Find all the buried unsats

    buried_unsat_acceptors = []
    buried_unsat_donors = []
    burial_cutoff = 0.01

    for seqpos in range(1, pose.size() + 1):
        for a in range(1, pose.residue(seqpos).nheavyatoms() + 1):
            a_id = rosetta.core.id.AtomID(a, seqpos)
            a_t = pose.residue(seqpos).atom_type(a)
          
            # Check for donors

            if a_t.is_acceptor():

                if atom_sasa[a_id] < burial_cutoff and not ((seqpos, a) in hb_acceptors):
                    buried_unsat_acceptors.append((seqpos, a))

            # Check for acceptors

            elif a_t.is_donor() and pose.residue(seqpos).name3() != 'PRO':
                is_bur_unsat = True

                for h in range(pose.residue(seqpos).attached_H_begin(a), pose.residue(seqpos).attached_H_end(a) + 1): 
                    h_id = rosetta.core.id.AtomID(h, seqpos)
               
                    if atom_sasa[h_id] > burial_cutoff or (seqpos, h) in hb_donors:
                        is_bur_unsat = False

                if is_bur_unsat:
                    buried_unsat_donors.append((seqpos, a))

    return buried_unsat_acceptors, buried_unsat_donors

def get_backrub_ensemble_consensus_buhs_for_each_atom(pose):
    '''Get the list of numbers of ensemble consensus buried unsatisfied
    hbonds for each atom.
    '''
    br_mover = rosetta.protocols.backrub.BackrubMover()
    br_mover.set_max_atoms(4)
    
    buried_unsat_acceptors, buried_unsat_donors = get_buhs_for_each_atom(pose)

    # Update the number of segments
    
    tmp_pose = pose.clone()
    br_mover.apply(tmp_pose)

    # Iterate throught all segments

    for i in range(1, br_mover.num_segments() + 1):
       
        # For each segment, generate 5 structures
        
        for j in range(5):
            tmp_pose = pose.clone()
            
            br_mover.set_next_segment_id(i)
            br_mover.apply(tmp_pose)
            
            tmp_buried_unsat_acceptors, tmp_buried_unsat_donors = get_buhs_for_each_atom(tmp_pose)

            # Only keep the consensus buried unsats

            new_buried_unsat_acceptors = []
            for p in buried_unsat_acceptors:
                if p in tmp_buried_unsat_acceptors:
                    new_buried_unsat_acceptors.append(p)

            buried_unsat_acceptors = new_buried_unsat_acceptors

            new_buried_unsat_donors = []
            for p in buried_unsat_donors:
                if p in tmp_buried_unsat_donors:
                    new_buried_unsat_donors.append(p)

            buried_unsat_donors = new_buried_unsat_donors

    return buried_unsat_acceptors, buried_unsat_donors

def get_bb_buried_unsats(pose, buried_unsats):
    '''Get the buried unsats for only backbones.'''
    bb_buried_unsats = []

    for res, a in buried_unsats:
        if pose.residue(res).atom_name(a).strip() in ['N', 'O']:
            bb_buried_unsats.append((res, a))

    return bb_buried_unsats

def get_hbond_metrics(pose, bb_remodeled_residues, designable_residues, movable_residues):
    '''Get all hbond related metrics. Return a dictionary
    of the metrics.
    '''
    # Get the buried unsatisfied hydrogen bonds for each atom

    buried_unsat_acceptors, buried_unsat_donors = get_buhs_for_each_atom(pose)
    buried_unsats = buried_unsat_acceptors + buried_unsat_donors
    bb_buried_unsats = get_bb_buried_unsats(pose, buried_unsats)

    buried_unsat_acceptors_backrub, buried_unsat_donors_backrub = get_backrub_ensemble_consensus_buhs_for_each_atom(pose)
    buried_unsats_backrub = buried_unsat_acceptors_backrub + buried_unsat_donors_backrub
    bb_buried_unsats_backrub = get_bb_buried_unsats(pose, buried_unsats_backrub)

    # Get the metrics

    hbond_metrics = {}

    hbond_metrics['bb_buried_unsat_for_bb_remodeled_residues'] = sum(1 for x in bb_buried_unsats if x[0] in bb_remodeled_residues) 
    hbond_metrics['buried_unsat_for_designable_residues'] = sum(1 for x in buried_unsats if x[0] in designable_residues) 
    hbond_metrics['buried_unsat_for_movable_residues'] = sum(1 for x in buried_unsats if x[0] in movable_residues) 
    hbond_metrics['buried_unsat_for_all_residues'] =  len(buried_unsats)

    hbond_metrics['bb_backrub_ensemble_buhs_for_bb_remodeled_residues'] = sum(1 for x in bb_buried_unsats_backrub if x[0] in bb_remodeled_residues) 
    hbond_metrics['backrub_ensemble_consensus_buhs_for_designable_residues'] = sum(1 for x in buried_unsats_backrub if x[0] in designable_residues) 
    hbond_metrics['backrub_ensemble_consensus_buhs_for_movable_residues'] = sum(1 for x in buried_unsats_backrub if x[0] in movable_residues) 
    hbond_metrics['backrub_ensemble_consensus_buhs_for_all_residues'] =  len(buried_unsats_backrub)

    return hbond_metrics
