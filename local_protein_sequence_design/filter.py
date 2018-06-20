import os
import json

import numpy as np

import pyrosetta
from pyrosetta import rosetta

from local_protein_sequence_design import site_settings
from local_protein_sequence_design import fragment_quality_analysis
from local_protein_sequence_design import IO 
from local_protein_sequence_design.basic import *


def get_residue_selector_for_residues(residues):
    '''Get a residue selector for a given set of residues.'''
    return rosetta.core.select.residue_selector.ResidueIndexSelector(','.join(list(str(i) for i in residues)))

def residues_average_energy(pose, residues):
    '''Return the average energy of a set of residues.'''
    sfxn = rosetta.core.scoring.get_score_function()
    sfxn(pose)

    return np.mean(list(pose.energies().residue_total_energy(i) for i in residues))

def residues_max_energy(pose, residues):
    '''Return the max energy of a set of residues.'''
    sfxn = rosetta.core.scoring.get_score_function()
    sfxn(pose)

    return max(list(pose.energies().residue_total_energy(i) for i in residues))

def get_num_buried_unsatisfied_hbonds(pose, residues):
    '''Get the number of buried unsatisfied hbonds
    for each of the given set of residues.
    '''
    bupc = rosetta.protocols.simple_pose_metric_calculators.BuriedUnsatisfiedPolarsCalculator(
            'default', 'default')

    sfxn = rosetta.core.scoring.get_score_function()
    sfxn(pose)
    
    buhs_for_each_res = json.loads(bupc.get('residue_bur_unsat_polars', pose))

    return sum(buhs_for_each_res[i - 1] for i in residues)

def get_holes_score_for_residues(pose, residues):
    '''Get the holes score for a list of residues.'''
    ss = site_settings.load_site_settings()
    dalphaball = os.path.abspath(ss['dalphaball'])
    rosetta.basic.options.set_file_option('holes:dalphaball', dalphaball)
   
    hf = rosetta.protocols.simple_filters.HolesFilter()
    hf.set_residue_selector(get_residue_selector_for_residues(residues))
    return hf.compute(pose) / sum(pose.residue(i).natoms() for i in residues)

def get_helix_complementarity_score(pose, residues):
    '''Get the helix complementarity score for a list of residues.
    NOTE: there mask be at least one helix in the given region.
    '''
    hscf = rosetta.protocols.denovo_design.filters.SSShapeComplementarityFilter()
    hscf.set_residue_selector(get_residue_selector_for_residues(residues))
    
    return hscf.compute(pose)

def get_fragment_quality_scores(pose, bb_remodeled_residues):
    '''Get the fragment quality scores given
    the backbone remodeled region.

    Return:
        bb_remodeled_worst_fragment_crmsd, bb_remodeled_mean_fragment_crmsd

    NOTE: This function should only be called when the
    current working directory is the directory of the design.
    Also the design.pdb.gz file must exist in the current directory.
    '''
    # Create a fasta file
    
    sequence = pose.sequence()
    IO.sequence_to_fasta_file('design.fasta', 'design', sequence)

    # Do fragment picking
    
    ss = site_settings.load_site_settings()
    fqa = fragment_quality_analysis.FragmentQualityAnalyzer(
           ss['runpsipred_single'], ss['csblast'], ss['blastpgp'], ss['placeholder_seqs'], ss['sparksx_path'],
           ss['fragment_picker'], ss['vall'], ss['fragment_quality_analysis_weights'],
           rosetta_database=ss['rosetta_database_fragment_picking'])

    fdf = fqa.pick_fragments('design.pdb.gz', 'design.fasta', 
            '.', query_pos=bb_remodeled_residues)

    crmsds = fragment_quality_analysis.FragmentQualityAnalyzer.get_position_crmsd(fdf) 

    # Remove temporary files

    os.remove('frags.fsc.200.9mers')

    return max(crmsds), np.mean(crmsds)

def generate_filter_scores(filter_info_file, pose, designable_residues, repackable_residues, bb_remodeled_residues):
    '''Generate the scores of filters and save the scores
    into the filter_info_file in json format.
    '''
    filter_scores = {}

    # Get the average energy of designable residues

    filter_scores['designable_residues_average_energy'] = residues_average_energy(pose, designable_residues)

    # Get the max energy of designable residues

    filter_scores['designable_residues_max_energy'] = residues_max_energy(pose, designable_residues)

    # Get the number of buried unsatisfied hbonds

    filter_scores['buried_unsat_for_designable_residues'] = get_num_buried_unsatisfied_hbonds(pose, designable_residues)

    # Get the local holes score

    filter_scores['local_holes_score'] = get_holes_score_for_residues(pose, designable_residues)

    # Get the backbone remodeled residues helixc complementarity score

    filter_scores['bb_remodeled_residues_helix_complementarity'] = get_helix_complementarity_score(pose, bb_remodeled_residues)

    # Get fragment quality analysis scores for the backbone remodeled region
    
    #bb_remodeled_worst_fragment_crmsd, bb_remodeled_mean_fragment_crmsd = get_fragment_quality_scores(pose, bb_remodeled_residues)
   
    #filter_scores['bb_remodeled_worst_fragment_crmsd'] = bb_remodeled_worst_fragment_crmsd
    #filter_scores['bb_remodeled_mean_fragment_crmsd'] = bb_remodeled_mean_fragment_crmsd

    with open(filter_info_file, 'w') as f:
        json.dump(filter_scores, f)

