#!/usr/bin/env python3

# python imports
import collections
import json
import numpy as np
import re
import sys

# rosetta imports
import pyrosetta
from pyrosetta import rosetta

# module imports
# import loop_helix_loop_reshaping as LHLR


def classify_loop_orientation(pose, loop_dict):
    """determine loop orientation, as defined in the publication:
    Lin, Y.-R. et al. Control over overall shape and size in de novo designed proteins.
    Proceedings of the National Academy of Sciences 112, E5478–E5485 (2015).

    Args:
        pose: rosetta pose object
        loop_dict: dictionary containing the key value pairs as defined in the fxn classify_loop_geometry

    Returns:
        character P(arallel) or A(ntiparallel) if the loop is between a sheet and helix
        character R(ight handed) or L(eft handed) if the loop is between two sheets
        character U(nclassified) if geometry does not fit within definitions of paper

    Raises:
        IndexError: alpha helix is less than 4 amino acids
    """
    unit_type = ''.join([loop_dict['start_ss'], loop_dict['stop_ss']])

    if unit_type == 'EH':
        c_alpha_xyz = pose.residue(loop_dict['start'] - 1).xyz('CA')
        c_beta_xyz = pose.residue(loop_dict['start'] - 1).xyz('CB')  # why not mention gly?

        # avg first 11 helix backbone atoms - why not multiple of 3?
        if loop_dict['stop_ss_length'] < 4:
            raise IndexError('Cannot lookup backbone atoms if secondary structure element is less than 4 amino acids')
        helix_backbone_xyz_list = []
        for helix_position in range(loop_dict['stop'] + 1, loop_dict['stop'] + 5):
            for backbone_atom in ['N', 'CA', 'C']:
                helix_backbone_xyz_list.append(pose.residue(helix_position).xyz(backbone_atom))
        del helix_backbone_xyz_list[-1]
        helix_avg_backbone = np.average(np.array(helix_backbone_xyz_list), axis=0)
        # helix_avg_backbone_xyz = rosetta.numeric.xyzVector_double_t(helix_avg_backbone)
        helix_avg_backbone_xyz = rosetta.numeric.xyzVector_double_t().assign(*helix_avg_backbone)
        # from rosetta documentation: Given two vectors (p1->p2 & p3->p4)
        angle_degrees = rosetta.numeric.angle_degrees(c_alpha_xyz, c_beta_xyz, c_alpha_xyz, helix_avg_backbone_xyz)

        if angle_degrees <= 80:
            orientation = 'P'
        elif angle_degrees >= 100:
            orientation = 'A'
        else:
            orientation = 'U'

    elif unit_type == 'HE':
        c_alpha_xyz = pose.residue(loop_dict['stop'] + 1).xyz('CA')
        c_beta_xyz = pose.residue(loop_dict['stop'] + 1).xyz('CB')

        # avg last 11 helix backbone atoms
        if loop_dict['start_ss_length'] < 4:
            raise IndexError('Cannot lookup backbone atoms if secondary structure element is less than 4 amino acids')
        helix_backbone_xyz_list = []
        for helix_position in range(loop_dict['start'] - 5, loop_dict['start'] - 1):
            for backbone_atom in ['N', 'CA', 'C']:
                helix_backbone_xyz_list.append(pose.residue(helix_position).xyz(backbone_atom))
        del helix_backbone_xyz_list[0]

        helix_avg_backbone = np.average(np.array(helix_backbone_xyz_list), axis=0)
        helix_avg_backbone_xyz = rosetta.numeric.xyzVector_double_t().assign(*helix_avg_backbone)

        # from rosetta documentation: Given two vectors (p1->p2 & p3->p4)
        angle_degrees = rosetta.numeric.angle_degrees(c_alpha_xyz, c_beta_xyz, c_alpha_xyz, helix_avg_backbone_xyz)

        if angle_degrees <= 80:
            orientation = 'P'
        elif angle_degrees >= 100:
            orientation = 'A'
        else:
            orientation = 'U'

    elif unit_type == 'EE':
        start_sheet_xyz_vector = pose.residue(loop_dict['start'] - 1).xyz('C') \
                                 - pose.residue(loop_dict['start'] - 1).xyz('N')
        start_to_stop_xyz_vector = pose.residue(loop_dict['stop'] + 1).xyz('CA')\
                                   - pose.residue(loop_dict['start'] - 1).xyz('CA')
        c_alpha_xyz = pose.residue(loop_dict['stop'] + 1).xyz('CA')
        c_beta_xyz = pose.residue(loop_dict['stop'] + 1).xyz('CB')
        vector_othogonal_to_sheets = rosetta.numeric.cross(start_sheet_xyz_vector, start_to_stop_xyz_vector)
        direction = rosetta.numeric.dot(vector_othogonal_to_sheets, c_beta_xyz - c_alpha_xyz)
        if direction > 0:
            orientation = 'R'
        else:
            orientation = 'L'

    else:
        print('No json file written: loops can only be between two sheets or a helix and a sheet - '
              'two helices is not supported')
        return False

    return orientation


def check_geometry_compliance(loop_unit, loop_orientation, loop_abego_str):
    """use the 4 parameters to determine if the loop unit matches of the labeled geometries in the publication
    if the loop unit matches, compliance is True. if not, compliance is False

    Args:
        loop_unit: two character string corresponding to the secondary structure elements before and after loop
        (EE, EH, HE)
        loop_orientation: single character string corresponding to loop orientation (L, R, P, A)
        loop_abego_str: string of abego characters - one character for each residue in loop

    Returns:
        boolean True or False

    Raises:

    """
    # unit -> orientation -> abego
    published_loop_geometries = {
        'EE': {
            'L': ['GG', 'EA', 'AA', 'BG', 'AAA', 'BGA', 'GGA', 'BGG', 'AAAG'],
            'R': ['BAAB', 'AGAB', 'BAAE', 'BAAGB']
            },
        'EH': {
            'P': ['AB', 'BBB'],
            'A': ['BAB', 'GBB']
        },
        'HE': {
            'P': ['GB', 'BA', 'AGB', 'GBA', 'AGBA', 'BAAB'],
            'A': ['GBB', 'BAA', 'AGBB']
        }
    }

    if loop_orientation != 'U' and loop_abego_str in published_loop_geometries[loop_unit][loop_orientation]:
        return True
    else:
        return False


def classify_loop_geometry(pose):
    """aggregate metrics for each loop in pdb file and store in a
    dictionary with keys as loop number and values as metrics:
        start: first amino acid position of loop
        stop: last amino acid position of loop
        start_ss: dssp character of secondary structure element preceding loop
        start_ss_length: length of secondary structure element preceding loop
        stop_ss: dssp character of secondary structure element after loop
        stop_ss_length: length of secondary structure element after loop
        orientation: single character designation (P, A, R, L, U) based on nomenclature in publication
        abego_str: string of ABEGO characters - one for each residue in loop
        geometry_compliance: True or False based on matching to the labeled in

    Args:
        pose: rosetta pose object

    Returns:
        dictionary with keys as loops numbered from 1 to N and values as

    Raises:
        Exception: protein does not begin and end with loop

    """
    loop_geometry_dict = collections.defaultdict(dict)
    dssp_str = rosetta.core.scoring.dssp.Dssp(pose).get_dssp_secstruct()
    abego_list = list(rosetta.core.sequence.get_abego(pose))
    print(pose.sequence())
    print(dssp_str)
    print(''.join(abego_list))

    if dssp_str[0] != 'L' or dssp_str[-1] != 'L':
        raise Exception('protein does not begin and end with loop')

    # finditer returns SRE_Match object
    loop_index_tuples = [match.span() for match in re.finditer('L+', dssp_str)]
    for i, loop_index_tuple in enumerate(loop_index_tuples):
        # ignore C and N terminus loops
        if loop_index_tuple[0] == 0 or loop_index_tuple[1] == len(dssp_str):
            continue
        start_index = loop_index_tuple[0]
        stop_index = loop_index_tuple[1]
        loop_geometry_dict[i]['start'] = start_index + 1
        loop_geometry_dict[i]['stop'] = stop_index
        loop_geometry_dict[i]['start_ss'] = dssp_str[start_index - 1]
        loop_geometry_dict[i]['stop_ss'] = dssp_str[stop_index]
        loop_geometry_dict[i]['start_ss_length'] = loop_index_tuple[0] - loop_index_tuples[i - 1][1]
        loop_geometry_dict[i]['stop_ss_length'] = loop_index_tuples[i + 1][0] - loop_index_tuple[1]

    for loop_index, loop_dict in loop_geometry_dict.items():
        orientation = classify_loop_orientation(pose, loop_dict)
        if orientation:
            loop_geometry_dict[loop_index]['orientation'] = orientation
        else:
            return False
        loop_geometry_dict[loop_index]['abego_str'] = ''.join(abego_list[loop_dict['start'] - 1: loop_dict['stop']])
        loop_geometry_dict[loop_index]['geometry_compliance'] = check_geometry_compliance(
            ''.join([loop_dict['start_ss'], loop_dict['stop_ss']]),
            loop_dict['orientation'],
            loop_dict['abego_str']
        )
    return loop_geometry_dict


if __name__ == '__main__':
    """script for aggregating structural information for loops in pdb file
    currently in development - may integrate functions into loop_helix_loop_reshaping 
    
    takes pdb files as command line argument
    writes json files containing structural the generated structural information
    """

    pyrosetta.init(options='-mute all')
    input_pdbs = sys.argv[1:]

    compliance_counter = collections.Counter()
    for input_pdb in input_pdbs:
        # adding print statements to help track down unexpected secondary structure elements
        print(input_pdb)
        pose = rosetta.core.import_pose.pose_from_file(input_pdb)
        loop_geometry_dict = classify_loop_geometry(pose)
        if loop_geometry_dict:
            with open('{0}.json'.format(input_pdb.split('/')[-1].split('.')[0]), 'w') as o:
                json.dump(loop_geometry_dict, o)
            for loop_number, loop_dict in loop_geometry_dict.items():
                loop_unit = ''.join([loop_dict['start_ss'], loop_dict['stop_ss']])
                compliance_counter[loop_unit] += loop_dict['geometry_compliance']
    with open('summary_geometry_compliance.json', 'w') as o:
        json.dump(compliance_counter, o)

    # # for debugging
    # pyrosetta.init()
    # input_pdb = '/Users/nwhoppe/PycharmProjects/loop_helix_loop_reshaping/test_inputs/2lv8_cleaned.pdb'
    #
    # pose = rosetta.core.import_pose.pose_from_file(input_pdb)
    # loop_geometry_dict = classify_loop_geometry(pose)
    # json.dump(loop_geometry_dict, '{0}.json'.format(input_pdb.split('/')[-1].split('.')[0]))
