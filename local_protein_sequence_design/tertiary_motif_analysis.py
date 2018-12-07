import os
import subprocess

from local_protein_sequence_design import site_settings


def parse_column_to_list(file_str, column, return_type=float):
    """return second column of file as list of column entries casted as return_type

    Args:
        file_str: string of file - not the file buffer itself
        column: zero-indexed integer column number
        return_type: basic data type
    Returns:
        list of return_type values

    Raises:
        ValueError: column entry cannot be casted to specified data type
        TypeError: column number passed is not an int
        IndexError: column index is out of range
    """
    value_list = []
    with open(file_str, 'r') as f:
        for line in f:
            value_list.append(return_type(line.rstrip().split()[column]))
    return value_list


def calc_tertiary_motif_scores(pdb_file):
    """Uses software from Grigoryan lab at Dartmouth to calculate tertiary motif metrics:
    https://grigoryanlab.org/terms/termanal.php

    Args:
        pdb_file: pdb formatted file, can be gzipped

    Returns:
        dsc50_list: tertiary design score for each residue
        ssc50_list: tertiary structure score for each residue

    Raises:

    """
    ss = site_settings.load_site_settings()
    python2_7 = os.path.abspath(ss['python2.7'])
    termanal = os.path.abspath(ss['termanal'])

    terminal_cmd = [python2_7, termanal, '--p', pdb_file]
    subprocess.check_call(terminal_cmd)

    if pdb_file.endswith('.gz'):
        prefix = pdb_file.split('.gz')[0]
    else:
        prefix = pdb_file

    dsc50_list = parse_column_to_list('{0}.dsc50'.format(prefix), 1)
    ssc50_list = parse_column_to_list('{0}.ssc50'.format(prefix), 1)

    return dsc50_list, ssc50_list
