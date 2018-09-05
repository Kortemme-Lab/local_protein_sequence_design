#!/usr/bin/env python3
'''Generate DNA sequences for a given set of designs.
Usage:
    ./generate_dna_sequnces.py design_path1 [design_path2 ...]
'''

import sys
import os

import pyrosetta
from pyrosetta import rosetta

def reverse_translate(aa_seq):
    '''Reverse translate an amino acid sequence into a DNA sequence.'''
    ecoli_reverse_translate={   # (fold)
        'F': ['TTT', 'TTC'], 
        'L': ['CTG', 'TTA', 'TTG', 'CTT', 'CTC', 'CTA'], 
        'I': ['ATT', 'ATC', 'ATA'], 
        'M': ['ATG'], 
        'V': ['GTG', 'GTT', 'GTC', 'GTA'], 
        'S': ['AGC', 'TCT', 'TCC', 'TCG', 'AGT', 'TCA'], 
        'P': ['CCG', 'CCA', 'CCT', 'CCC'], 
        'T': ['ACC', 'ACG', 'ACT', 'ACA'], 
        'A': ['GCG', 'GCC', 'GCA', 'GCT'], 
        'Y': ['TAT', 'TAC'], 
        'H': ['CAT', 'CAC'], 
        'Q': ['CAG', 'CAA'], 
        'N': ['AAC', 'AAT'], 
        'K': ['AAA', 'AAG'], 
        'D': ['GAT', 'GAC'], 
        'E': ['GAA', 'GAG'], 
        'C': ['TGC', 'TGT'], 
        'W': ['TGG'], 
        'R': ['CGT', 'CGC', 'CGG', 'CGA', 'AGA', 'AGG'], 
        'G': ['GGC', 'GGT', 'GGG', 'GGA'], 
        '.': ['TAA', 'TGA', 'TAG']
    }

    dna_seq = ''
    for aa in aa_seq:
        dna_seq += ecoli_reverse_translate[aa][0]

    return dna_seq

def generate_dna_sequnces(design_paths, append_stop_codon=True, append_to_100_aa=True):
    '''Generate the design sequences. Save into
    a dna_sequences.tsv file
    '''
    # Update the design paths to remove paths that does not have the pdb file

    design_paths = [d for d in design_paths if os.path.exists(os.path.join(d, 'design.pdb.gz'))]

    dna_sequences = []

    for d in design_paths:
        pose = rosetta.core.import_pose.pose_from_file(os.path.join(d, 'design.pdb.gz'))
        aa_seq = pose.sequence()
        if append_stop_codon:
            aa_seq += '.'

        # Append .GS to get a sequence with 100 AA

        if append_to_100_aa and len(aa_seq) < 100:
            diff = 100 - len(aa_seq)
            seq_append = '.GS'* diff
            aa_seq += seq_append[:diff]

        dna_sequences.append(reverse_translate(aa_seq))

    # Write to a tsv file

    with open('dna_sequences.tsv', 'w') as f:
        f.write('gene name\tFASTA_seq\n')

        for i in range(len(design_paths)):
            f.write(design_paths[i] + '\t' + dna_sequences[i] + '\n')

if __name__ == '__main__':
    design_paths = sys.argv[1:]

    pyrosetta.init()
    generate_dna_sequnces(design_paths)



