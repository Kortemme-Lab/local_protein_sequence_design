#!/usr/bin/env python3

import numpy as np
import pandas as pd
from sklearn.cluster import AgglomerativeClustering
from scipy.spatial.distance import squareform
import sys

import pyrosetta
from pyrosetta import rosetta

from calculate_bb_remodeled_region_rmsd import calculate_bb_remodeled_region_rmsds


def hierarchical_cluster(distance_matrix, n_clusters=10):
    """cluster distance matrix hierarchically
    :param distance_matrix: n x n symmetrical matrix w/ zero diagonal, entry Mij is distance between row i and col j
    :return: cluster classification for each entry 1 to n
    """

    # convert distance matrix into vector
    distance_vector = squareform(np.array(distance_matrix))
    hierarchical_clustering = AgglomerativeClustering(n_clusters=n_clusters, affinity='precomputed').fit(distance_vector)
    print(hierarchical_clustering.labels_)
    print(len(hierarchical_clustering.labels_))
    return hierarchical_clustering.labels_


if __name__ == '__main__':
    unique_design_tsv = sys.argv[1]
    data_path_prefix = '/netapp/home/nwhoppe/local_protein_sequence_design/data/sequence_design_for_LHL_reshaping_iter1_EHEE_rd1_0284_layer15/'

    design_df = pd.read_csv(unique_design_tsv, sep='\t', header=0, index_col='design_id')
    design_path_list = [data_path_prefix + str(design_id) for design_id in design_df.index.values]

    pyrosetta.init()
    rmsd_matrix = calculate_bb_remodeled_region_rmsds(design_path_list)
    cluster_labels = hierarchical_cluster(rmsd_matrix)
    for design_id in design_df.index:
        cluster_id = cluster_labels[design_path_list.index(data_path_prefix + str(design_id))]
        design_df.loc[design_id, 'cluster_id'] = cluster_id

    design_df.to_csv('_clusters.tsv'.format(unique_design_tsv.split('.tsv')[0]), sep='\t')
