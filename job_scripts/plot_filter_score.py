#!/usr/bin/env python3
'''Plot the distribution of a given filter score
Usage:
    ./plot_filter_score.py summary_table_file score_term
'''
import sys
import json

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def plot_filter_score(data_frame, score_term):
    '''Plot the statistics of a given filter score term'''
    # Read the scores

    scores = []

    for i in data_frame.index:
        score = data_frame[score_term][i]
        
        if not np.isnan(score):
            scores.append(score)

    # Make the plot

    median = np.median(scores)
    percentile_low = np.percentile(scores, 10)
    percentile_up = np.percentile(scores, 90)
    upper_cut = median + 1.5 * (percentile_up - percentile_low)
    lower_cut = median - 1.5 * (percentile_up - percentile_low)

    num_bins = 100
    bin_width = (upper_cut - lower_cut) / num_bins

    bins = [lower_cut + bin_width * i for i in range(num_bins)]

    hist, bin_edges = np.histogram(scores, bins=bins)

    plt.bar(bin_edges[0:-1], hist, width=bin_width)
    plt.xlabel(score_term)
    plt.ylabel('count')
    plt.show()

if __name__ == '__main__':
    summary_table_file = sys.argv[1]
    score_term = sys.argv[2]

    df = pd.read_csv(summary_table_file, sep='\t')

    plot_filter_score(df, score_term)
