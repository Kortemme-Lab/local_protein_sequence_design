#!/usr/bin/env python3

import json
import pandas as pd
import sys

if __name__ == "__main__":
"""quick script to score each design in a filter summary table using Gabe's logisitic regression model"""
    filter_tsv = sys.argv[1]
    log_reg_model_file = sys.argv[2]

    with open(log_reg_model_file, 'r') as f:
        log_reg_dict = json.load(f)

    filter_df = pd.read_csv(filter_tsv, sep='\t', header=0, index_col='design_id')
    output_df = filter_df.copy()
    for design_id in filter_df.index:
        row = filter_df.loc[design_id, :]
        prediction_score = 0
        for metric, metric_dict in log_reg_dict.items():
            mean = metric_dict['mean']
            scale = metric_dict['scale']
            weight = metric_dict['weight']
            prediction_score += ((row[metric] - mean) / scale) * weight
        output_df.loc[design_id, "length"] = len(row['sequence']) - 1
        output_df.loc[design_id, "logistic_regression_score"] = prediction_score

    output_df.to_csv('{0}_logreg.tsv'.format(filter_tsv.split('.tsv')[0]), sep='\t')
