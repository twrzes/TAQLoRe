#!/usr/bin/env python3

# -*- coding: utf-8 -*-

# Importing modules
import os
import sys
import numpy as np
import pandas as pd
import tqdm
import gc
import csv
import argparse
import scipy
import scipy.stats

# This function is to read the file and to get the DataFrame with unweighted counts
# Input:
#   - input_file_path: Path to the file with reads and their corresponding transcripts (total_summary.txt file)
# Output: DataFrame with unweighted counts without 'Total' and 'N/A' rows and columns
def FileRead(input_file_path):
    # Reading the file to DataFrame
    to_return = pd.read_csv(input_file_path, sep="\t", header=0)
    return to_return

# This function is used to normalize the counts using TMM
# Input:
#   - input_downweighted_df: DataFrame with downweighted counts
# Output: DataFrame with normalized read counts
def TMMNormalization(input_downweighted_df):
    # Adding small number to prevent doing log(0)
    final_df = input_downweighted_df.copy()
    # Output DataFrames
    first_df = pd.DataFrame()
    geom_mean_df = pd.DataFrame()
    # Iterate on the sample names list to obtain a DataFrame of downweighted read counts divided by geometric mean of a transcript
    for i in range(final_df.shape[0]):
        # Obtain list of counts for a particular transcript
        transcript_counts = list(final_df.iloc[i,2:])
        # Obtain geometric mean
        transcript_geom_mean = scipy.stats.mstats.gmean(transcript_counts)
        # If the geometric mean of the transcript is equal 0, omit the transcript
        if transcript_geom_mean == 0:
            del transcript_counts, transcript_geom_mean
            continue
        # Creating output DataFrames
        transcript_geom_mean_df = pd.DataFrame({0: [final_df.iloc[i,0]], 1: [final_df.iloc[i,1]] , 2: [transcript_geom_mean]})
        geom_mean_df = geom_mean_df.append(transcript_geom_mean_df)
        del transcript_geom_mean_df
        # Obtain DataFrame with downweighted read counts multiplied by geometric mean of a transcript
        temp_df = pd.DataFrame([final_df.iloc[i,:]])
        temp_df.iloc[:,2:] = temp_df.iloc[:,2:]/transcript_geom_mean
        # Append the row to the final DataFrame
        first_df = first_df.append(temp_df)
        del transcript_counts, transcript_geom_mean, temp_df
    del i
    geom_mean_df.columns = ['#Transcript_splicing_pattern', 'Transcript_ID', 'geometric_mean']
    # Read sample names from column names of DataFrame
    iteration_list = list(final_df.columns)[2:]
    # Iterate on the DataFrame sample name columns with normalization factors and normalize the read counts by multiplification of each read count by median of normalization factor per sample
    sample_factor_df = pd.DataFrame()
    for i in range(len(iteration_list)):
        # Creating output DataFrame with sample median of normalized values
        norm_sample_factor_df = pd.DataFrame({0: [iteration_list[i]], 1: [scipy.median(list(first_df[iteration_list[i]]))]})
        sample_factor_df = sample_factor_df.append(norm_sample_factor_df)
        del norm_sample_factor_df
        final_df[iteration_list[i]] = final_df[iteration_list[i]]*(scipy.median(list(first_df[iteration_list[i]])))
    del i, first_df, iteration_list
    sample_factor_df.columns = ['#sample_id', 'normalized_values_median']
    return final_df, geom_mean_df, sample_factor_df

# Test paths
"""
in_path = '/tgac/workarea/group-eg/project_Capture/data/CACNA1C/July_17/CACNA1C_metagene_novel_exons/Results/filtered_downweighted_read_counts/table.downweighted_read_counts.sum_100.min_reads_per_sample_1.min_samples_2.txt'
"""

if __name__ == '__main__':
    # Import the input files using argparse module
    parser = argparse.ArgumentParser(description='The purpose of this script is to obtain the normalized read counts and normalization factors (TMM normalization) for all the reads')
    parser.add_argument("in_path", type=str, help="Path to the file with downweighted reads.")
    parser.add_argument("out_normalized_tmm_counts", type=str, help="Path to the output file with TMM-normalized read counts")
    parser.add_argument("out_geom_mean", type=str, help="Path to the output file with geometric mean of downweighted read counts for each transcript")
    parser.add_argument("out_norm_median", type=str, help="Path to the output file with median of normalized downweighted read counts for each sample")
    args = parser.parse_args()
    # Creating DataFrame from file with downweighted read counts
    input_df = FileRead(args.in_path)
    # Normalizing the downweighted counts using TMM normalization
    normalized_df, geometric_mean_df, normalized_sample_factors_df = TMMNormalization(input_df)
    del input_df
    gc.collect()
    # Writing the output files with normalized read counts and factors
    normalized_df.to_csv(args.out_normalized_tmm_counts, sep="\t", header=True, index=False, quoting=csv.QUOTE_NONE)
    geometric_mean_df.to_csv(args.out_geom_mean, sep="\t", header=True, index=False, quoting=csv.QUOTE_NONE)
    normalized_sample_factors_df.to_csv(args.out_norm_median, sep="\t", header=True, index=False, quoting=csv.QUOTE_NONE)
