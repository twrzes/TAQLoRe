#!/usr/bin/env python3

# -*- coding: utf-8 -*-

# Importing modules
import os
import sys
import tqdm
import gc
import argparse

# This function assess whether the sum of reads in list is above threshold
# Input:
#   - input_list: List with floats of downweighted read counts for each sample
#   - sum_threshold: Threshold for the sum of downweighted read counts
# Output: True if sum of reads above or equal to threshold, False if otherwise
def SumReadsAboveThreshold(input_list, sum_threshold):
    # Sum all floats in list
    sum_of_reads = sum(input_list)
    # If sum is above or equal to threshold, return True
    if sum_of_reads >= float(sum_threshold):
        return True
    # If sum of reads is below threshold, return False
    else:
        return False

# This function assess whether there is at least n samples with at least k reads
# Input:
#   - input_list: List with floats of downweighted read counts for each sample
#   - read_count_threshold: Threshold for the number of reads present in at least specified number of samples
#   - sample_count_threshold: Threshold for number of samples having at least specified number of reads
# Output: True of there are at least n samples with at least k reads, False if otherwise
def NumReadsInNumSamplesAboveThreshold(input_list, read_count_threshold, sample_count_threshold):
    # Create a list with floats equal to or above read count threshold
    reads_above_read_count_threshold = [x for x in input_list if x >= float(read_count_threshold)]
    # If number of samples in this list is above sample threshold, return True
    if len(reads_above_read_count_threshold) >= float(sample_count_threshold):
        return True
    # Otherwise return False
    else:
        return False

# This function creates a list with lines being written to the output file. It also does the filtering using sub-functions
# Input:
#   - input_table_path: Path to the table with unfiltered isoforms
#   - num_reads_sum_threshold: Threshold for sum of reads in all samples
#   - num_reads_sample_threshold: Threshold for number of reads in chosen number of libraries (previous argument)
#   - num_samples_threshold: Threshold for number of samples having at least chosen number of reads (previous argument)
# Output: List with lines to write in the output file
def FilteringReads(input_table_path, num_reads_sum_threshold, num_reads_sample_threshold, num_samples_threshold):
    # Final list with lines meeting filtering criteria
    final_list_lines_to_return = []
    # Opening file handle
    f = open(input_table_path, 'r')
    # Iterating over lines in the file
    for line in f:
        # If line starts with hash, strip it and append to the list - it is a header
        if line.startswith('#'):
            final_list_lines_to_return.append(line.strip())
        # If line does not start with hash, do the filtering
        else:
            # Strip and split line by tab character
            line_sep = line.strip().split('\t')
            # Create a list of read counts for all samples as floats
            read_counts_list_floats = line_sep[2:]
            read_counts_list_floats = [float(x) for x in read_counts_list_floats]
            # Checking if sum of reads is above threshold
            sum_reads_above_threshold = SumReadsAboveThreshold(read_counts_list_floats, num_reads_sum_threshold)
            # Checking if in at least specified number of samples you have at least this number of reads
            reads_in_samples_above_threshold = NumReadsInNumSamplesAboveThreshold(read_counts_list_floats, num_reads_sample_threshold, num_samples_threshold)
            # If both thresholds are True, append stripped line to list of lines
            if sum_reads_above_threshold and reads_in_samples_above_threshold:
                final_list_lines_to_return.append(line.strip())
            # Delete unnecessary variables
            del line_sep, read_counts_list_floats, sum_reads_above_threshold, reads_in_samples_above_threshold
    # Close file handle
    f.close()
    # Delete unnecessary variables
    del f, line
    # Return final list
    return final_list_lines_to_return

# This function writes all isoforms above specified threshold (in list of lines) to the output file
# Input:
#   - list_of_lines: List of all lines to be written to output file
#   - output_file_path: Path to the output file
# Output: File with downweighted read counts for isoforms above specified thresholds
def WriteOutputFiles(list_of_lines, output_file_path):
    # Open file handle
    g = open(output_file_path, 'w')
    # Create a string of lines (join list by newline characters)
    string_to_write = '\n'.join(list_of_lines)
    string_to_write = string_to_write + '\n'
    # Write everything to output file
    g.write(string_to_write)
    # Close file handle
    g.close()
    # Delete unnecessary variables
    del g, string_to_write

# Testing arguments
"""
args = parser.parse_args(['/tgac/workarea/group-eg/project_Capture/data/CACNA1C/July_17/CACNA1C_metagene_novel_exons/Results/reads_splicing_patterns_removed_exons', '.included.txt', 'out_table_counts_all_isoforms.txt'])
"""


if __name__ == '__main__':
    # Import the input files using argparse module
    parser = argparse.ArgumentParser(description='This script can be used to filter out isoform above specified threshold. Thresholds are set for minimum number of total reads per isoform and per number of samples having at least specified number of reads')
    parser.add_argument('in_downweighted_table_path', type=str, help='Path to the file with unfiltered downweighted read counts.')
    parser.add_argument('filter_sum_reads', type=str, help='This specifies to keep isoforms having at least this number of reads across all sample.')
    parser.add_argument('filter_num_reads_in_samples', type=str, help='This specifies to keep isoforms having at least this number of reads in at least filter_num_samples samples.')
    parser.add_argument('filter_num_samples', type=str, help='This specifies to keep isoforms having at least filter_num_reads_in_samples reads in at least this number of samples.')
    parser.add_argument('out_filtered_table_path', type=str, help='Path to the output table with isoforms meeting filtering criteria.')
    args = parser.parse_args()
    # List of all lines to write in the output file (above thresholds)
    lines_to_write_list = FilteringReads(args.in_downweighted_table_path, args.filter_sum_reads, args.filter_num_reads_in_samples, args.filter_num_samples)
    # Write all elements in list to the output file
    WriteOutputFiles(lines_to_write_list, args.out_filtered_table_path)
    # Delete unnecessary variables
    del lines_to_write_list
    gc.collect()
