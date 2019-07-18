#!/usr/bin/env python3

# -*- coding: utf-8 -*-

# Importing modules
import os
import sys
import numpy as np
import gc
import copy
import argparse
from natsort import natsorted


# This function creates a dictionary of all known paths of exon-exon junction
# Input: Path to the input file with mappings
# Output: Dictionary with start node as key and list of nodes this edge connects to
def KnownEdgesGraphDict(path_to_input_file):
    # Creating intermediate dictionary of ENSEMBL transcript IDs as keys and list of exons as values
    dict_transcript_exons = {}
    # Exon number
    exon_number = 0
    # Open file handle
    f = open(path_to_input_file, 'r')
    # Iterate over lines in input file
    for line in f:
        # If line does not start with hash
        if not line.startswith('#'):
            # Increment exon number
            exon_number += 1
            # Split line by tab character
            line_sep = line.strip().split('\t')
            # If it is UTR, we do not need it in the dictionary
            if line_sep[-1] == 'UTR':
                del line_sep
                continue
            # If an exon is coding
            if 'Yes' in line_sep[6]:
                # Split transcript IDs by semicolon
                transcript_id = line_sep[4].split(';')
                # Split coding status of each isoform containing exon by semicolon
                coding_status = line_sep[6].split(';')
                # For transcript in the list of transcript id
                for i in range(len(transcript_id)):
                    # If a transcript is coding, add an exon to dictionary list
                    if coding_status[i] == 'Yes':
                        # With .setdefault method, append the exon_number to the list, if list is not there because transcript_id is not a key in the dictionary, create this key with empty list and subsequently append
                        dict_transcript_exons.setdefault(transcript_id[i], []).append(str(exon_number))
                # Delete unnecessary variables
                del transcript_id, coding_status, i
            # Delete unnecessary variables
            del line_sep
    # Close file
    f.close()
    # Delete unnecessary variables
    del exon_number, f
    # Create a final dictionary to return containing all possible known paths
    dict_known_exon_paths = {}
    # Iterate over each key in dictionary
    for transcript_id in dict_transcript_exons.keys():
        # Obtain list of all junctions
        list_of_all_junctions = dict_transcript_exons[transcript_id]
        # Iterate over list of all junctions
        for i in range(len(list_of_all_junctions)-1):
            # Start node is ith element in list
            start_node = list_of_all_junctions[i]
            # End node is i+1 element in list
            end_node = list_of_all_junctions[i+1]
            # If a key is not in final dictionary, create it with empty list
            if start_node not in dict_known_exon_paths.keys():
                dict_known_exon_paths[start_node] = []
            # Append end node to the list of nodes if it is not already there
            if end_node not in dict_known_exon_paths[start_node]:
                dict_known_exon_paths[start_node].append(end_node)
            # Delete unnecessary variables
            del start_node, end_node
        # Delete unnecessary variables
        del list_of_all_junctions, i
    # Delete unnecessary variables
    del dict_transcript_exons, transcript_id
    # Return final dictionary
    return dict_known_exon_paths

# Generating a dictionary with transcript IDs as keys and average and median expression as sub-keys, as well with expression the whole line in the file
# Input: Path to file with normalized (TMM) expression
# Output: Dictionary with mean and median expression of isoforms, as well as whole line in the row of file
def AverageMedianSamplesExpressionTMMDict(path_to_file_with_expression):
    # Dictionary to return
    dict_mean_median_tmm = {}
    # String with sample names to return
    string_with_sample_names = ''
    # Open file
    f = open(path_to_file_with_expression, 'r')
    # Iterate on lines in file
    for line in f:
        # If line starts with hash
        if line.startswith('#'):
            # Split string by tab character
            line_sep = line.strip().split('\t')
            # string_with_sample_names
            string_with_sample_names = line_sep[2:]
            string_with_sample_names = '\t'.join(string_with_sample_names)
        # If line does not start with hash
        else:
            # Split line by tab character
            line_sep = line.strip().split('\t')
            # Transcript ID
            transcript_id = line_sep[1]
            # Expression list
            expression_list = line_sep[2:]
            # Convert values to floats
            expression_list = [float(x) for x in expression_list]
            # Calculate mean
            mean_tmm = str(round(np.mean(expression_list), 6))
            # Calculate median
            median_tmm = str(round(np.median(expression_list), 6))
            # Calculate sum
            sum_tmm = str(round(sum(expression_list), 6))
            del expression_list, line_sep
            # Add to dictionary
            dict_mean_median_tmm[transcript_id] = {}
            dict_mean_median_tmm[transcript_id]['average'] = mean_tmm
            dict_mean_median_tmm[transcript_id]['median'] = median_tmm
            dict_mean_median_tmm[transcript_id]['sum'] = sum_tmm
            dict_mean_median_tmm[transcript_id]['whole_line'] = line.strip()
            del transcript_id, mean_tmm, median_tmm, sum_tmm
    # Close file
    f.close()
    # Delete variables not needed anymore
    del line, f
    # Return dictionary
    return dict_mean_median_tmm, string_with_sample_names

# This function looks for all possible paths in a graph
# The code of this function has been copied from: https://www.python.org/doc/essays/graphs/
# Input:
#   - graph: Dictionary with known exon-exon junctions
#   - start: Start node
#   - end: End node
#   - path: List of path the start and end node are connected
# Output: List of lists with all possible paths between start and end node
def FindAllPathsInGraph(graph, start, end, path=[]):
    # Add starting node to path
    path = path + [start]
    # If the start and end nodes are the same, return the start node as a list
    if start == end:
        return [path]
    # If graph has no start node, return empty list
    if start not in graph.keys():
        return []
    # Paths as an empty list
    paths = []
    # For end node in start position in the known junctions dict
    for node in graph[start]:
        # If node is not in the list
        if node not in path:
            # Recursively obtain a path one step after using node as a start node and path with added start node as path
            newpaths = FindAllPathsInGraph(graph, node, end, path)
            # For each new path, append this path to the list to return
            for newpath in newpaths:
                paths.append(newpath)
    # Return list of paths
    return paths


# This function iterates on the annotation file (with coding status) and calculate missing exons/novel junctions
# Input:
#   - dict_graph_known_junctions: Dictionary with all known exon-exon junctions in the form of the graph
#   - annotation_coding_status_file_path: Path to file with annotations (together with coding status)
# Output: Dictionary with missing exons/novel exon-exon junctions
def MissingExonsNovelJunctionsDict(dict_graph_junctions, annotation_coding_status_file_path):
    # Create an output dictionary
    missing_exons_novel_junctions_dict = {}
    # Open file handle
    f = open(annotation_coding_status_file_path, 'r')
    # Iterate over lines in file
    for line in f:
        # Omit header lines
        if not line.startswith('#'):
            # Strip and split line by header
            line_sep = line.strip().split('\t')
            # If transcript is already known or not annotated as coding, skip to the next line
            if line_sep[1] != 'Novel' or line_sep[5] != 'coding':
                del line_sep
                continue
            # Transcript ID
            transcript_id = line_sep[2]
            # Splicing pattern string
            splicing_pattern_str = line_sep[0]
            # Split splicing pattern by underscore
            splicing_pattern_list = splicing_pattern_str.split('_')
            # Missing exon list
            missing_exons_list = []
            # Novel junctions list
            novel_junctions_list = []
            # Iterate over list to extract start and stop node
            for i in range(len(splicing_pattern_list)-1):
                # Start node
                start_node = splicing_pattern_list[i]
                # End node
                end_node = splicing_pattern_list[i+1]
                # Obtain a list of paths
                list_of_paths = FindAllPathsInGraph(dict_graph_junctions, start_node, end_node)
                # If list of paths is empty, it means that there is no path that exists in the known junction dictionary
                # Therefore, start and end node should be added to novel junction list
                if len(list_of_paths)==0:
                    novel_junction_string = start_node + '_' + end_node
                    novel_junctions_list.append(novel_junction_string)
                    del novel_junction_string
                else:
                    # If any sub-list in list contains start and end node than it means that path exist
                    # If not then there is at least one exon that is missing
                    if not any([len(x)==2 for x in list_of_paths]):
                        # For each sub-list in list of paths
                        # List of missing exons for a specific edge
                        specific_edge_missing_exons_list = []
                        for sublist in list_of_paths:
                            missing_exons = [x for x in sublist if x!=start_node and x!=end_node]
                            # Extend a list of missing exons
                            specific_edge_missing_exons_list.extend(missing_exons)
                            del missing_exons
                        del sublist
                        # Remove duplicates from list
                        seen = set()
                        specific_edge_missing_exons_list = [x for x in specific_edge_missing_exons_list if x not in seen and not seen.add(x)]
                        del seen
                        # Extend missing exon list
                        missing_exons_list.extend(specific_edge_missing_exons_list)
                        del specific_edge_missing_exons_list
                # Delete unnecessary variables
                del start_node, end_node, list_of_paths
            # Delete unnecessary variables
            del i, splicing_pattern_str, splicing_pattern_list
            # Sorting all lists
            missing_exons_list.sort(key=int)
            novel_junctions_list.sort()
            # If there is no missing exon
            if len(missing_exons_list)==0:
                missing_exons_str = 'N/A'
                missing_exons_num_str = '0'
            # Else, join missing exons by semicolon and count them
            else:
                missing_exons_num_str = str(len(missing_exons_list))
                missing_exons_str = ';'.join(missing_exons_list)
            # If there is no novel junction
            if len(novel_junctions_list)==0:
                novel_junctions_str = 'N/A'
                novel_junctions_num_str = '0'
            # Else, join missing exons by semicolon and count them
            else:
                novel_junctions_num_str = str(len(novel_junctions_list))
                novel_junctions_str = ';'.join(novel_junctions_list)
            if len(missing_exons_list)==0 and len(novel_junctions_list)==0:
                novel_combination = 'novel_combination'
            else:
                novel_combination = 'N/A'
            # Delete unnecessary variables
            del missing_exons_list, novel_junctions_list
            # If transcript ID is already in the dictionary, exit script
            if transcript_id in missing_exons_novel_junctions_dict.keys():
                print('TRANSCRIPT ID ALREADY IN THE DICTIONARY!!!')
                print(transcript_id)
                sys.exit(1)
            # Add missing exons and novel junctions to the dictionary
            missing_exons_novel_junctions_dict[transcript_id] = {}
            missing_exons_novel_junctions_dict[transcript_id]['missing_exons_num'] = missing_exons_num_str
            missing_exons_novel_junctions_dict[transcript_id]['missing_exons_names'] = missing_exons_str
            missing_exons_novel_junctions_dict[transcript_id]['novel_junctions_num'] = novel_junctions_num_str
            missing_exons_novel_junctions_dict[transcript_id]['novel_junctions_names'] = novel_junctions_str
            missing_exons_novel_junctions_dict[transcript_id]['novel_combination'] = novel_combination
            # Delete unnecessary variables
            del line_sep, transcript_id, missing_exons_num_str, missing_exons_str, novel_junctions_str, novel_junctions_num_str, novel_combination
    # Close the file handle
    f.close()
    # Delete unnecessary variables
    del f, line
    # Return the dictionary
    return missing_exons_novel_junctions_dict

# Create a list with lines to write in the output file
# Columns in the file are:
#   splice pattern, transcript ID, number of missing exons, names of missing exons, number of novel splice patterns, names of novel splice patterns, information if a splice pattern is a novel combination of existing exon paths, average expression, median expression, sum expression, expression values per sample
# Input:
#   - expression_dict: Dictionary with expression values
# Output: List of lines to write in the output table
def LinesOutputTableList(expression_dict, exons_junct_dict):
    # List of lines to output
    list_of_all_lines = []
    # Create a list of all possible transcript ID names
    list_of_transcript_ids = list(exons_junct_dict.keys())
    # Sort the list by the transcript number
    list_of_transcript_ids = natsorted(list_of_transcript_ids)
    # Iterate over list of transcript IDs
    for transcript_id in list_of_transcript_ids:
        # Obtain a line from expression dictionary
        line_expr_dict = expression_dict[transcript_id]['whole_line']
        # Split the line by tab character
        line_expr_dict_split_list = line_expr_dict.split('\t')
        # Splice pattern is the first element in the list
        splice_pattern = line_expr_dict_split_list[0]
        # Expression values for all samples are from column 2
        expression_all_samples = line_expr_dict_split_list[2:]
        # Join the expression of all samples by tab character
        expression_all_samples = '\t'.join(expression_all_samples)
        # Create a list of all columns per line
        list_of_all_columns_per_line = [splice_pattern, transcript_id, exons_junct_dict[transcript_id]['missing_exons_num'], exons_junct_dict[transcript_id]['missing_exons_names'], exons_junct_dict[transcript_id]['novel_junctions_num'], exons_junct_dict[transcript_id]['novel_junctions_names'], exons_junct_dict[transcript_id]['novel_combination'], expression_dict[transcript_id]['average'], expression_dict[transcript_id]['median'], expression_dict[transcript_id]['sum'], expression_all_samples]
        # Join all elements in this list per tab character
        list_of_all_columns_per_line = '\t'.join(list_of_all_columns_per_line)
        # Append the line to the list of lines
        list_of_all_lines.append(list_of_all_columns_per_line)
        # Delete unnecessary variables
        del line_expr_dict, line_expr_dict_split_list, splice_pattern, expression_all_samples, list_of_all_columns_per_line
    # Delete unnecessary variables
    del list_of_transcript_ids, transcript_id
    # Return list of lines
    return list_of_all_lines

# This function generates an output file
# Input:
#   - input_list: List with all the strings to write in the file
#   - header_bool: True if write a header, False if not
#   - header_string: String to write as header
#   - path_to_output_file: Path to output file
# Output: Output file
def WriteOutputFile(input_list, header_bool, header_string, path_to_output_file):
    # Opening output file handle
    g = open(path_to_output_file, 'w')
    # If there is header, write header string to output file
    if header_bool:
        g.write(header_string)
    # Join list by newlines, add newline at the end of file
    string_to_write = '\n'.join(input_list)
    string_to_write = string_to_write + '\n'
    # Writing contents of output file
    g.write(string_to_write)
    g.close()
    del g, string_to_write

# Test paths
"""
args = parser.parse_args(['/tgac/workarea/group-eg/project_Capture/data/CACNA1C/July_17/CACNA1C_metagene_novel_exons/Meta_gene_files/meta_gene_genomic_exon_coordinates_UTR.txt', '/tgac/workarea/group-eg/project_Capture/data/CACNA1C/July_17/CACNA1C_metagene_novel_exons/Results/tmm_normalized_filtered_downweighted_read_counts/table.tmm_norm.sum_100.min_reads_per_sample_1.min_samples_2.txt', '/tgac/workarea/group-eg/project_Capture/data/CACNA1C/July_17/CACNA1C_metagene_novel_exons/Results/coding_length_extra_annotation/annotation_table.sum_100.min_reads_per_sample_1.min_samples_2.coding_length.txt'])
"""

if __name__ == '__main__':
    # Import the input files using argparse module
    parser = argparse.ArgumentParser(description='This script is used to generate a table with annotation of exon skipping events in coding transcripts (and their expression measures).')
    parser.add_argument("in_mappings", type=str, help="Path to the file with exonic positions in genome and meta-gene.")
    parser.add_argument("in_expression_tmm", type=str, help="Path to the file with TMM-normalized expression values for all isoforms.")
    parser.add_argument("in_annotations", type=str, help="Path to the file with all splice patterns and their annotations (including coding status).")
    parser.add_argument("out_info_table", type=str, help="Path to table with information about missing exons/novel junctions.")
    args = parser.parse_args()
    # Dictionary with paths of exon-exon junctions
    dict_known_junctions = KnownEdgesGraphDict(args.in_mappings)
    # Creating dictionary with average and median expression of each isoforms
    dict_mean_median_all_samples_expression, string_with_sample_names = AverageMedianSamplesExpressionTMMDict(args.in_expression_tmm)
    # Dictionary with missing exons/novel junctions
    dict_exon_miss_novel_junct = MissingExonsNovelJunctionsDict(dict_known_junctions, args.in_annotations)
    del dict_known_junctions
    gc.collect()
    # Create list of lines to output file
    list_of_lines_output_table = LinesOutputTableList(dict_mean_median_all_samples_expression, dict_exon_miss_novel_junct)
    del dict_mean_median_all_samples_expression, dict_exon_miss_novel_junct
    gc.collect()
    # Write output table
    WriteOutputFile(list_of_lines_output_table, True, str('#Splicing_pattern\tTranscript_ID\tNum_of_missing_exons\tMissing_exon_names\tNum_of_novel_junctions\tNovel_junction_names\tNew_combination_of_junctions\tAverage_expression\tMedian_expression\tSum_expression\t' + string_with_sample_names + '\n'), args.out_info_table)
    del list_of_lines_output_table
    gc.collect()
