#!/usr/bin/env python3

# -*- coding: utf-8 -*-


# Importing modules
import os
import sys
import tqdm
import gc
import argparse
import copy


# This function generates dictionary with transcript ID and exon ID of the existing exons (in 1nt window)
# Input: Path to the file with all genomic positions of exons (with exon and transcript IDs)
# Output: Dictionary with genomic nucleotide positions as keys and transcript IDs and exon IDs as
def ParseExistingExonsFile(path_to_file):
    # Creating output dictionary
    to_return_dict = {}
    # Opening file with known exons
    f = open(path_to_file, 'r')
    # Read all lines of the file
    lines = f.readlines()
    # Closing and deleting the file handle
    f.close()
    del f
    # Iterating over lines in the file
    for i in tqdm.trange(len(lines)):
        # If the line starts with gene ID ('ENS' string)
        # This allows to omit header
        if lines[i].startswith('ENS'):
            # Stript and split the line
            line_split = lines[i].strip().split('\t')
            # Transcript ID
            transcript_id = line_split[3]
            # Exon ID
            exon_id = line_split[7]
            # Exon start
            exon_start = int(line_split[8])
            # Exon stop
            exon_stop = int(line_split[9])
            # Coding sequence start
            # Sometimes this string is empty because exon is not coding - then change start and stop to 0
            try:
                coding_start = int(line_split[13])
            except IndexError:
                coding_start = 0
            # Coding sequence stop
            try:
                coding_stop = int(line_split[14])
            except IndexError:
                coding_stop = 0
            # Iterate through the exons and create add transcript/exon IDs to the dictionary
            for j in range(exon_start, exon_stop+1):
                # Checking if key exists in the dictionary
                # If the key exists in the dictionary
                if j in to_return_dict:
                    # Append the transcript ID to list of transcript IDs
                    to_return_dict[j]['transcript_id'].append(transcript_id)
                    # Append the exon ID to list of exon IDs
                    to_return_dict[j]['exon_id'].append(exon_id)
                    # If the position is coding, append 'Yes' to coding list, else append 'No'
                    if j>= coding_start and j<=coding_stop:
                        to_return_dict[j]['coding_status'].append('Yes')
                    else:
                        to_return_dict[j]['coding_status'].append('No')
                # If key does not exist in the dictionary
                else:
                    # Create a dictionary within dictionary
                    to_return_dict[j] = {}
                    # Create a list with transcript_id
                    to_return_dict[j]['transcript_id'] = [transcript_id]
                    # Create a list with exon_id
                    to_return_dict[j]['exon_id'] = [exon_id]
                    # Create a list with coding status
                    # If the position is coding, put 'Yes' in coding list, else put 'No'
                    if j>= coding_start and j<=coding_stop:
                        to_return_dict[j]['coding_status'] = ['Yes']
                    else:
                        to_return_dict[j]['coding_status'] = ['No']
            # Delete unnecessary variables
            del j, line_split, transcript_id, exon_id, exon_start, exon_stop, coding_start, coding_stop
    # Delete unnecessary variables
    del i, lines
    # Return the output dictionary
    return to_return_dict

# This function adds the coordinates for novel exons to the dictionary from last file
# Input:
#   - path_to_file: Path to the file with all novel exons coordinates
#   - known_exons_dict_1nt: Dictionary containing 1nt-coordinates for known exons
# Output: Dictionary with added novel exons
def ParseNovelExonsFile(path_to_file, known_exons_dict_1nt):
    # Creating a 1:1 copy of dictionary with known exons
    to_return_dict = copy.deepcopy(known_exons_dict_1nt)
    # Opening file with novel exons
    f = open(path_to_file, 'r')
    # Read all lines of the file
    lines = f.readlines()
    # Closing and deleting the file handle
    f.close()
    del f
    # Iterating over lines in the file
    for i in tqdm.trange(len(lines)):
        # Stript and split the line
        line_split = lines[i].strip().split('\t')
        # Transcript ID
        transcript_id = line_split[3]
        # Exon ID
        exon_id = line_split[3]
        # Exon start
        exon_start = int(line_split[1])
        # Exon stop
        exon_stop = int(line_split[2])
        # Iterate through the exons and add transcript/exon IDs to the dictionary
        for j in range(exon_start, exon_stop+1):
            # If a coordinate is already in the dictionary, exit the script because something is very wrong!
            if j in to_return_dict.keys():
                break
                continue
            else:
                # Create a dictionary within dictionary
                to_return_dict[j] = {}
                # Create a list with transcript_id
                to_return_dict[j]['transcript_id'] = [transcript_id]
                # Create a list with exon_id
                to_return_dict[j]['exon_id'] = [exon_id]
                # Create a list with unknown coding status
                to_return_dict[j]['coding_status'] = ['NA']
        # Delete unnecessary variables
        del j, line_split, transcript_id, exon_id, exon_start, exon_stop
    # Delete unnecessary variables
    del i, lines
    # Return the output dictionary
    return to_return_dict

# This function creates the list of lines with 1-based coordinates of all the exons for meta-gene, with genomic and meta-gene coordinates
# Input:
#   - exons_dictionary: Dictionary containing all exons
#   - strand: strand of a gene
# Output: List with all rows of the file
def CreateOutputList(exons_dictionary, strand):
    # Output string
    output_list = []
    # List containing all coordinates
    if strand == '+':
        coordinates_list = sorted(list(exons_dictionary.keys()))
    elif strand == '-':
        coordinates_list = sorted(list(exons_dictionary.keys()), reverse=True)
    else:
        print('WRONG STRAND - CHECK IT!!!')
        sys.exit(1)
    # Meta-gene start
    meta_gene_start = 1
    # Meta-gene stop
    meta_gene_stop = 1
    # Genomic start
    genomic_start = coordinates_list[0]
    # Iterate over list of coordinates
    for i in tqdm.trange(len(coordinates_list)):
        # Last iteration needs special handling, therefore it is at the beginning of this loop
        if coordinates_list[i] == coordinates_list[-1]:
            # Genomic stop
            genomic_stop = coordinates_list[i]
            # Genomic stop string
            genomic_stop_str = str(genomic_stop)
            # Meta-gene stop string
            meta_gene_stop_str = str(meta_gene_stop)
            # String of semicolon-delimited trainscript IDs
            transcript_ids_str = exons_dictionary[genomic_start]['transcript_id']
            transcript_ids_str = ";".join(transcript_ids_str)
            # String of semicolon-delimited exon IDs
            exon_ids_str = exons_dictionary[genomic_start]['exon_id']
            exon_ids_str = ";".join(exon_ids_str)
            # String of semicolon-delimited coding status of the exon
            coding_status_str = exons_dictionary[genomic_start]['coding_status']
            coding_status_str = ";".join(coding_status_str)
            # Generating a line string
            if strand == '+':
                line_str = [str(meta_gene_start), meta_gene_stop_str, str(genomic_start), genomic_stop_str, transcript_ids_str, exon_ids_str, coding_status_str]
            elif strand == '-':
                line_str = [str(meta_gene_start), meta_gene_stop_str, genomic_stop_str, str(genomic_start), transcript_ids_str, exon_ids_str, coding_status_str]
            else:
                print('WRONG STRAND - CHECK IT!!!')
                sys.exit(1)
            line_str = "\t".join(line_str)
            # Append the line to the output list
            output_list.append(line_str)
            # Deleting unnecessary variables
            del genomic_stop, genomic_stop_str, meta_gene_stop_str, transcript_ids_str, exon_ids_str, coding_status_str, line_str
        else:
            # Compare dictionary for genomic start and i-th iteration
            if exons_dictionary[genomic_start] == exons_dictionary[coordinates_list[i]]:
                meta_gene_stop += 1
                continue
            # If dictionaries are different or this is the last iteration
            else:
                # Genomic stop
                genomic_stop = coordinates_list[i-1]
                # Genomic stop string
                genomic_stop_str = str(genomic_stop)
                # Meta-gene stop string
                meta_gene_stop_str = str(meta_gene_stop-1)
                # String of semicolon-delimited trainscript IDs
                transcript_ids_str = exons_dictionary[genomic_start]['transcript_id']
                transcript_ids_str = ";".join(transcript_ids_str)
                # String of semicolon-delimited exon IDs
                exon_ids_str = exons_dictionary[genomic_start]['exon_id']
                exon_ids_str = ";".join(exon_ids_str)
                # String of semicolon-delimited coding status of the exon
                coding_status_str = exons_dictionary[genomic_start]['coding_status']
                coding_status_str = ";".join(coding_status_str)
                # Generating a line string
                if strand == '+':
                    line_str = [str(meta_gene_start), meta_gene_stop_str, str(genomic_start), genomic_stop_str, transcript_ids_str, exon_ids_str, coding_status_str]
                elif strand == '-':
                    line_str = [str(meta_gene_start), meta_gene_stop_str, genomic_stop_str, str(genomic_start), transcript_ids_str, exon_ids_str, coding_status_str]
                line_str = "\t".join(line_str)
                # Append the line to the output list
                output_list.append(line_str)
                # Deleting unnecessary variables
                del genomic_stop, genomic_stop_str, meta_gene_stop_str, transcript_ids_str, exon_ids_str, coding_status_str, line_str
                # Changing iterators' values
                genomic_start = coordinates_list[i]
                meta_gene_start = meta_gene_stop
                meta_gene_stop += 1
    # Deleting unnecessary variables
    del i, coordinates_list, meta_gene_start, meta_gene_stop, genomic_start
    # Returning output list
    return output_list

# This function creates output file
# Input:
#   - path_to_output_file: Path to output file
#   - list_with_lines: List containing all lines to write to the output file
def WriteOutputFile(path_to_output_file, list_with_lines):
    # Creating a string to write to the output file
    to_write_str = "\n".join(list_with_lines)
    to_write_str = to_write_str + "\n"
    # Opening output file handle
    g = open(path_to_output_file, 'w')
    # Writing the output file
    g.write(to_write_str)
    # Closing the file and deleting unnecessary variables
    g.close()
    # Deleting all variables
    del to_write_str, g

# Testing paths
"""
in_known_exons = '/tgac/workarea/group-eg/project_Capture/data/CACNA1C/July_17/CACNA1C_metagene_novel_exons/Input/CACNA1C_gene_transcripts_exons_coding.txt'
in_novel_exons = '/tgac/workarea/group-eg/project_Capture/data/CACNA1C/July_17/CACNA1C_metagene_novel_exons/Input/novel_exons_coordinates_1nt.txt'
out_file = '/tgac/workarea/group-eg/project_Capture/data/CACNA1C/July_17/CACNA1C_metagene_novel_exons/test.txt'
"""

if __name__ == '__main__':
    # Import the input files using argparse module
    parser = argparse.ArgumentParser(description='This script generates a table with genomic and meta-gene positions of all exons. All coordinates are disjoin, meaning that you will have separate lines for different exons, transcripts and coding status.')
    parser.add_argument("in_known_exons", type=str, help='Path to the file containing known exons (downloaded from ENSEMBL v90)')
    parser.add_argument("in_novel_exons", type=str, help='Path to the file containing novel exons (BED3 format with 1-based coordinates)')
    parser.add_argument("param_strand", type=str, help='Strand of a gene')
    parser.add_argument("out_file", type=str, help='Path to the output file containing a table with meta-gene and genomic coordinates of all known and novel exons')
    args = parser.parse_args()
    # Dictionary with 1nt positions as keys, with transcript/exon IDs and coding status
    dict_exons_transcripts_1nt = ParseExistingExonsFile(args.in_known_exons)
    gc.collect()
    # Creating a dictionary containing nucleotide positions (1nt per dict key) of known and novel exons
    dict_with_novel_exons_transcripts_1nt = ParseNovelExonsFile(args.in_novel_exons, dict_exons_transcripts_1nt)
    del dict_exons_transcripts_1nt
    gc.collect()
    # List containing all rows for the final file
    output_list = CreateOutputList(dict_with_novel_exons_transcripts_1nt, args.param_strand)
    del dict_with_novel_exons_transcripts_1nt
    gc.collect()
    # Writing to the output file
    WriteOutputFile(args.out_file, output_list)
    del output_list
    gc.collect()
