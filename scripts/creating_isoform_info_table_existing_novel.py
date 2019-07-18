#!/usr/bin/env python3

# -*- coding: utf-8 -*-

# Importing modules
import os
import sys
import tqdm
import gc
import argparse

# This function iterates over a file with all mappings and creates two dictionaries - dictionary with all exon numbers with novel exons and dictionary with splice patterns and ENSEMBL Transcript IDs (ENST)
# Input: Path to file with all mappings
# Output:
#   - dict_novel_exons_to_return: Dictionary with exon numbers (equal to line numbers in file with mappings) as keys and novel exon names as values
def NovelExonsSplicePatternsDict(path_to_file_with_mappings):
    # Creating dictionary with novel exons
    dict_novel_exons_to_return = {}
    # Opening file handle
    f = open(path_to_file_with_mappings, 'r')
    # Creating line number variable
    line_num = 0
    # Iterating over lines in file
    for line in f:
        # If line does not start with hash
        if not line.startswith('#'):
            # Increment line number
            line_num += 1
            # Strip and split lines by tab character
            line_sep = line.strip().split('\t')
            # If last column is 'UTR', then continue to next iteration
            if line_sep[-1] == 'UTR':
                del line_sep
                continue
            # Obtain list of transcript IDs
            # In this line you also have novel exon IDs
            transcript_novel_exons = line_sep[4]
            # If this line contains 'novel_exon' string it means that it is a novel exon and it should be added to a dictionary of novel exon, under a line_num key
            if 'novel_exon' in transcript_novel_exons:
                dict_novel_exons_to_return[str(line_num)] = transcript_novel_exons
            # Delete unnecessary variables
            del line_sep, transcript_novel_exons
    # Close the file handle
    f.close()
    # Delete unnecessary variables
    del f, line_num, line
    return dict_novel_exons_to_return

# This function iterates over lines in downweighted read counts and creates a list with lines to output file
# The whole output file will be a table with following lines:
# | #Splicing_pattern | Annotated_novel |  Transcript_ID  | Novel_exon_status | Novel_exon_IDs |
# |:-----------------:|:---------------:|:---------------:|:-----------------:|:--------------:|
# |     1_2_3_4_5     |    Annotated    |    ENST012345   |         No        |       N/A      |
# |    2_3_100_112    |      Novel      | CACNA1C_novel-1 |        Yes        |  novel_exon123 |
# Input:
#   - new_exons_dict: Dictionary with line numbers in mapping file containing novel exons
#   - path_to_filtered_downweighted_read_counts: Path to the file with filtered downweighted read counts
# Output: List with lines to write to the output file
def SplicePatternInfoTable(new_exons_dict, path_to_filtered_downweighted_read_counts):
    # Creating output list
    output_lines_list_to_return = []
    # Opening file handle
    f = open(path_to_filtered_downweighted_read_counts, 'r')
    # Iterating over lines in file
    for line in f:
        # Omit header line in file
        if not line.startswith('#'):
            # Strip and split lines by tab character
            line_sep = line.strip().split('\t')
            # Splice pattern is at the beginning of the line
            splice_pattern = line_sep[0]
            # Transcript ID is the second column in the line
            transcript_id = line_sep[1]
            # If transcript ID contains 'ENST' then it is an annotated transcript
            if 'ENS' in transcript_id:
                annotated_novel_status = 'Annotated'
            # Else if there is a 'CACNA1C_novel' string in transcript ID, the transcript is novel
            elif '_novel' in transcript_id:
                annotated_novel_status = 'Novel'
            # Else quit the script because you did something wrong!
            else:
                print(transcript_id)
                print('WRONG TRANSCRIPT ID IN THE SECOND COLUMN OF THE FILE WITH DOWNWEIGHTED COUNTS!')
                sys.exit(1)
            # Split splice pattern by underscore
            splice_pattern_split = splice_pattern.split('_')
            # Novel exons numeric IDs present in splice pattern
            novel_exons_in_splice_pattern = [x for x in new_exons_dict.keys() if x in splice_pattern_split]
            # If length of this list is 0 which means that there are no novel exons in the transcript, set novel exon status as 'No' and novel exon IDs as 'N/A'
            if len(novel_exons_in_splice_pattern) == 0:
                novel_exon_status = 'No'
                novel_exon_ids = 'N/A'
            # If there is at least one novel exon, obtain its` name (semicolon-delimited) and change novel exon status to 'Yes'
            else:
                novel_exon_status = 'Yes'
                novel_exon_ids = [new_exons_dict[key] for key in novel_exons_in_splice_pattern]
                novel_exon_ids = ';'.join(novel_exon_ids)
            # Delete unnecessary variables
            del line_sep, splice_pattern_split, novel_exons_in_splice_pattern
            # Construct list with all lines
            output_line = [splice_pattern, annotated_novel_status, transcript_id, novel_exon_status, novel_exon_ids]
            output_line = '\t'.join(output_line)
            # Append line to the list of lines
            output_lines_list_to_return.append(output_line)
            # Delete variables not needed anymore
            del splice_pattern, annotated_novel_status, transcript_id, novel_exon_status, novel_exon_ids, output_line
    # Close the file handle
    f.close()
    # Delete unnecessary variables
    del f, line
    # Return list of lines
    return output_lines_list_to_return

# This function writes output file
# Input:
#   - list_of_lines: List of all lines to be written to output file
#   - output_file_path: Path to the output file
# Output: Output file
def WriteOutputFiles(list_of_lines, output_file_path):
    # Open file handle
    g = open(output_file_path, 'w')
    # Create a header string to write
    header_string_to_write = '#Splicing_pattern\tAnnotated_novel\tTranscript_ID\tNovel_exon_status\tNovel_exon_IDs\n'
    # Write header to output file
    g.write(header_string_to_write)
    del header_string_to_write
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
args = parser.parse_args(['Meta_gene_files/meta_gene_genomic_exon_coordinates.txt', 'Results/filtered_downweighted_read_counts/table.downweighted_read_counts.sum_100.min_reads_per_sample_1.min_samples_2.txt'])
"""

if __name__ == '__main__':
    # Import the input files using argparse module
    parser = argparse.ArgumentParser(description='The purpose of this script is to generate a table with information about each splicing pattern, i.e. annotated/novel transcript, name of isoform (if known; if not known then the name will be CACNA1C_novel-number), presence/absence of novel exons and IDs of novel exons (if present; if not then N/A).')
    parser.add_argument('in_mappings', type=str, help='Path to the file with 1-based positions of each exon (with separated coding/non-coding parts of the exons to another lines in file) in genome, metagene, transcript IDs, exon IDs and coding status of the sequence.')
    parser.add_argument('in_downweighted_filtered_counts', type=str, help='Path to the file with filtered downweighted read counts.')
    parser.add_argument('out_annotation', type=str, help='Path to the output file with desired information about all splicing patterns.')
    args = parser.parse_args()
    # Dictionaries with novel exons and splicing patterns of known transcripts
    novel_exon_dict = NovelExonsSplicePatternsDict(args.in_mappings)
    # List with all lines to the output file
    out_lines_list = SplicePatternInfoTable(novel_exon_dict, args.in_downweighted_filtered_counts)
    del novel_exon_dict
    gc.collect()
    # Writing output table with all information about splice patterns
    WriteOutputFiles(out_lines_list, args.out_annotation)
    del out_lines_list
    gc.collect()
