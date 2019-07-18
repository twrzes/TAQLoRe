#!/usr/bin/env python3

# -*- coding: utf-8 -*-

# Importing modules
import os
import sys
import pybedtools as pb
import numpy as np
import tqdm
import gc
import argparse

# This function creates a dictionary with novel exon name as key and start and end as value
# Input: Path to the BED file with 'broad' and 'narrow' positions of novel exons
# Output: Dictionary with novel exon name as key and start and end as value
def NovelExonStartStopDict(path_to_narrow_novel_exons_file):
    # Creating dictionary to return
    dict_novel_exons_start_stop = {}
    # Open file handle
    f = open(path_to_narrow_novel_exons_file, 'r')
    # Iterate over lines in a file
    for line in f:
        # Strip and split line by tab character
        line_sep = line.strip().split('\t')
        # Obtain novel exon name
        novel_exon_id = line_sep[4]
        # If 9th or 10th column has question mark, use third and fourth as start and stop, respectively
        if any([line_sep[8]=='?', line_sep[9]=='?']):
            start_pos = line_sep[2]
            stop_pos = line_sep[3]
        else:
            start_pos = line_sep[8]
            stop_pos = line_sep[9]
        # Add to the dictionary
        dict_novel_exons_start_stop[novel_exon_id] = {}
        dict_novel_exons_start_stop[novel_exon_id]['start'] = start_pos
        dict_novel_exons_start_stop[novel_exon_id]['stop'] = stop_pos
        # Remove unnecessary variables
        del line_sep, novel_exon_id, start_pos, stop_pos
    # Close file
    f.close()
    # Remove unnecessary variables
    del f, line
    # Return dictionary with novel exons
    return dict_novel_exons_start_stop

# This function creates a dictionary (1-based) of all mappings
# Input:
#   - path_to_input_file: Path to the file with all mappings
#   - path_to_genome_fasta: Path to genomic FASTA file
#   - chrom: Chromosome (must match chromosome in FASTA file)
#   - dictionary_novel_exons_start_stop: Dictionary with start and stop positions of novel exons
# Output: Dictionary of following structure:
# Dictionary:
#   - key: string with exon number (row in file with mappings)
#       - key: 'start'
#           - value: start position (0-based)
#       - key: 'stop'
#           - value: stop position
#       - key: 'coding_status'
#           - value: 'novel_exon', 'Yes', 'No'
#       - key: 'exon_name'
#           - value: 'novel_exon1, 'novel_exon2', etc.; or 'coding_line_num1', 'non_coding_line_num10', etc.
#       - key: 'sequence'
#           - value: sequence of exon
def MappingsDictionary(path_to_input_file, path_to_genome_fasta, chrom, dictionary_novel_exons_start_stop):
    # Creating dictionary to return
    dict_all_mappings = {}
    # Exon number
    exon_number = 0
    # Open file handle
    f = open(path_to_input_file, 'r')
    # Iterate over lines in input file
    tq = tqdm.tqdm(total=os.path.getsize(path_to_input_file), unit='B', unit_scale=True)
    for line in f:
        tq.update(len(line))
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
            # Make a sub-dictionary
            dict_all_mappings[str(exon_number)] = {}
            # If an exon is a novel one
            if 'novel_exon' in line_sep[4]:
                # Coding status
                coding_status = 'novel_exon'
                # Exon name
                exon_name = line_sep[4]
                # Start position
                start_pos = dictionary_novel_exons_start_stop[line_sep[4]]['start']
                # Stop position
                stop_pos = dictionary_novel_exons_start_stop[line_sep[4]]['stop']
            # If an exon is not a novel exon
            else:
                # Start position
                start_pos = line_sep[2]
                # Stop position
                stop_pos = line_sep[3]
                # If coding status is No
                if not 'Yes' in line_sep[6]:
                    coding_status = 'No'
                    exon_name = 'non_coding_line_num' + str(exon_number)
                # If coding status is Yes
                else:
                    coding_status = 'Yes'
                    exon_name = 'coding_line_num' + str(exon_number)
            # Obtain a space separated string with chromosome, start (0-based) and end
            bed_string = [chrom, str(int(start_pos)-1), stop_pos]
            bed_string = ' '.join(bed_string)
            # Convert coordinates to BedTool
            a = pb.BedTool(bed_string, from_string=True)
            del bed_string
            # Extract sequence for a particular region
            genomic_sequence = a.sequence(fi=path_to_genome_fasta)
            del a
            # Opening a file with sequence
            genomic_sequence_string = open(genomic_sequence.seqfn).read()
            del genomic_sequence
            # Split a sequence by newline character
            genomic_sequence_list = genomic_sequence_string.strip().split('\n')
            del genomic_sequence_string
            # Remove FASTA headers from sequences
            genomic_sequence_list = [x for x in genomic_sequence_list if not x.startswith('>')]
            # Sanity check - must be of length 1
            if len(genomic_sequence_list) != 1:
                print('FASTA HAS MORE THAN ONE SEQUENCE - CHECK IT!!!')
                print(exon_number)
                sys.exit(1)
            # Join all elements of list by nothing
            genomic_sequence_fasta = ''.join(genomic_sequence_list)
            del genomic_sequence_list
            # Convert to capitals
            genomic_sequence_fasta = genomic_sequence_fasta.upper()
            # Add start position to dictionary
            dict_all_mappings[str(exon_number)]['start'] = start_pos
            # Add stop position to dictionary
            dict_all_mappings[str(exon_number)]['stop'] = stop_pos
            # Adding exon name to dictionary
            dict_all_mappings[str(exon_number)]['exon_name'] = exon_name
            # Add coding status to dictionary
            dict_all_mappings[str(exon_number)]['coding_status'] = coding_status
            # Add sequence to the dictionary
            dict_all_mappings[str(exon_number)]['sequence'] = genomic_sequence_fasta
            # Delete unnecessary variables
            del start_pos, stop_pos, exon_name, coding_status, genomic_sequence_fasta
    # Close file handle
    f.close()
    # Delete unnecessary variables
    del f, exon_number, tq, line
    # Return dictionary
    return dict_all_mappings

# Generating a dictionary with transcript IDs as keys and average and median expression as sub-keys
# Input: Path to file with normalized (TMM) expression
# Output: Dictionary with mean and median expression of isoforms
def AverageMedianExpressionTMMDict(path_to_file_with_expression):
    # Dictionary to return
    dict_mean_median_tmm = {}
    # Open file
    f = open(path_to_file_with_expression, 'r')
    # Iterate on lines in file
    for line in f:
        # If line does not start with hash
        if not line.startswith('#'):
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
            del expression_list, line_sep
            # Add to dictionary
            dict_mean_median_tmm[transcript_id] = {}
            dict_mean_median_tmm[transcript_id]['average'] = mean_tmm
            dict_mean_median_tmm[transcript_id]['median'] = median_tmm
            del transcript_id, mean_tmm, median_tmm
    # Close file
    f.close()
    # Delete variables not needed anymore
    del line, f
    # Return dictionary
    return dict_mean_median_tmm

# Generating two lists - one containing FASTA sequences of all transcripts containing novel/non-coding exons and second containing table with annotation of these transcripts
# Input:
#   - path_to_annotation_table_with_coding_status: Path to annotation table for each splicing pattern with coding status
#   - dict_mappings: Dictionary with all mappings
#   - dict_exprs_tmm_norm_mean_median: Dictionary with mean and median of all isoforms
# Output:
#   - List with lines of table with transcript information table
#   - List with lines for FASTA file
def LinesInfoTableFastaFileLists(path_to_annotation_table_with_coding_status, dict_mappings, dict_exprs_tmm_norm_mean_median):
    # List with lines of info table
    list_info_table_lines = []
    # List with FASTA lines for transcript containing novel exons
    list_fasta_lines_novel_exons = []
    # List with FASTA lines for transcript containing non-coding exons
    list_fasta_lines_non_coding_exons = []
    # Open file handle
    f = open(path_to_annotation_table_with_coding_status, 'r')
    # Iterate over lines in input file
    tq = tqdm.tqdm(total=os.path.getsize(path_to_annotation_table_with_coding_status), unit='B', unit_scale=True)
    for line in f:
        tq.update(len(line))
        # Omit header
        if not line.startswith('#'):
            # Split line by tab character
            line_sep = line.strip().split('\t')
            # If a transcript contains novel/non_coding exon
            if any([line_sep[5]=='novel_exon', line_sep[5]=='non_coding_exon']):
                # Isoform splicing pattern
                splicing_pattern = line_sep[0]
                # Transcript ID
                transcript_id = line_sep[2]
                # Split the splicing pattern to list of exons numeric IDs
                splicing_pattern_list = splicing_pattern.split('_')
                # Coding status
                coding_status = line_sep[5]
                # FASTA sequence list
                fasta_sequence_list = []
                # Novel/non-coding exons list
                novel_non_coding_exons_list = []
                # Iterate over exon IDs in splicing pattern list to obtain FASTA sequence
                for exon_num_id in splicing_pattern_list:
                    # Append FASTA sequence of each exon to list
                    fasta_sequence_list.append(dict_mappings[exon_num_id]['sequence'])
                    # If exon is novel/non-coding, append its' name to list of novel/non-coding exon names
                    if dict_mappings[exon_num_id]['coding_status'] != 'Yes':
                        novel_non_coding_exons_list.append(dict_mappings[exon_num_id]['exon_name'])
                # Delete unnecessary variables
                del exon_num_id, splicing_pattern_list
                # Join list with exon names by semicolon
                novel_non_coding_exons_str = ';'.join(novel_non_coding_exons_list)
                del novel_non_coding_exons_list
                # Join sequence of transcript by no character
                fasta_sequence_str = ''.join(fasta_sequence_list)
                del fasta_sequence_list
                # Length of transcript
                transcript_length = len(fasta_sequence_str)
                # Transcript length module 3
                transcript_len_mod_three = str(transcript_length % 3)
                transcript_length = str(transcript_length)
                # Average expression
                expression_average = dict_exprs_tmm_norm_mean_median[transcript_id]['average']
                # Median expression
                expression_median = dict_exprs_tmm_norm_mean_median[transcript_id]['median']
                # Final FASTA entry
                fasta_entry = '>' + transcript_id + '\n' + fasta_sequence_str + '\n'
                del fasta_sequence_str
                # If there is novel_exon, add fasta entry to the list for novel_exon
                if 'novel_exon' in novel_non_coding_exons_str:
                    list_fasta_lines_novel_exons.append(fasta_entry)
                if 'non_coding' in novel_non_coding_exons_str:
                    list_fasta_lines_non_coding_exons.append(fasta_entry)
                # Sanity check
                if not any(['novel_exon' in novel_non_coding_exons_str, 'non_coding' in novel_non_coding_exons_str]):
                    print('YOU DID SOMETHING WRONG WITH NON-CODING EXONS!')
                    print(line_sep)
                    sys.exit(1)
                del fasta_entry
                # Create a list for info table
                # Columns in the table are:
                # Splice_pattern, transcript_id, coding_status, non_coding_exons_names, transcript_length, transcript_length_mod_3, average_expression, median_expression
                info_table_line = [splicing_pattern, transcript_id, coding_status, novel_non_coding_exons_str, transcript_length, transcript_len_mod_three, expression_average, expression_median]
                # Join list by tab characters
                info_table_line = '\t'.join(info_table_line)
                # Append to the list of info lines
                list_info_table_lines.append(info_table_line)
                # Delete unnecessary variables
                del splicing_pattern, transcript_id, coding_status, novel_non_coding_exons_str, transcript_length, transcript_len_mod_three, expression_average, expression_median, info_table_line
            # Delete unnecessary variables
            del line_sep
    # Close file
    f.close()
    # Remove unnecessary variables
    del f, line, tq
    # Return lists
    return list_info_table_lines, list_fasta_lines_novel_exons, list_fasta_lines_non_coding_exons

# This function generates an output file
# Input:
#   - input_list: List with all the strings to write in the file
#   - header_bool: True if write a header, False if not
#   - header_string: String to write as header
#   - path_to_output_file: Path to output file
# Output: Output annotation table of BED12 file
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
args = parser.parse_args(['/tgac/workarea/group-eg/project_Capture/data/CACNA1C/July_17/CACNA1C_metagene_novel_exons/Results/novel_exons_remap_genome_exact_positions/novel_exons_check_positions.bed', '/tgac/workarea/group-eg/project_Capture/data/CACNA1C/July_17/CACNA1C_metagene_novel_exons/Meta_gene_files/meta_gene_genomic_exon_coordinates_UTR.txt', '/tgac/workarea/group-eg/project_Capture/data/hg38_ERCC_genome.fa', 'chr12', '/tgac/workarea/group-eg/project_Capture/data/CACNA1C/July_17/CACNA1C_metagene_novel_exons/Results/tmm_normalized_filtered_downweighted_read_counts/table.tmm_norm.sum_100.min_reads_per_sample_1.min_samples_2.txt', '/tgac/workarea/group-eg/project_Capture/data/CACNA1C/July_17/CACNA1C_metagene_novel_exons/Results/coding_length_extra_annotation/annotation_table.sum_100.min_reads_per_sample_1.min_samples_2.coding_length.txt'])
"""

if __name__ == '__main__':
    # Import the input files using argparse module
    parser = argparse.ArgumentParser(description='This script is used to generate FASTA files and info table for transcripts having novel/non-coding exons.')
    parser.add_argument("in_narrow_novel_exons", type=str, help="Path to the file with narrowed start/stop positions of novel exons (based on coverage and canonical splice sites).")
    parser.add_argument("in_mappings", type=str, help="Path to the file with exonic positions in genome and meta-gene.")
    parser.add_argument("in_genomic_fasta", type=str, help="Path to FASTA file with genomic sequence.")
    parser.add_argument("chrom_string", type=str, help="Name of the chromosome where a gene is located.")
    parser.add_argument("in_expression_tmm", type=str, help="Path to the file with TMM-normalized expression values for all isoforms.")
    parser.add_argument("in_annotations", type=str, help="Path to the file with all splice patterns and their annotations (including coding status).")
    parser.add_argument("out_transcript_info", type=str, help="Path to table with info of non-coding transcripts.")
    parser.add_argument("out_fasta_novel_exons", type=str, help="Path to output FASTA file with sequences for transcripts containing novel exons.")
    parser.add_argument("out_fasta_non_coding_exons", type=str, help="Path to output FASTA file with sequences for transcripts containing non-coding exons.")
    args = parser.parse_args()
    # Dictionary with 'narrowed' positions of novel exons
    dict_narrowed_novel_exons = NovelExonStartStopDict(args.in_narrow_novel_exons)
    # Creating dictionary with start, end, coding status, exon name and sequence
    mappings_dict = MappingsDictionary(args.in_mappings, args.in_genomic_fasta, args.chrom_string, dict_narrowed_novel_exons)
    del dict_narrowed_novel_exons
    gc.collect()
    # Creating dictionary with average and median expression of each isoforms
    mean_median_expression_dict = AverageMedianExpressionTMMDict(args.in_expression_tmm)
    # Creating lists with FASTA sequences and table with transcripts containing novel/non-coding exons
    info_table_lines_list, fasta_lines_novel_exons_list, fasta_lines_non_coding_exons_list = LinesInfoTableFastaFileLists(args.in_annotations, mappings_dict, mean_median_expression_dict)
    del mappings_dict, mean_median_expression_dict
    gc.collect()
    # Write table with information about isoforms
    WriteOutputFile(info_table_lines_list, True, '#Splicing_pattern\tTranscript_ID\tCoding_status\tNovel_non_coding_exon_names\tTranscript_length\tTranscript_length_mod_3\tAverage_expression\tMedian_expression\n', args.out_transcript_info)
    # Write FASTA file with sequences for transcripts containing novel exons
    WriteOutputFile(fasta_lines_novel_exons_list, False, None, args.out_fasta_novel_exons)
    # Write FASTA file with sequences for transcripts containing non-coding exons
    WriteOutputFile(fasta_lines_non_coding_exons_list, False, None, args.out_fasta_non_coding_exons)
    del info_table_lines_list, fasta_lines_novel_exons_list, fasta_lines_non_coding_exons_list
    gc.collect()
