#!/usr/bin/env python3

# -*- coding: utf-8 -*-

# Importing modules
import os
import sys
import pybedtools as pb
import tqdm
import gc
import argparse

# This function creates a dictionary (0-based) of all mappings
# Input:
#   - path_to_input_file: Path to the file with all mappings
#   - path_to_genome_fasta: Path to genomic FASTA file
#   - chrom: Chromosome (must match chromosome in FASTA file)
# Output: Dictionary of following structure:
# Dictionary:
#   - key: string with exon number (row in file with mappings)
#       - key: 'start'
#           - value: start position (0-based)
#       - key: 'stop'
#           - value: stop position
#       - key: 'coding_status'
#           - value: 'novel_exon', 'Yes', 'No'
#       - key: 'sequence'
#           - value: sequence of exon
def MappingsDictionary(path_to_input_file, path_to_genome_fasta, chrom, strand):
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
            # Start position (0-based)
            start_pos = str(int(line_sep[2])-1)
            # Stop position
            end_pos = line_sep[3]
            # Coding status
            if 'novel_exon' in line_sep[4]:
                coding_status = 'novel_exon'
            elif not 'Yes' in line_sep[6]:
                coding_status = 'No'
            else:
                coding_status = 'Yes'
            # If coding status is 'Yes', obtain a sequence
            if coding_status == 'Yes':
                # Obtain a space separated string with chromosome, start and end
                bed_string = [chrom, start_pos, end_pos, '.', '.', strand]
                bed_string = ' '.join(bed_string)
                # Convert coordinates to BedTool
                a = pb.BedTool(bed_string, from_string=True)
                del bed_string
                # Extract sequence for a particular region
                genomic_sequence = a.sequence(fi=path_to_genome_fasta, s=True)
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
            # If coding status is not 'Yes', put None in genomic_sequence_fasta
            else:
                genomic_sequence_fasta = None
            # Add start position to dictionary
            dict_all_mappings[str(exon_number)]['start'] = start_pos
            # Add stop position to dictionary
            dict_all_mappings[str(exon_number)]['stop'] = end_pos
            # Add coding status to dictionary
            dict_all_mappings[str(exon_number)]['coding_status'] = coding_status
            # Add sequence to the dictionary
            dict_all_mappings[str(exon_number)]['sequence'] = genomic_sequence_fasta
            # Delete unnecessary variables
            del start_pos, end_pos, coding_status, genomic_sequence_fasta
    # Close file handle
    f.close()
    # Delete unnecessary variables
    del f, exon_number, tq, line
    # Return dictionary
    return dict_all_mappings

# This function is to assess coding status of a particular splicing pattern
# If isoform is coding, it returns string with sequence
# Input:
#   - list_splicing_pattern: List with numeric IDs of exons
#   - dict_exonic_num_ids: Dictionary with mappings
# Output: Coding status of a particular isoform
def CodingStatusAssessFASTASequence(list_splicing_pattern, dict_exonic_num_ids):
    # Final coding status is None
    final_coding_status = None
    # Sequence string to output
    sequence_string_to_output = None
    # Create a sequence list
    sequence_list = []
    # Length of coding isoform
    coding_isoform_length_string = 'N/A'
    # Iterate over list of exonic num IDs
    for exonic_num_id in list_splicing_pattern:
        # Obtain a coding status of exon
        coding_status_in_dict = dict_exonic_num_ids[exonic_num_id]['coding_status']
        # If coding status is 'Yes', append a sequence to a sequence list
        if coding_status_in_dict == 'Yes':
            sequence_list.append(dict_exonic_num_ids[exonic_num_id]['sequence'])
            del coding_status_in_dict
        # If coding status is No, final_coding_status is 'non_coding_exon'
        elif coding_status_in_dict == 'No':
            final_coding_status = 'non_coding_exon'
            del coding_status_in_dict
            break
        # Else exit script
        else:
            print('THERE IS WEIRD CODING STATUS IN THE DICTIONARY - CHECK IT!!')
            print(exonic_num_id)
            sys.exit(1)
    del exonic_num_id
    # If coding status is None
    if not final_coding_status:
        # Join the sequence list
        sequence_string = ''.join(sequence_list)
        # Sequence string to output if isoform is coding
        sequence_string_to_output_temp = sequence_string
        # Delete last three characters
        sequence_string = sequence_string[:-3]
        # If length of string module 3 is not 0, the isoform is frame_disruption
        if len(sequence_string) % 3 != 0:
            final_coding_status = 'frame_disruption'
        else:
            # Divide sequence into list with three-character strings
            n=3
            sequence_codons = [sequence_string[i:i+n] for i in range(0, len(sequence_string), n)]
            del n
            # If there is a stop codon, final coding status is 'premature_stop_codon'
            if any(['TAA' in sequence_codons, 'TAG' in sequence_codons, 'TGA' in sequence_codons]):
                final_coding_status = 'premature_stop_codon'
            # Else the sequence is coding
            else:
                final_coding_status = 'coding'
                sequence_string_to_output = sequence_string_to_output_temp
                coding_isoform_length_string = str(len(sequence_string_to_output))
            del sequence_codons
        del sequence_string, sequence_string_to_output_temp
    del sequence_list
    # If final_coding_status is still None, return error
    if not final_coding_status:
        print('ERROR IN CODING STATUS - CHECK IT!')
        print(list_splicing_pattern)
        sys.exit(1)
    else:
        # Return coding status
        return final_coding_status, sequence_string_to_output, coding_isoform_length_string

# This function is to obtain RGB-coded colours for different coding status
# Input: Coding status
# Output: RGB string with colour
def ItemRgbMatcher(string_with_coding_status):
    if string_with_coding_status == 'novel_exon':
        item_rgb_string = '255,0,0'
    elif string_with_coding_status == 'non_coding_exon':
        item_rgb_string = '255,165,0'
    elif string_with_coding_status == 'frame_disruption':
        item_rgb_string = '0,0,255'
    elif string_with_coding_status == 'premature_stop_codon':
        item_rgb_string = '0,255,0'
    elif string_with_coding_status == 'coding':
        item_rgb_string = '0,0,0'
    else:
        print('WRONG CODING STATUS IN RGB GENERATION - CHECK IT')
        print(string_with_coding_status)
        sys.exit(1)
    # Return item RGB
    return item_rgb_string

# This function is to generate line in BED12 format
# Input:
#   - chromosome: string with chromosome
#   - spl_patt_list: List with splicing pattern
#   - cod_stat_str: Coding status of isoform
#   - trans_id_str: ID of transcript
#   - mapp_all_dict: Mappings dictionary
# Output: BED12 line
def BED12LineCreator(chromosome, spl_patt_list, cod_stat_str, trans_id_str, mapp_all_dict, strand):
    # List of tuples with start and stop for each exon
    start_stop_tuples_list = []
    # Iterate over splice patterns list
    for exon_number in spl_patt_list:
        # Start position
        start_position = int(mapp_all_dict[exon_number]['start'])
        # Stop position
        stop_position = int(mapp_all_dict[exon_number]['stop'])
        # Create a tuple with start and stop
        tuple_to_append = (start_position, stop_position)
        # Append tuple to list
        start_stop_tuples_list.append(tuple_to_append)
        del start_position, stop_position
    del exon_number
    # Sort list of tuples
    start_stop_tuples_list.sort()
    # Obtain start and stop positions in transcript
    start_pos_transcript = min([x[0] for x in start_stop_tuples_list])
    stop_pos_transcript = max([x[1] for x in start_stop_tuples_list])
    # Create a list with start positions relative to start and block sizes
    block_start_list = []
    block_size_list = []
    # Iterate over list of tuples to create a block size and block starts
    for blocks in start_stop_tuples_list:
        # Block start
        iteration_block_start = blocks[0] - start_pos_transcript
        # Block size
        iteration_block_size = blocks[1] - blocks[0]
        # Append to lists
        block_start_list.append(str(iteration_block_start))
        block_size_list.append(str(iteration_block_size))
        # Delete unnecessary variables
        del iteration_block_size, iteration_block_start
    # Delete unnecessary variables
    del blocks, start_stop_tuples_list
    # Block count
    block_count = len(block_start_list)
    # Convert blocks lists to comma-delimited strings
    block_start_str = ','.join(block_start_list)
    block_size_str = ','.join(block_size_list)
    del block_start_list, block_size_list
    # Obtain itemRgb
    itemRgb_str = ItemRgbMatcher(cod_stat_str)
    # Final list with all columns
    final_list_bed12_line = [chromosome, str(start_pos_transcript), str(stop_pos_transcript), trans_id_str, '0', strand, str(start_pos_transcript), str(stop_pos_transcript), itemRgb_str, str(block_count), block_size_str, block_start_str]
    # Delete unnecessary variables
    del start_pos_transcript, stop_pos_transcript, itemRgb_str, block_count, block_size_str, block_start_str
    # Join list by tab characters
    line_to_return = '\t'.join(final_list_bed12_line)
    return line_to_return

# This function creates a line to write in the annotation table with coding_status column
# Input:
#   - line_string: Line string of annotation table
#   - coding_status_string: Coding status string to add at the end of the line
# Output: Line with added coding status
def AnnotationTableLineCreator(line_string, coding_status_string):
    # Strip the line
    to_return = line_string.strip()
    # Add the coding status after tab character
    to_return = to_return + '\t' + coding_status_string
    # Return table string
    return to_return

# This function generates two output lists - one having lines of annotated table with coding status of isoform and second list with BED12 lines
# Input:
#   - chrom_str: String with chromosome
#   - path_to_annotation_table: Path to annotation table for each splice isoform
#   - dict_mappings: Dictionary with all mappings
# Output:
#   - List with lines of table with coding status
#   - List with lines for bed12 file
def LinesBed12CodingStatusLists(chrom_str, path_to_annotation_table, dict_mappings, strand):
    # List with lines of annotation table
    list_annotation_table_coding_status = []
    # List with lines of BED12 file
    list_bed12_file = []
    # List with FASTA lines for coding sequences
    list_fasta_lines = []
    # Open file handle
    f = open(path_to_annotation_table, 'r')
    # Iterate over lines in input file
    tq = tqdm.tqdm(total=os.path.getsize(path_to_annotation_table), unit='B', unit_scale=True)
    for line in f:
        tq.update(len(line))
        # Omit header
        if not line.startswith('#'):
            # Split line by tab character
            line_sep = line.strip().split('\t')
            # Isoform splicing pattern
            splicing_pattern = line_sep[0]
            # Transcript ID
            transcript_id = line_sep[2]
            # Split the splicing pattern to list of exons numeric IDs
            splicing_pattern_list = splicing_pattern.split('_')
            del splicing_pattern
            # Isoform sequence is None
            isoform_sequence = None
            # Coding length is N/A
            coding_length = 'N/A'
            # If the isoform has novel exon, change coding_status to 'novel_exon'
            if line_sep[3] == 'Yes':
                coding_status = 'novel_exon'
            else:
                # Obtain coding status
                coding_status, isoform_sequence, coding_length = CodingStatusAssessFASTASequence(splicing_pattern_list, dict_mappings)
            del line_sep
            # If isoform sequence is not None, append transcript id as header and sequence in another line
            if isoform_sequence:
                string_to_append = '>' + transcript_id + '\n' + isoform_sequence + '\n'
                list_fasta_lines.append(string_to_append)
            del isoform_sequence
            # Obtain BED12 line
            bed12_final_line = BED12LineCreator(chrom_str, splicing_pattern_list, coding_status, transcript_id, dict_mappings, strand)
            # Create line of annotation
            annotation_final_line = AnnotationTableLineCreator(line, coding_status)
            # Adding length of coding sequence to annotation table
            annotation_final_line = annotation_final_line + '\t' + coding_length
            del coding_length
            # Append all lines to lists
            list_annotation_table_coding_status.append(annotation_final_line)
            list_bed12_file.append(bed12_final_line)
            # Delete unnecessary variables
            del coding_status, transcript_id, splicing_pattern_list, bed12_final_line, annotation_final_line
    # Close the file handle
    f.close()
    # Delete unnecessary variables
    del f, tq, line
    # Return annotation list and BED12 list
    return list_annotation_table_coding_status, list_bed12_file, list_fasta_lines

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
args = parser.parse_args(['/tgac/workarea/group-eg/project_Capture/data/CACNA1C/July_17/CACNA1C_metagene_novel_exons/Meta_gene_files/meta_gene_genomic_exon_coordinates_UTR.txt', '/tgac/workarea/group-eg/project_Capture/data/hg38_ERCC_genome.fa', 'chr12', '/tgac/workarea/group-eg/project_Capture/data/CACNA1C/July_17/CACNA1C_metagene_novel_exons/Results/annotation_table_filtered_downweighted_read_counts/annotation_table.sum_100.min_reads_per_sample_1.min_samples_2.txt'])
"""

if __name__ == '__main__':
    # Import the input files using argparse module
    parser = argparse.ArgumentParser(description='This script is used to generate a BED12 file, FASTA file with coding sequences and annotation of isoform table with coding status of isoforms.')
    parser.add_argument("in_mappings", type=str, help="Path to the file with exonic positions in genome and meta-gene.")
    parser.add_argument("in_genomic_fasta", type=str, help="Path to FASTA file with genomic sequence.")
    parser.add_argument("chrom_string", type=str, help="Name of the chromosome where a gene is located.")
    parser.add_argument("strand_string", type=str, help="Strand of the gene.")
    parser.add_argument("in_annotations", type=str, help="Path to the file with all splice patterns and their annotations.")
    parser.add_argument("out_annotation_table", type=str, help="Path to table with coding status and length of coding sequences.")
    parser.add_argument("out_bed12", type=str, help="Path to output BED12 file.")
    parser.add_argument("out_fasta", type=str, help="Path to output FASTA file with sequence of coding transcripts.")
    args = parser.parse_args()
    # Creating dictionary with start, end, coding status and sequence
    mappings_dict = MappingsDictionary(args.in_mappings, args.in_genomic_fasta, args.chrom_string, args.strand_string)
    # BED-like file with genomic position and all exon/transcript IDs
    annotation_list, bed12_list, fasta_sequences_list = LinesBed12CodingStatusLists( args.chrom_string, args.in_annotations, mappings_dict, args.strand_string)
    del mappings_dict
    gc.collect()
    # Write annotation file
    WriteOutputFile(annotation_list, True, '#Splicing_pattern\tAnnotated_novel\tTranscript_ID\tNovel_exon_status\tNovel_exon_IDs\tCoding_status\tCoding_isoform_length\n', args.out_annotation_table)
    # Write BED12 file
    WriteOutputFile(bed12_list, True, 'track name="all_isoforms" description="all_isoforms" visibility=2 itemRgb="On"\n', args.out_bed12)
    # Write FASTA file with coding sequences
    WriteOutputFile(fasta_sequences_list, False, None, args.out_fasta)
    del annotation_list, bed12_list, fasta_sequences_list
    gc.collect()
