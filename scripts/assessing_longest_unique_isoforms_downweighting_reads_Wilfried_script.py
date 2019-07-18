#!/usr/bin/env python3

# -*- coding: utf-8 -*-

# Importing modules
import os
import sys
import tqdm
import gc
import argparse
import pathlib
import copy

# This function iterates over the file with matching barcodes and sequencing runs, and creates a dictionary with matched barcodes and sample names
# Input: Path to the file with barcode-to-sample mappings
#   - path_to_file
#   - barcode: String with barcode (e.g. 'barcode01')
# Output: Dict with sequencing run dates as keys, barcodes as sub-keys and sample names as strings in sub-dictionaries
def IndividualBrainRegionMatcherDict(path_to_file):
    # Create output dictionary
    output_dict = {}
    # Create output list
    output_list = []
    # Open file handle
    f = open(path_to_file, 'r')
    # Iterate on lines of file
    for line in f:
        # Omit header line
        if not line.startswith('#'):
            # Split line by tab character
            line_sep = line.strip().split('\t')
            # First column is sequencing run
            sequencing_run = line_sep[0]
            # Second column is barcode
            barcode = line_sep[1]
            # Third column is sample name
            sample_name = line_sep[2]
            # Append sample name to output list
            output_list.append(sample_name)
            # If key is not in dictionary, create it
            if sequencing_run not in output_dict:
                output_dict[sequencing_run] = {}
            # If barcode is already is sub-dictionary, raise an error
            if barcode in output_dict[sequencing_run]:
                print("BARCODE ALREADY IN DICTIONARY!")
                print(barcode)
                sys.exit(1)
            # Add key to sub-dictionary
            output_dict[sequencing_run][barcode] = sample_name
            # Remove unnecessary variables
            del line_sep, sequencing_run, barcode, sample_name
    # Close file handle
    f.close()
    # Delete unnecessary variables
    del line, f
    # Return dictionary
    return output_dict, output_list

# This function obtains a list of files with splicing patterns
# Input:
#   - dir_path: Path to the directory with files
#   - file_distinct: string to distinguish files of interest from other files in possible subdirectories
# Output: List of files and sample names (in the same order, i.e. first element in a file name list corresponds to the first element in the sample names list)
def FileList(dir_path, file_distinct):
    path_to_return = [] # Create empty list
    # Append all the files found in a specific directory (and its' subdirectories) to the list
    for dir_path, subdirs, files in os.walk(dir_path):
        for name in files:
            path_to_return.append(str(pathlib.PurePath(dir_path, name)))
        # This break ensures that the list of files consists only the top directory
        break
    path_to_return = [x for x in path_to_return if file_distinct in x] # Filter out specific filenames and sort them
    path_to_return.sort()
    # Create a list of sample names by parsing the filenames
    return path_to_return # Return file paths

# This function reads each file in file list and extract unique splicing patterns for each read
# Input: List with file paths
# Output: List containing all unique splicing patterns
def AllUniqueIsoformsList(file_list):
    # Creating list to return
    splice_pattern_list_to_return = []
    # Iterating over file list
    for file_path in file_list:
        # Open file handle
        f = open(file_path, 'r')
        # Iterate over contents of file
        for line in f:
            # If line does not start with hash
            if not line.startswith('#'):
                # Strip and split line by tab-character
                line_sep = line.strip().split('\t')
                # Obtaining splicing pattern
                splicing_pattern = line_sep[1]
                # Adding splicing pattern to final list
                splice_pattern_list_to_return.append(splicing_pattern)
                # Delete variables not needed anymore
                del line_sep, splicing_pattern
        # Close the file handle
        f.close()
        # Delete variables not needed anymore
        del f, line
    # Delete variables not needed anymore
    del file_path
    # Remove duplicated patterns from list
    seen = set()
    splice_pattern_list_to_return = [x for x in splice_pattern_list_to_return if x not in seen and not seen.add(x)]
    del seen
    # Sorting list by length (descending)
    splice_pattern_list_to_return.sort(key=len, reverse=True)
    # Return list with all splice patterns
    return splice_pattern_list_to_return

# This function generates a list with the longest unique splicing patterns
# Input: List with all splicing patterns
# Output: Dictionary with all splicing pattern containing a list of all longest isoforms for this specific splicing pattern
# Workflow - for each splicing pattern:
#   - Choose remove the isoform from list of all isoforms
#   - Match all isoforms contained within splicing pattern
#   - Match splicing pattern contained within isoforms
#   - If there is no splicing patterns contained within any of isoforms, the splicing pattern is the longest isoform
#       - In such case, add to the dictionary - all isoforms within splicing pattern as key and splicing pattern as the longest isoform in list within each key (appending to the list)
def LongestSplicePatternClumpingList(list_with_all_splicing_patterns):
    # Creating dictionary with the longest unique isoforms
    longest_unique_splicing_patterns_dict_to_return = {}
    # Iterating over the list of splicing patterns
    for i in tqdm.trange(len(list_with_all_splicing_patterns)):
        splicing_pattern = list_with_all_splicing_patterns[i]
        # Obtain list with all splicing patterns with this particular splicing_pattern removed from list
        list_without_particular_splicing_pattern = [x for x in list_with_all_splicing_patterns if x!=splicing_pattern]
        # Obtain all patterns being within this particular splicing pattern
        matching_x_within_splicing_patterns = [x for x in list_without_particular_splicing_pattern if splicing_pattern.startswith(str(x + "_")) or splicing_pattern.endswith(str("_" + x)) or str("_" + x + "_") in splicing_pattern]
        # Obtain all patterns where splicing_pattern is within these patterns
        matching_splicing_patterns_within_x = [x for x in list_without_particular_splicing_pattern if x.startswith(str(splicing_pattern + "_")) or x.endswith(str("_" + splicing_pattern)) or str("_" + splicing_pattern + "_") in x]
        # If length of matching_splicing_patterns_within_x is 0, then the splicing_pattern is the longest one
        if len(matching_splicing_patterns_within_x) == 0:
            # Sanity check: the key should not exist in the dictionary
            # If the key exists, quit immediately
            if splicing_pattern in longest_unique_splicing_patterns_dict_to_return.keys():
                print("SPLICING PATTERN ALREADY IN THE DICTIONARY")
                print(splicing_pattern)
                sys.exit(1)
            # Add the longest isoform to the dictionary
            longest_unique_splicing_patterns_dict_to_return[splicing_pattern] = [splicing_pattern]
            # Iterating over a list of patterns within splice_pattern
            for splicing_pattern_item in matching_x_within_splicing_patterns:
                # Try appending the splicing pattern into a list in the dictionary
                try:
                    longest_unique_splicing_patterns_dict_to_return[splicing_pattern_item].append(splicing_pattern)
                # When the key is not here, create it
                except KeyError:
                    longest_unique_splicing_patterns_dict_to_return[splicing_pattern_item] = [splicing_pattern]
        # Sanity checks
        else:
            # Splicing pattern must be a key in dictionary
            if splicing_pattern not in longest_unique_splicing_patterns_dict_to_return.keys():
                print("SANITY CHECK FAILED!")
                print("SPLICING PATTERN NOT IN DICTIONARY")
                print(item_in_splicing_pattern)
                print(splicing_pattern)
                sys.exit(1)
            for item_in_splicing_pattern in matching_splicing_patterns_within_x:
                # All of splicing_patterns within x must be already keys in dictionary
                if item_in_splicing_pattern not in longest_unique_splicing_patterns_dict_to_return.keys():
                    print("SANITY CHECK FAILED!")
                    print("LONGEST SPLICING PATTERN NOT IN DICTIONARY!")
                    print(item_in_splicing_pattern)
                    print(splicing_pattern)
                    sys.exit(1)
        # Deleting unnecessary variables
        del list_without_particular_splicing_pattern, matching_x_within_splicing_patterns, matching_splicing_patterns_within_x, splicing_pattern
    # Deleting unnecessary variables
    del i
    # Return dictionary
    return longest_unique_splicing_patterns_dict_to_return


# This function creates empty dictionary with longest isoforms as keys, sample names as subkeys and counts equal to 0
# Input:
#   - dictionary_longest_isoforms: Dictionary containing all isoforms as keys and list of the longest isoforms matching this key
#   - list_of_samples: List of all samples (sorted by sequencing date and barcode)
# Output: Dictionary with longest isoforms as keys, sample names as subkeys in sub-dictionary and counts equal to 0
def CreatingEmptyDictionaryZeroCounts(dictionary_longest_isoforms, list_of_samples):
    # Creating dictionary to return
    to_return_dict = {}
    # Iterating over keys in dictionary
    for splicing_pattern in dictionary_longest_isoforms.keys():
        # Iterating over all longest isoform in list for a particular splicing pattern in dictionary
        for longest_isoform in dictionary_longest_isoforms[splicing_pattern]:
            # If key does not exist
            if longest_isoform not in to_return_dict.keys():
                # Create sub-dictionary
                to_return_dict[longest_isoform] = {}
            # Iterate through all samples in sample list
            for sample_name in list_of_samples:
                # If key does not exist, create it with 0
                if sample_name not in to_return_dict[longest_isoform].keys():
                    to_return_dict[longest_isoform][sample_name] = 0
                # Sanity check - if count is not equal to zero or key does not exist, exit
                if to_return_dict[longest_isoform][sample_name] != 0:
                    print("EMPTY DICTIONARY HAS SOMETHING DIFFERENT THAT ZERO INTEGER!")
                    sys.exit(1)
            # Delete unnecessary variables
            del sample_name
        # Delete unnecesarry variables
        del longest_isoform
    # Delete variables not needed anymore
    del splicing_pattern
    # Return empty dictionary
    return to_return_dict


# This function iterates on all files and creates a dictionary with all downweighted reads
# Input:
#   - empty_dictionary_zero_counts: Dictionary with zero counts, with longest isoforms as keys, sample names as subkeys and 0 counts
#   - input_file_list: List of all filenames to iterate on
#   - dictionary_longest_isoforms: Dictionary with all splicing patterns as keys and list of the longest isoform for each key
# Output: Dictionary containing downweighted read counts (structure the same as for zero-count dictionary)
def DownweightedReadCounts(empty_dictionary_zero_counts, input_file_list, dictionary_longest_isoforms, dict_mappings):
    # Output dictionary - copy of dictionary with zero counts
    dictionary_downweighted_to_return = copy.deepcopy(empty_dictionary_zero_counts)
    # Iterating over file list
    for i in tqdm.trange(len(input_file_list)):
        file_name = input_file_list[i]
        # Opening file handle
        f = open(file_name, 'r')
        # Sequencing run date
        sequencing_run_date = file_name.split('/')[-1].split('.')[0]
        # Barcode
        barcode = file_name.split('/')[-1].split('.')[1]
        # String with sample name
        sample_name_str = dict_mappings[sequencing_run_date][barcode]
        # Iterating over lines in file
        for line in f:
            # If line does not start with hash
            if not line.startswith('#'):
                # Strip and split line by tab-character
                line_sep = line.strip().split('\t')
                # Obtaining splicing pattern
                splicing_pattern = line_sep[1]
                # Obtaining list with all longest isoforms from dictionary
                longest_isoforms_splicing_pattern = dictionary_longest_isoforms[splicing_pattern]
                # Calculating the downweighted read count (6 decimals)
                downweighted_read_count = round(1/len(longest_isoforms_splicing_pattern), 6)
                # Iterating over list of longest isoforms
                for isoform in longest_isoforms_splicing_pattern:
                    # Increase count for a specific isoform by downweighted_read_count
                    dictionary_downweighted_to_return[isoform][sample_name_str] += downweighted_read_count
                # Delete unnecessary variables
                del line_sep, splicing_pattern, longest_isoforms_splicing_pattern, downweighted_read_count, isoform
        # Close file
        f.close()
        # Delete unnecessary variables
        del file_name, f, sequencing_run_date, barcode, sample_name_str, line
    # Delete unnecessary variables
    del i
    # Return output dictionary
    return dictionary_downweighted_to_return

# This function iterates over a file with all mappings and creates a dictionary with splice patterns and ENSEMBL Transcript IDs (ENST)
# Input: Path to file with all mappings
# Output: Dictionary with splice patterns as keys and ENSEMBL Transcript IDs as values
def SplicePatternsKnownTranscriptsDicts(path_to_file_with_mappings):
    # Creating dictionary with splice patterns of exisitng transcripts
    dict_splice_patterns_existing_transcripts = {}
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
            # If this line contains 'novel_exon' string then continue
            if 'novel_exon' in transcript_novel_exons:
                del line_sep, transcript_novel_exons
                continue
            # Else there are existing transcripts, so append line_num to list for each transcript ID (semicolon-delimited)
            else:
                # Split the list with ENSEMBL transcript IDs by semicolon
                existing_transcripts = transcript_novel_exons.split(';')
                # For every ENSEMBL transcript ID, add the line_num to the list in the dictionary of splice patterns
                for transcript_id in existing_transcripts:
                    # With .setdefault method, append the line_num to the list, if list is not there because transcript_id is not a key in the dictionary, create this key with empty list and subsequently append
                    dict_splice_patterns_existing_transcripts.setdefault(transcript_id, []).append(str(line_num))
                # Delete unnecessary variables
                del transcript_id, existing_transcripts
            # Delete unnecessary variables
            del line_sep, transcript_novel_exons
    # Close the file handle
    f.close()
    # Delete unnecessary variables
    del f, line_num, line
    # For every list in sub-dictionary convert a list of exon number to underscore-delimited string
    for ensembl_transcript_id in dict_splice_patterns_existing_transcripts:
        dict_splice_patterns_existing_transcripts[ensembl_transcript_id] = "_".join(dict_splice_patterns_existing_transcripts[ensembl_transcript_id])
    # Invert keys and items in dictionary of splice patterns, so every key is item and item is key
    dict_splice_patterns_existing_transcripts = {v: k for k,v in dict_splice_patterns_existing_transcripts.items()}
    # Return dictionary witb exisiting transcripts (in this order)
    return dict_splice_patterns_existing_transcripts

# This function writes output table to the specified output file
# Input:
#   - gene_name: Gene name
#   - input_dictionary: Input dictionary with downweighted read counts
#   - known_transcripts_dict: Dictionary with splicing patterns of known transcripts
#   - list_of_names_individuals_brain_regions: List of all sample names (sorted)
#   - output_file_path: Path to the output file
# Output: Tab-delimited file with downweighted read counts for each sample (sorted by read and by barcode)
def WriteOutputFiles(gene_name, input_dictionary, known_transcripts_dict, list_of_names_individuals_brain_regions, output_file_path):
    # Opening output file handle
    g = open(output_file_path, 'w')
    # Adding header
    header_to_write = '#Transcript_splicing_pattern\tTranscript_ID\t' + '\t'.join(list_of_names_individuals_brain_regions) + '\n'
    # Writing header
    g.write(header_to_write)
    # List with all strings
    to_write_list = []
    # Removing header string
    del header_to_write
    # Number for novel transcript
    novel_transcript_num = 0
    # Iterating over list of exons (sorted as integers)
    exon_id_list = list(input_dictionary.keys())
    exon_id_list.sort(key=len, reverse=True)
    for exon_num_id in exon_id_list:
        # Obtaining tab-delimited string with read counts strings
        string_with_read_counts = [str(input_dictionary[exon_num_id][x]) for x in list_of_names_individuals_brain_regions]
        string_with_read_counts = '\t'.join(string_with_read_counts)
        # If a splice pattern is a key in the dictionary with known splice patterns, obtain transcript ID and transcript status
        if exon_num_id in known_transcripts_dict.keys():
            transcript_id = known_transcripts_dict[exon_num_id]
        # If a splice pattern is not in dictionary of known splice patterns, transcript ID is 'CACNA1C_novel-' + number and transcript status is 'Novel'
        else:
            # Increment novel transcript number for one
            novel_transcript_num += 1
            transcript_id = gene_name + '_novel-' + str(novel_transcript_num)
        # Final string with exon number
        final_string = str(exon_num_id) + '\t' + transcript_id + '\t' + string_with_read_counts
        # Append the final string to final list
        to_write_list.append(final_string)
        # Delete variables that are not needed anymore
        del string_with_read_counts, final_string
    # Delete variables that are not needed anymore
    del exon_num_id, exon_id_list, novel_transcript_num
    # Joining the list by new-line character
    to_write_string = '\n'.join(to_write_list)
    del to_write_list
    to_write_string = to_write_string + '\n'
    # Write all the lines
    g.write(to_write_string)
    # Close the file and delete all variables
    g.close()
    del g, to_write_string



# Testing arguments
"""
args = parser.parse_args(['/tgac/workarea/group-eg/project_Capture/data/CACNA1C/July_17/CACNA1C_metagene_novel_exons/Results/reads_splicing_patterns_removed_exons', '.included.txt', 'out_table_counts_all_isoforms.txt'])
"""


if __name__ == '__main__':
    # Import the input files using argparse module
    parser = argparse.ArgumentParser(description='This script assess the longest isoform and down-weight reads mapping to multiple longest isoforms. It reports downweighted counts for each sample.')
    parser.add_argument('in_barcode_sample_name', type=str, help='Path to the file with barcode-to-sample-name mappings.')
    parser.add_argument('in_gene_name', type=str, help='Gene name.')
    parser.add_argument('in_splice_patterns_dir', type=str, help='Directory with files containing reads and their splicing patterns')
    parser.add_argument('in_distinct_string', type=str, help='String allowing to distinguish files of interest from other files.')
    parser.add_argument('in_mappings', type=str, help='Path to the file with 1-based positions of each exon (with separated coding/non-coding parts of the exons to another lines in file) in genome, metagene, transcript IDs, exon IDs and coding status of the sequence.')
    parser.add_argument('out_table', type=str, help='Path to the output file with downweighted read counts')
    args = parser.parse_args()
    # Create dictionary with barcode-to-name mappings and list with sample names
    barcode_sample_names_dict, sample_name_list_sorted = IndividualBrainRegionMatcherDict(args.in_barcode_sample_name)
    # List of all paths with files containing reads and their corresponding splicing patterns
    file_path_list = FileList(args.in_splice_patterns_dir, args.in_distinct_string)
    # List with all splicing patterns among all files
    all_patterns_list = AllUniqueIsoformsList(file_path_list)
    # Dictionary with all splicing patterns and list of splicing patterns they should be associated to
    splicing_pattern_dictionary = LongestSplicePatternClumpingList(all_patterns_list)
    del all_patterns_list
    gc.collect()
    # Creating empty dictionary with zero counts for each longest isoform
    dictionary_zero_counts_all_longest_isoforms = CreatingEmptyDictionaryZeroCounts(splicing_pattern_dictionary, sample_name_list_sorted)
    # Dictionary with the longest isoforms and sample names
    downweighted_counts_dict = DownweightedReadCounts(dictionary_zero_counts_all_longest_isoforms, file_path_list, splicing_pattern_dictionary, barcode_sample_names_dict)
    del dictionary_zero_counts_all_longest_isoforms, file_path_list, splicing_pattern_dictionary
    gc.collect()
    # Splice patterns of known isoforms
    dict_known_isoforms_splice_patterns = SplicePatternsKnownTranscriptsDicts(args.in_mappings)
    # Writing output table
    WriteOutputFiles(args.in_gene_name, downweighted_counts_dict, dict_known_isoforms_splice_patterns, sample_name_list_sorted, args.out_table)
    del downweighted_counts_dict, sample_name_list_sorted, dict_known_isoforms_splice_patterns
    gc.collect()
