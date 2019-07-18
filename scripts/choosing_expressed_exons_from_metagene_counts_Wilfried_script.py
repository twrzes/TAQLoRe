#!/usr/bin/env python3

# -*- coding: utf-8 -*-

# Importing modules
import os
import sys
import tqdm
import gc
import argparse
import pathlib

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


# This function generates the list of files with total read counts covering each exon for a particular sample
# Input:
#   - dir_path: Path to the directory with files with total read counts covering each exon for a particular sample
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

# This function reads exon counts and matches sample name (by sequencing run and barcode) to the brain region and individual
# Input:
#   - directory_with_files: Directory where files with exon counts are stored
#   - string_to_choose_file: String that distinguishes files of interest from other files
#   - dict_mappings: Dictionary with barcode-to-name mappings
# Output: Dictionary with exon numbers as keys and subdictionary - sample names (sequencing run, individual ID, brain region) as keys and read numbers (integers) as values
def ReadExonsCountsMatchedSampleNames(directory_with_files, string_to_choose_file, dict_mappings):
    # Obtaining list of files
    list_of_files = FileList(directory_with_files, string_to_choose_file)
    # Returning only files with 'barcode' string
    list_of_files = [x for x in list_of_files if 'barcode' in x]
    # Creating output dictionary
    to_return_dict = {}
    # Iterating over list of files
    for filename in list_of_files:
        # Obtaining sequencing run date
        sequencing_run_date = filename.split('/')[-1].split('_exon_counts_')[0].split('.')[0]
        # Barcode
        barcode_num = filename.split('/')[-1].split('_exon_counts_')[0].split('.')[1]
        # Convertin sequencing run and barcode number to sample name (used as key in dictionary)
        id_brain_region = dict_mappings[sequencing_run_date][barcode_num]
        # Opening file handle
        f = open(filename, 'r')
        # Iterating over lines in file
        for line in f:
            # Splitting the file by tab
            line_sep = line.split('\t')
            # Exon number
            exon_num  = line_sep[0]
            # Read count - integer
            read_count = int(line_sep[1])
            # Creating exon number sub-dictionary
            if exon_num not in to_return_dict.keys():
                to_return_dict[exon_num] = {}
            # Adding to the dictionary
            to_return_dict[exon_num][id_brain_region] = read_count
            # Deleting variables not needed anymore
            del line_sep, exon_num, read_count
        # Closing file handle
        f.close()
        # Deleting variables not needed anymore
        del line, f, sequencing_run_date, barcode_num, id_brain_region
    # Deleting variables not needed anymore
    del filename, list_of_files
    # Returning output Dictionary
    return to_return_dict

# This function checks if in the list there is at least n individuals
# Input:
#   - input_list: List with sample names
#   - individual_list: List with individual IDs
#   - threshold_individual: Threshold for number of individuals
#   - threshold_tissues: Threshold_for number of tissues per individual
# Output: True of number of individuals is above or equal to threshold, False if it is not above threshold
def IndividualsTissuesThresholdCheck(input_list, individual_list, threshold_individual, threshold_tissues):
    # Starting boolean to return is False
    to_return_bool = False
    # Converting the list to have only individual IDs
    input_individuals = input_list
    input_individuals = [x.split('_')[1] for x in input_individuals]
    # Removing duplicated individual IDs
    seen = set()
    input_individuals = [x for x in input_individuals if x not in seen and not seen.add(x)]
    del seen
    # If there is at least threshold number of individuals, return True
    if len(input_individuals) >= int(threshold_individual):
        to_return_bool = True
    # Else if there are at least another_threshold libraries for the same individual, return True
    else:
        for individual in individual_list:
            list_of_specific_individual = [x for x in input_list if individual in x]
            if len(list_of_specific_individual) >= int(threshold_tissues):
                to_return_bool = True
                break
    return to_return_bool

# This function creates a list of individuals
# Input: Dictionary with all barcode-to-sample-mappings
# Output: List with all unique individual IDs
def IndividualsIDList(dict_mappings):
    # Output list
    output_list = []
    # For key in dictionary
    for sequencing_run in dict_mappings:
        # For all barcodes
        for barcode_id in dict_mappings[sequencing_run]:
            # Append sample name to list
            output_list.append(dict_mappings[sequencing_run][barcode_id])
        del barcode_id
    del sequencing_run
    # Obtain only individual IDs
    output_list = [x.split('_')[1] for x in output_list]
    # Remove duplicates
    seen = set()
    output_list = [x for x in output_list if x not in seen and not seen.add(x)]
    del seen
    # Return list of individuals
    output_list.sort()
    return output_list


# This function filters exons for specified threshold and creates two dictionaries - with included and excluded exons
# Input:
#   - dict_mappings: Dictionary with barcode-to-name mappings
#   - input_dict: Dictionary with all exon read counts
#   - num_of_samples: String - at least this number of individuals must have at least specified number of reads
#   - num_of_tissues: String - at least this number of tissues within the same individual must have at least specified number of reads
#   - num_of_reads: String - exon is valid if it is covered by at least this number of reads
# Output:
#   - include_dict: Dictionary with included exons (above threshold)
#   - exclude_dict: Dictionary with excluded exons (below threshold)
def ReadCountsFilter(dict_mappings, input_dict, num_of_individuals, num_of_tissues, num_of_reads):
    # List of individuals
    list_of_individuals = IndividualsIDList(dict_mappings)
    # Creating dictionary with exons covered by reads above threshold
    include_dict = {}
    # Creating dictionary with exons covered by reads below threshold
    exclude_dict = {}
    # Creating exon number list
    exon_number_list = list(input_dict.keys())
    exon_number_list.sort(key=int)
    # Iterating over exon number list
    for exon_number in exon_number_list:
        # Obtaining list of samples with read counts above threshold
        list_of_samples_above_read_count_threshold = [x for x in input_dict[exon_number].keys() if input_dict[exon_number][x]>=int(num_of_reads)]
        # Obtaining indicator whether an exon fulfill thresholding criteria or not
        to_include_exon = IndividualsTissuesThresholdCheck(list_of_samples_above_read_count_threshold, list_of_individuals, num_of_individuals, num_of_tissues)
        # If any threshold is fulfilled then write a dictionary to include dictionary, else write to dictionary to exclude
        if to_include_exon:
            include_dict[exon_number] = input_dict[exon_number]
        else:
            exclude_dict[exon_number] = input_dict[exon_number]
        # Delete variables not needed anymore
        del list_of_samples_above_read_count_threshold, to_include_exon
    # Delete variables not needed anymore
    del exon_number, exon_number_list, list_of_individuals
    # Return dictionaries with exons above and below threshold, respectively
    return include_dict, exclude_dict

# This function writes output table to the specified output file
# Input:
#   - input_dictionary: Dictionary with read counts for exon and sample
#   - list_of_names_individuals_brain_regions: List of all sample names (sorted)
#   - output_file_path: Path to the output file
# Output: Tab-delimited file with read counts for each sample (sorted by read and by barcode)
def WriteOutputFiles(input_dictionary, list_of_names_individuals_brain_regions, output_file_path):
    # Opening output file handle
    g = open(output_file_path, 'w')
    # Adding header
    header_to_write = '#Exon_number\t' + '\t'.join(list_of_names_individuals_brain_regions) + '\n'
    # Writing header
    g.write(header_to_write)
    # List with all strings
    to_write_list = []
    # Removing header string
    del header_to_write
    # Iterating over list of exons (sorted as integers)
    exon_id_list = list(input_dictionary.keys())
    exon_id_list.sort(key=int)
    for exon_num_id in exon_id_list:
        # Obtaining tab-delimited string with read counts strings
        string_with_read_counts = [str(input_dictionary[exon_num_id][x]) for x in list_of_names_individuals_brain_regions]
        string_with_read_counts = '\t'.join(string_with_read_counts)
        # Final string with exon number
        final_string = str(exon_num_id) + '\t' + string_with_read_counts
        # Append the final string to final list
        to_write_list.append(final_string)
        # Delete variables that are not needed anymore
        del string_with_read_counts, final_string
    # Delete variables that are not needed anymore
    del exon_num_id, exon_id_list
    # Joining the list by new-line character
    to_write_string = '\n'.join(to_write_list)
    del to_write_list
    to_write_string = to_write_string + '\n'
    # Write all the lines
    g.write(to_write_string)
    # Close the file and delete all variables
    g.close()
    del g, to_write_string


# Testing paths
"""
args = parser.parse_args(['/tgac/workarea/group-eg/project_Capture/data/CACNA1C/July_17/CACNA1C_metagene_novel_exons', 'exon_counts_0.7', '2', '2', '5', 'test.include', 'test.exclude'])

"""

if __name__ == '__main__':
    # Import the input files using argparse module
    parser = argparse.ArgumentParser(description='This script obtains a list of exons expressed in i individuals or j brain tissues of the same individual, at the level of at least k reads.')
    parser.add_argument('in_barcode_sample_name', type=str, help='Path to the file with barcode-to-sample-name mappings.')
    parser.add_argument('in_dir', type=str, help='Path to the directory where files with exon counts are stored.')
    parser.add_argument('in_distinct_string', type=str, help='String allowing to distinguish files of interest from other files.')
    parser.add_argument('num_individuals', type=str, help='String with number of individuals (i).')
    parser.add_argument('num_libraries', type=str, help='String with number of brain regions (libraries) for the same individual (j).')
    parser.add_argument('read_count_threshold', type=str, help='String with read count threshold (k).')
    parser.add_argument('out_file_included', type=str, help='Path to the output file containing exon counts for exons fulfilling filtering criteria.')
    parser.add_argument('out_file_excluded', type=str, help='Path to the output file containing exon counts for exons removed after filtering.')
    args = parser.parse_args()
    # Create dictionary with barcode-to-name mappings and list with sample names
    barcode_sample_names_dict, sample_name_list_sorted = IndividualBrainRegionMatcherDict(args.in_barcode_sample_name)
    # Obtaining dictionary with all samples as keys, exon numbers as sub-keys and read counts as values
    read_count_dict = ReadExonsCountsMatchedSampleNames(args.in_dir, args.in_distinct_string, barcode_sample_names_dict)
    # Creating dictionary with exons above and below threshold
    dict_exons_above_threshold, dict_exons_below_threshold = ReadCountsFilter(barcode_sample_names_dict, read_count_dict, args.num_individuals, args.num_libraries, args.read_count_threshold)
    del read_count_dict
    gc.collect()
    # Writing output table for exons above threshold
    WriteOutputFiles(dict_exons_above_threshold, sample_name_list_sorted, args.out_file_included)
    # Writing output table for exons below threshold
    WriteOutputFiles(dict_exons_below_threshold, sample_name_list_sorted, args.out_file_excluded)
    # Deleting the working space variables
    del sample_name_list_sorted, dict_exons_above_threshold, dict_exons_below_threshold
    gc.collect()
