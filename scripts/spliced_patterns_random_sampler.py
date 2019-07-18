#!/usr/bin/env python3

# -*- coding: utf-8 -*-


# Importing modules
import sys
import argparse
import os
import random

# This function is to read all sequences into list of striped lines
# Input: Path to the input file
# Output: List with striped lines
def CreatingListAllLines(input_path):
    # Open the file handle
    f = open(input_path, 'r')
    # Read all rows to the list
    list_to_return = f.readlines()
    # Strip all rows
    list_to_return = [x.strip() for x in list_to_return]
    # Close the file and delete unnecessary variables
    f.close()
    del f
    # Return list
    return list_to_return

# This function generates random list of lines (of length equal to num_elements)
# Input:
#   - input_list: List with all striped lines
#   - num_elements: number of elements to obtain
# Output: List with 'num_elements' random lines
def ListLinesRandomizer(input_list, num_elements):
    # Randomly choose num_elements items from list
    to_return_list = random.sample(input_list, int(num_elements))
    # Return the list
    return to_return_list

# This function writes n reads with splice patterns to the output file
# Input:
#   - input_list: List with randomly chosen lines
#   - output_file_path: Path to the output file
def WriteOutputFile(input_list, output_file_path):
    # Open output file handler
    g = open(output_file_path, 'w')
    # Join all elements of list by newline character and add newline character at the end of string
    string_to_write = '\n'.join(input_list)
    string_to_write = string_to_write + '\n'
    # Write the string to output file
    g.write(string_to_write)
    # Close file handle
    g.close()
    # Remove all the variables
    del string_to_write, g


# Test paths
"""
in_file = "2017_01_13_CACNA1C_barcode09.fa"
n_items = '1760'
out_file = "test.test"
"""

if __name__ == '__main__':
    # Import the input files using argparse module
    parser = argparse.ArgumentParser(description='This script is used to randomly choose n sequences from FASTA file and to write it to the output file.')
    parser.add_argument("in_file", type=str, help="Path to the input FASTA file")
    parser.add_argument("n_items", type=str, help="Number of sequences to report")
    parser.add_argument("out_file", type=str, help="Path to the output FASTA file with randomly chosen n sequences")
    args = parser.parse_args()
    # Creating input dictionary with all FASTA sequences
    list_all_lines = CreatingListAllLines(args.in_file)
    # Creating list with n keys
    rand_lines_list = ListLinesRandomizer(list_all_lines, args.n_items)
    # Creating output file
    WriteOutputFile(rand_lines_list, args.out_file)
