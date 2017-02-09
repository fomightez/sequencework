#! /usr/bin/env python

# extract_conserved_from_alignments.py by Damian Menning and Wayne Decatur
# ver 0.2
#
#*******************************************************************************
# USES Python 2.7
# PURPOSE: Takes a file of aligned sequences in FASTA format and extracts the
# residues that are 100% conserved in all sequences. As it scans the sequences,
# it writes to a file the residue if it is the same or a blank line if it is not.
#
#
#
# Dependencies:
# Biopython
#
# v.0.1. Basics all in place.
# v.0.2. Added a setting for a minimum requirement for the length of the
# consensus in order for a block of consensus residues to get printed to
# the output file.
#
# To do:
#
#
# TO RUN:
# Example,
# Enter on the command line of your terminal, the line
#-----------------------------------
# python extract_conserved_from_alignments.py
#-----------------------------------
#
#
#*******************************************************************************
##################################
#  USER ADJUSTABLE VALUES        #
file_name = "test_alignment.fas"

minimum_length_of_consensus_block = 1 # won't add sequence ragment to the output
# file unless matches in a series at least this length. Default is 1 so that
# each block of matching sequences, no matter how short, are saved to output.

##################################
#

#
#*******************************************************************************
#**********************END USER ADJUSTABLE VARIABLES****************************










#*******************************************************************************
#*******************************************************************************
###DO NOT EDIT BELOW HERE - ENTER ANY VALUES ABOVE###

import os
import sys
from Bio import SeqIO





###---------------------------HELPER FUNCTIONS---------------------------------###

def generate_output_file_name(file_name):
    '''
    Takes a file name as an argument and returns string for the name of the
    output file. The generated name is based on the original file
    name.

    Specific example
    ================
    Calling function with
        ("Gomphonemataceae 18S cut.fas")
    returns
        "Gomphonemataceae 18S cut_conserved.fas"
    '''
    main_part_of_name, file_extension = os.path.splitext(
        file_name) #from http://stackoverflow.com/questions/541390/extracting-extension-from-filename-in-python
    if '.' in file_name:  #I don't know if this is needed with the os.path.splitext method but I had it before so left it
        if minimum_length_of_consensus_block != 1:
            return main_part_of_name + "_conserved_" + str(
                minimum_length_of_consensus_block) +"bps" + file_extension
        else:
            return main_part_of_name + "_conserved" + file_extension
    else:
        if minimum_length_of_consensus_block != 1:
            return file_name + "_conserved_" + str(
                minimum_length_of_consensus_block) +"bps.fas"
        else:
            return file_name + "_conserved.fas"



def make_dict_of_alignments(fasta_seq_file):
    '''
    takes a file of aligned sequence records in fasta format and reads each one
    placing itin a dictionary for easy access later using the id as a key.
    returns the dictionary it produced.
    also returns a string of the first sequence identified to use as the measure
    of length and for establishing a default consensus sequence for comparison.
    '''
    # initialize needed variables
    aligned_seq_dict = {}
    example_sequence = "NotDefined"
    for seq_record in SeqIO.parse(fasta_seq_file, "fasta"):
        # extract the needed info and make dictionary
        name, sequence = str(seq_record.id), str(seq_record.seq)
        aligned_seq_dict[name] = sequence
        # grab first aligned sequence encountered as representative
        if example_sequence == "NotDefined":
            example_sequence = str(seq_record.seq)
    return aligned_seq_dict, example_sequence

###--------------------------END OF HELPER FUNCTIONS---------------------------###











###-------------------------Main body of script------------------------------###
###----------------------GET FILE AND PREPARE TO PARSE-----------------------###

# generate output file name based on input file name
output_file_name = generate_output_file_name(file_name)
# inititalize output file handler
output_file_handler = open(output_file_name, "w")

#give user some stats and feeback
sys.stderr.write("Reading in your file...")

# Need dictionary of the aligned sequences for accessing and will use
# a representatice sequence from there to extract details such as length and
# default consensus.
alignments_dict, representative_aligned_seq = make_dict_of_alignments(file_name)

#initialize a string that will be the current sequence block
current_sequence_for_possible_output = ""

# initialize a string for storing the consensus. This is optional in addition to
# the outout file being made
consensus_string = ""


# Iterate over the representative sequence string using each character for
# comparison. Idx holds the current index position and will make calling each
# position in each sequence easy.
# NOTE that the `else` is indented correctly and belongs to the `for` loop. It
# will only be performed if the for loop completes without a `break` as
# described at https://docs.python.org/2/tutorial/controlflow.html
for idx, extracted_residue in enumerate(representative_aligned_seq):
    for each_key in alignments_dict:
        if alignments_dict[each_key][idx] != extracted_residue:
            # Now that chain of matches has ended, evaluate if chain series long
            # enough to save to output file
            if (len(current_sequence_for_possible_output) >= minimum_length_of_consensus_block):
                #write conserved series to output file
                output_file_handler.write(current_sequence_for_possible_output)
                # A line break is written to the line output file
                output_file_handler.write("\n")
            #reset the string collection for blocks to save to the output file
            current_sequence_for_possible_output = ""
            # consensus_string gets a blank character added to it
            consensus_string += " "
            break
    else:
        # if matches in all sequences, save to a growing sequence that may possibly
        # be written to output if length of sequence chunck meets minimum length
        # when the chain of matches ends.
        current_sequence_for_possible_output += extracted_residue
        # consensus string gets this residue added to it
        consensus_string += extracted_residue

# close the file handling stream
output_file_handler.close()


#give user some stats and feeback
sys.stderr.write("\nConcluded. \n")
sys.stderr.write("\nA file named '"+
    output_file_name+"' was saved with the results "+
    "\nin same directory as the input file.\n\n")
sys.stderr.write("\nConsensus sequence = \n"+
    consensus_string +".\n\n")



#*******************************************************************************
###-***********************END MAIN PORTION OF SCRIPT***********************-###
#*******************************************************************************
