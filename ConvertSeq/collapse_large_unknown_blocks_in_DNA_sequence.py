#!/usr/bin/env python
# collapse_large_unknown_blocks_in_DNA_sequence.py
__author__ = "Wayne Decatur" #fomightez on GitHub
__license__ = "MIT"
__version__ = "0.1.0"


# collapse_large_unknown_blocks_in_DNA_sequence.py by 
# Wayne Decatur
# ver 0.1
#
#*******************************************************************************
# Verified compatible with both Python 2.7 and Python 3.6; written initially in 
# Python 3. 
#
#
# PURPOSE: Takes a DNA sequence in FASTA format and reduces the large blocks of
# unknown nucleotides, represented by repated `N`s, to be shorter. This is meant 
# to be used for sequences that appear to have arbitrarily-sized blocks of 
# unknown nucleotides so the blocks are short and easier to align / or visualize 
# how the sub-elements should align in multiple sequence alignment in cases 
# where it looks like the block may not be as large as arbitarily made. (Or to
# help assess if that may be the case.)  
# The output sequence remains in FASTA format.
# Also works if the provided FASTA file is a multi-FASTA, i.e., contains 
# multiple sequences in FASTA format in the one file. All sequences in the file
# will have the large blocks of repeated Ns collapsed to a small size. Typcially
# this would then be submitted to software that performs a MULTIPLE SEQUENCE
# ALIGNMENT to see if some conerved elements are better aligned now, indicating
# the perhaps arbitrary blocks were making it hard to see the adjacent placement
# while keeping in mind these new small blocks of `N`s are indeed 'fuzzy'.
#
#
#
# Written to run from command line or imported into/pasted/loaded inside a 
# Jupyter notebook cell.
#
#
#
#
#
#
#
#
#
#
#
#
#
# Dependencies beyond the mostly standard libraries/modules:
# Biopython
#
#
# VERSION HISTORY:
# v.0.1. basic working version.
#
#
# To do:
#
#
#
#
# TO RUN:
# Examples,
# Enter on the command line of your terminal, the line
#-----------------------------------
# python collapse_large_unknown_blocks_in_DNA_sequence.py sequence.fa
#-----------------------------------
#
# Issue `collapse_large_unknown_blocks_in_DNA_sequence.py -h` for details.
# 
#
#
# To use this after importing/pasting or loading into a cell in a Jupyter 
# notebook, specify at least the sequence file (or results as a string) in the 
# call to the main function similar to below:
# collapse_large_unknown_blocks_in_DNA_sequence("sequence.fa")
#-or-
# To specify OPTIONAL different number to collapse the large blocks of unknown 
# nucleotides down to, specify  `num_of_repeatedNs_to_collapse_to` like on next 
# line when calling the main function:
# collapse_large_unknown_blocks_in_DNA_sequence("sequence.fa",num_of_repeatedNs_to_collapse_to=10)
#-or-
# To specify OPTIONAL suffix to add to the names of generated files, specify 
# `suffix_for_saving` like on next line when providing alignment as string:
# extract_regions_from_clustal_alignment("sequence.fa",suffix_for_saving="_condensed")
# 
#
# 
#
'''
CURRENT ACTUAL CODE FOR RUNNING/TESTING IN A NOTEBOOK WHEN IMPORTED/LOADED OR 
PASTED IN ANOTHER CELL:
collapse_large_unknown_blocks_in_DNA_sequence("sequence.fa")
'''
#
#
#*******************************************************************************
#





#*******************************************************************************
##################################
#  USER ADJUSTABLE VALUES        #

##################################
#

num_of_repeatedNs_to_collapse_to = 5 # Number of repeated `N` characters in a 
# row to collapse down to. For example, if set to `5` than collapse to five if 
# have a block of six or more unknown nucleotides in a row. Can be changed when
# script or main function called

suffix_for_saving = "_col" #to be used for naming the output automatically
# when running script from command line to act on an input file

#
#*******************************************************************************
#**********************END USER ADJUSTABLE VARIABLES****************************






















#*******************************************************************************
#*******************************************************************************
###DO NOT EDIT BELOW HERE - ENTER VALUES ABOVE###

import sys
import os
import re
from Bio import SeqIO
from Bio.Seq import Seq 
from Bio.SeqRecord import SeqRecord 
from Bio.Alphabet import generic_dna



###---------------------------HELPER FUNCTIONS-------------------------------###


def generate_output_file_name(file_name,suffix_for_saving):
    '''
    Takes a file name as an argument and returns string for the name of the
    output file. The generated name is based on the original file
    name.


    Specific example
    =================
    Calling function with
        ("sequence.fa", "_col")
    returns
        "sequence_col.fa"
    '''
    main_part_of_name, file_extension = os.path.splitext(
        file_name) #from 
    #http://stackoverflow.com/questions/541390/extracting-extension-from-filename-in-python
    if '.' in file_name:  #I don't know if this is needed with the os.path.splitext method but I had it before so left it
        return main_part_of_name + suffix_for_saving  + file_extension
    else:
        return file_name + suffix_for_saving + ".fa"

def get_start_n_ends_for_match_to_pattern(pattern_obj,a_string):
    '''
    Takes a compiled pattern object and returns a list of the start 
    and end indices in tuple form of the occurences of the pattern in
    the provided string.
    
    Returns a list of tuples where first value is the start and
    second is the end of the span of the pattern match.
    '''
    start_end_tuples = []
    for m in pattern_obj.finditer(a_string):
        start_end_tuples.append(m.span()) # from https://stackoverflow.com/a/250306/8508004
    return start_end_tuples

###--------------------------END OF HELPER FUNCTIONS-------------------------###
###--------------------------END OF HELPER FUNCTIONS-------------------------###



#*******************************************************************************
###------------------------'main' function of script--------------------------##

def collapse_large_unknown_blocks_in_DNA_sequence(sequence,
    num_of_repeatedNs_to_collapse_to = num_of_repeatedNs_to_collapse_to,
    suffix_for_saving = suffix_for_saving):
    '''
    Main function of script. 
    Takes a DNA sequence in FASTA format and reduces the large blocks of 
    unknown sequences, represented by N, to be shorter.
    The sequence remains in FASTA format.

    `num_of_repeatedNs_to_collapse_to` and `suffix_for_saving` optional when
    calling the main function from a Jupyter cell after import or pasting in 
    another cell.
    '''
    # make a records list because individual entries in this will get altered 
    # later as needed
    records = []  #intialize a list of records for the present FASTA file
    for record in SeqIO.parse(sequence, "fasta"):
        records.append(record)

    # make a pattern object since it will be used often
    pat_obj = re.compile(
        "N{{{},}}".format(num_of_repeatedNs_to_collapse_to+1), re.I) # this will 
    # match greedily to repeats that are of length one more than 
    # `num_of_repeatedNs_to_collapse_to` or greater; based on on page 
    # 469 (appendix 2) of Practical Computing for Biologists;
    # bracket handling for combining regex use of brackets to control numbers to 
    # match combined with use of `.format()` to provide the number as a variable 
    # comes from https://stackoverflow.com/a/5466478/8508004  

    # go through parsed records and act if necessary
    modf_seq = None #reset here allows use in notebook cell when running that 
    # cell multiple times
    for indx, record in enumerate(records):
        # first check if any occurences of more than the 
        # `num_of_repeatedNs_to_collapse_to` unknown nucleotides (N's) in a row. 
        # If so,then process.
        if "N"*(num_of_repeatedNs_to_collapse_to+1) in record.seq:
            # process seq since it has at least one instance needing collapsing
            unk_block_locations = get_start_n_ends_for_match_to_pattern(
                pat_obj,str(record.seq))
            #print(record.id) # ONLY FOR DEBUGGING
            #print (unk_block_locations) # ONLY FOR DEBUGGING
            
            # Initially to be able to modify SeqRecord / Seq objects I learned 
            # how to go from a Seq to a MutableSeq object and back to Seq and 
            # replace the record with the update Seq, see 
            # `working_out_removing_Nblocks_from_fasta_with_biopython.ipynb` and 
            # `replacing and editing sequences with biopython.md`. And I thought 
            # I'd need that approach here but I think think here since not 
            # truncating / slicing or replacing certain nucleotide or range in 
            # one step, best to do modifying at the string level and then just 
            # convert the finish string back to a Seq object in one step.
            
            # use the block locations to rebuild the string, collapsing where 
            # appropriate. Process the blocks in reverse order to how they are 
            # in the sequence so that the indices for blocks earlier in the 
            # sequence aren't changed as parts of the string are recombined.
            modf_seq = str(record.seq)
            # Verified this even works when a block of unknowns is very start of 
            # sequence.
            for start_n_ends in reversed(unk_block_locations):
                #print(start_n_ends[0]) # ONLY FOR DEBUGGING
                #print (seq_to_mod[:start_n_ends[0]]) # ONLY FOR DEBUGGING
                modf_seq = (modf_seq[:start_n_ends[0]] + 
                              "N"*num_of_repeatedNs_to_collapse_to + 
                              modf_seq[start_n_ends[1]:]) # based on
                              # https://stackoverflow.com/a/41753038/8508004
            # replace the appropriate record with the modified sequence, tagging
            # the description line
            records[indx] = SeqRecord(
                Seq(modf_seq, generic_dna), 
                id=record.id, description=record.description+" COLLAPSED")#based
            # on https://www.biostars.org/p/48797/ and `.ungap()` method, see
            # https://github.com/biopython/biopython/issues/1511 , and `description`
            # from what I've seen for `id` plus https://biopython.org/wiki/SeqIO
            #print (records[indx]) # ONLY FOR DEBUGGING

    # Replace the FASTA file with the modified records
    output_file_name = generate_output_file_name(sequence,suffix_for_saving)
    SeqIO.write(records,output_file_name, "fasta");
    # Feedback
    sys.stderr.write("\n\n*****************DONE**************************\n")
    if modf_seq:
        sys.stderr.write("Long blocks of unknown nucleotides collapsed to "
            "length of {}, and \nsaved in FASTA format as '{}'.".format(
            num_of_repeatedNs_to_collapse_to, output_file_name))
    else:
        sys.stderr.write("No blocks of unknown nucleotides of length greater "
            "than {} \ndetected in {}, and so no changes made.".format(
            num_of_repeatedNs_to_collapse_to, sequence))
    sys.stderr.write("\n*****************DONE**************************\n")



###--------------------------END OF MAIN FUNCTION----------------------------###
###--------------------------END OF MAIN FUNCTION----------------------------###










#*******************************************************************************
###------------------------'main' section of script---------------------------##
def main():
    """ Main entry point of the script """
    # placing actual main action in a 'helper'script so can call that easily 
    # with a distinguishing name in Jupyter notebooks, where `main()` may get
    # assigned multiple times depending how many scripts imported/pasted in.
    kwargs = {}
    kwargs['num_of_repeatedNs_to_collapse_to']= num_of_repeatedNs_to_collapse_to
    kwargs['suffix_for_saving'] = suffix_for_saving
    collapse_large_unknown_blocks_in_DNA_sequence(sequence,**kwargs)
    # using https://www.saltycrane.com/blog/2008/01/how-to-use-args-and-kwargs-in-python/#calling-a-function
    # to build keyword arguments to pass to the function above
    # (see https://stackoverflow.com/a/28986876/8508004 and
    # https://stackoverflow.com/a/1496355/8508004 
    # (maybe https://stackoverflow.com/a/7437238/8508004 might help too) for 
    # related help). Makes it easy to add more later.





if __name__ == "__main__" and '__file__' in globals():
    """ This is executed when run from the command line """
    # Code with just `if __name__ == "__main__":` alone will be run if pasted
    # into a notebook. The addition of ` and '__file__' in globals()` is based
    # on https://stackoverflow.com/a/22923872/8508004
    # See also https://stackoverflow.com/a/22424821/8508004 for an option to 
    # provide arguments when prototyping a full script in the notebook.
    ###-----------------for parsing command line arguments-------------------###
    import argparse
    parser = argparse.ArgumentParser(prog=
        'collapse_large_unknown_blocks_in_DNA_sequence.py',
        description="collapse_large_unknown_blocks_in_DNA_sequence.py \
        Takes a DNA sequence in FASTA format and reduces the large blocks of \
        unknown nucleotides, represented by numerous uninterrupted `N`s, to be \
        shorter. This is meant to \
        be used for sequences that appear to have arbitrarily-sized blocks of \
        unknown nucleotides so the blocks are short and easier to align or to \
        visualize how they should align in \
        multiple sequence alignment, for cases where it looks like the block \
        may not be as large as arbitarily made. Or to help assess if that is \
        indeed the case. \
        **** Script by Wayne Decatur   \
        (fomightez @ github) ***")

    parser.add_argument("sequence_file", help="Name of sequence file to \
        use as input to produce a file where the repeated unknown \
        nucleotidess, such as `NNNNNNNNN`, are reduced to `NNNNN`. Can be a \
        multi-FASTA file, i.e., multiple sequences in FASTA format in one \
        file. All sequences included will be modified\
        .", metavar="SEQUENCE_FILE")

    parser.add_argument('-ct', '--collapse_to', action='store', type=int, 
    default= num_of_repeatedNs_to_collapse_to, help="OPTIONAL: Adjust the \
    number of repeated `N`s to collapse down to. Any block larger than \
    that will be represented as this number of `N` in a row. \
    If none provided, '{0}' will be used, meaning the blocks larger than {0} \
    will become '{1}'.".format(num_of_repeatedNs_to_collapse_to,
    "N"*num_of_repeatedNs_to_collapse_to))
    
    parser.add_argument('-os', '--output_suffix', action='store', type=str, 
    default= suffix_for_saving, help="OPTIONAL: Set a suffix for including in \
    file name of output. \
    If none provided, '{}' will be used.".format(suffix_for_saving))




    #I would also like trigger help to display if no arguments provided because 
    # need at least one for url
    if len(sys.argv)==1:    #from http://stackoverflow.com/questions/4042452/display-help-message-with-python-argparse-when-script-is-called-without-any-argu
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    sequence = args.sequence_file
    num_of_repeatedNs_to_collapse_to = args.collapse_to
    suffix_for_saving = args.output_suffix


    main()

#*******************************************************************************
###-***********************END MAIN PORTION OF SCRIPT***********************-###
#*******************************************************************************
