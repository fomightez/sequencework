#!/usr/bin/env python
# get_specified_length_of_end_of_seq_from_FASTA.py
__author__ = "Wayne Decatur" #fomightez on GitHub
__license__ = "MIT"
__version__ = "0.1.0"


# get_specified_length_of_end_of_seq_from_FASTA.py by 
# Wayne Decatur
# ver 0.1
#
#*******************************************************************************
# Verified compatible with both Python 2.7 and Python 3.6; written initially in 
# Python 3. 
#
#
# PURPOSE: Takes a sequence file (FASTA-format), & a record id, and a number
# (integer), and extracts a sequence of specified length from the end of the
# indicated sequence. The number provided is what specifies the length 
# extracted. (The FASTA-formatted sequence file is assumed by default 
# to be a multi-FASTA, i.e., multiple sequences in the provided file, although 
# it definitely doesn't have to be. In case it is only a single sequence, the
# record id becomes moot, see below.) 
# A sequence string of the specified length will be returned. Redirect the 
# output to a file if using command line version and a file is needed.
#
# Note that if there is only one record in the specified sequence file, the 
# record id is moot and you can instead provide any string for that parameter 
# as it will be ignored. This makes the script more flexible in cases where 
# sequence files aren't complex as the user doesn't need to provide an actual 
# record id.
# 
# This script is meant to be used after you have performed a large alignment, 
# say of an entire chromosome, in order to have individual occurrences of 
# related segments fall linearly with where they match up along the span of the 
# sequence. Often due to large (seeming-to-be) arbitratrily-sized blocks of 
# repeated unknown nucleotides (which are often good to 'collapse', see 
# `collapse_large_unknown_blocks_in_DNA_sequence.py`) the 'ends' of regions 
# often fail to get extracted exactly right. This script is meant to help in 
# determining how best to clean up such instances. For example, in the obtained
# sequence, is there an 'end' that matches up better with the pattern of known
# 'ends' and should be added?
#
# It simply uses Python string slicing syntax to get the specified number of
# residues at the end of a sequence and so if that sequence has gap indicators
# such as "-", those will be included in the extracted sequence and count 
# towards the number of characters.
#
#
#
#
#
# 
#
# Written to run from command line or imported into/pasted/loaded inside a 
# Jupyter notebook cell. 
#
#
#
# Based on my script `get_seq_following_seq_from_FASTA.p` so it will inherit 
# any issues/bugs that one has.
#
#
#
#
#
#
#
# If you are Wayne, see  `GSD Extracting_XXXs_from_XXXX_orthologs.ipynb` 
# for impetus behind this script.
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
# To do:
#
#
#
#
# TO RUN:
# Examples,
# Enter on the command line of your terminal, the line
#-----------------------------------
# python get_specified_length_of_end_of_seq_from_FASTA.py seq.fa record_id amount_to_get > extracted_seq.fa
#-----------------------------------
#
# Issue `python get_specified_length_of_end_of_seq_from_FASTA.py -h` for details.
# 
#
#
# To use this after importing/pasting or loading into a cell in a Jupyter 
# notebook, call the main function similar to below:
# following_seq = get_specified_length_of_end_of_seq_from_FASTA("seq.fa", "Skluv" , 200)
# -or-
# when calling the main function you can include `filter_dashes = False` to 
# not remove the dashes from the supplied pattern sequence before comparing.
# This is the equivalent to using `--leave_dashes` option when calling the 
# script from the command line.
#
# 
#
'''
CURRENT ACTUAL CODE FOR RUNNING/TESTING IN A NOTEBOOK WHEN IMPORTED/LOADED OR 
PASTED IN ANOTHER CELL:
following_seq = get_specified_length_of_end_of_seq_from_FASTA("seq.fa", "Skluv", 200)
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




###---------------------------HELPER FUNCTIONS-------------------------------###





###--------------------------END OF HELPER FUNCTIONS-------------------------###
###--------------------------END OF HELPER FUNCTIONS-------------------------###



#*******************************************************************************
###------------------------'main' function of script--------------------------##

def get_specified_length_of_end_of_seq_from_FASTA(
    sequence_file, record_id, amount_to_get):
    '''
    Main function of script.
    Takes a sequence string, a sequence file (FASTA-format), and a record id and 
    extracts from the specified amount to get (provided as an integer) from the 
    end of the indicated sequence entry.
    The FASTA-formatted sequence file is assumed by default to be a 
    multi-FASTA, i.e., multiple sequences in the provided file, although it 
    definitely doesn't have to be. In case it is only a single sequence, the 
    record id becomes moot and users can provide anything for this parameter. 

    A sequence string of the specified length will be returned.
    '''
    # make a records list because need to adjust handling if only one record in 
    # sequence file because don't want record_id to matter in such a case & 
    # be sure to give feedback of that sort.
    records = []  #intialize a list of records for the FASTA file
    for record in SeqIO.parse(sequence_file, "fasta"):
        records.append(record)
    single_record = False
    if len(records) == 1:
        record_id = records[0].id
        single_record = True
        # feedback
        sys.stderr.write("Single sequence with id of '{}' provided in the "
            "sequence file.\nIt will be "
            "used to extract the last "
            "{} characters in it.\n\n".format(record_id, amount_to_get))
    
    # If more than one, go through parsed records and collect the record to act 
    # on.
    record_to_mine  = None
    if single_record:
        record_to_mine = records[0]
    else:
        for record in records:
            if record.id == record_id:
                record_to_mine = record

    # Now extract the sequence of specified length from the record.
    if len(record_to_mine) < amount_to_get:
        sys.stderr.write("Note that the sepecified number of residues "
        "to get, {}, exceeds the\nlength of the specified record, which is {} "
        "residues in length.\nThe entire sequence has been "
        "returned.".format(amount_to_get, len(record_to_mine)))
        return str(record_to_mine.seq)
    else:
        return str(record_to_mine.seq[-amount_to_get:])

    



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
    #kwargs['another_arg'] = another_arg
    result = get_specified_length_of_end_of_seq_from_FASTA(
        sequence_file, record_id, amount_to_get,**kwargs)
    # using https://www.saltycrane.com/blog/2008/01/how-to-use-args-and-kwargs-in-python/#calling-a-function
    # to build keyword arguments to pass to the function above
    # (see https://stackoverflow.com/a/28986876/8508004 and
    # https://stackoverflow.com/a/1496355/8508004 
    # (maybe https://stackoverflow.com/a/7437238/8508004 might help too) for 
    # related help). Makes it easy to add more later.
    print (result) #redirect stdout if need a file is needed






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
        'get_specified_length_of_end_of_seq_from_FASTA',
        description="get_specified_length_of_end_of_seq_from_FASTA \
        takes a sequence file (FASTA-format), & a record id, and a number \
        (integer), and extracts a sequence of specified length from the end of \
        the indicated sequence. The number provided is what specifies the \
        length extracted. (The FASTA-formatted \
        sequence file is assumed by default to be \
        a multi-FASTA, i.e., multiple sequences in the provided file, \
        although it definitely doesn't have to be. In case it is only a \
        single sequence, the record id becomes moot, see below.) A sequence \
        string of the specified length will be returned. Redirect the output \
        to a file if that is what is needed. \
        **** Script by Wayne Decatur   \
        (fomightez @ github) ***")

    parser.add_argument("sequence_file", help="Name of sequence file to \
        use as input. Must be FASTA format. Can be a \
        multi-FASTA file, i.e., multiple sequences in FASTA format in one \
        file.", metavar="SEQUENCE_FILE")
    parser.add_argument("record_id", help="Specific identifier of sequence \
        entry in sequence file to mine. If the provided sequence file only \
        contains one sequence, that sequence will be mined and what is provided \
        for this parameter will be ignored. In other words, if the sequence \
        file is not a multi-FASTA file, you don't need to determine the \
        identifier and can instead just enter `blahblah` or any other \
        nonsensical string in this spot.", metavar="RECORD_ID")
    parser.add_argument("amount_to_get", type=int, help="Number (integer) of \
        residues \
        to retrieve from end. ", 
        metavar="NUMBER_TO_GET")






    #I would also like trigger help to display if no arguments provided because 
    # need at least one for url
    if len(sys.argv)==1:    #from http://stackoverflow.com/questions/4042452/display-help-message-with-python-argparse-when-script-is-called-without-any-argu
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    sequence_file = args.sequence_file
    record_id = args.record_id
    amount_to_get = args.amount_to_get


    main()

#*******************************************************************************
###-***********************END MAIN PORTION OF SCRIPT***********************-###
#*******************************************************************************
