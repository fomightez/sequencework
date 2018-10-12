#!/usr/bin/env python
# get_seq_following_seq_from_multiFASTA.py
__author__ = "Wayne Decatur" #fomightez on GitHub
__license__ = "MIT"
__version__ = "0.1.0"


# get_seq_following_seq_from_multiFASTA.py by 
# Wayne Decatur
# ver 0.1
#
#*******************************************************************************
# Verified compatible with both Python 2.7 and Python 3.6; written initially in 
# Python 3. 
#
#
# PURPOSE: Takes a sequence pattern string, a sequence file (FASTA-format), & a 
# record id, and extracts a sequence of specified size following the sequence 
# pattern. (The FASTA-formatted sequence file is assumed by default to be a 
# multi-FASTA, i.e., multiple sequences in the provided file, although it 
# definitely doesn't have to be. In case it is only a single sequence, the
# record id becomes moot, see below.) 
# A sequence string of the specified length will be returned. Redirect the 
# output to a file if using command line version and a file is needed.
#
# The provided sequence pattern will be matched regardless of case, as both the
# input sequence and pattern to search will be converted to lowercase. Beyond 
# being insensitive of the case, REGULAR EXPRESSION SEARCH TERM SYNTAX IS 
# ACCEPTABLE in the provided sequence pattern.
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
# sequence. Often due to  large (seeming-to-be) arbitratrily-sized blocks of 
# repeated unknown nucleotides (which are often good to 'collapse', see 
# `collapse_large_unknown_blocks_in_DNA_sequence.py`) the 'ends' of regions 
# often fail to get extracted exactly right. This script is meant to help in 
# determining how best to clean up such instances. For example, in the obtained
# sequence, is there an 'end' that matches up better with the pattern of known
# 'ends' and should be added?
#
# It is designed to handle/filter gaps ('dashes') in the provided sequence 
# patterns. The idea being that the known sequence ends may be manually 
# extracted from sequence alignments. This way the user is not wasting time 
# removing the gap indications / dashes from the collected text lines. The 
# default handling of removing the gaps to ignore them can be overriden. 
# The idea is that maybe you'll have a multiple sequence alignment 
# file saved as FASTA with dashes, i.e., aligned FASTA file format and may want 
# to use this script.  (The caveat is that number of residues to get will then 
# be counting the gaps / dashes too.)
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
#
#
#
#
#
#
#
#
# If you are Wayne, see 
# `Gathering XXXXXXX from aligned XXXXXXXXX genomes collected from the XXXXX genomes project.ipynb` 
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
# - check regex works
#
#
#
# TO RUN:
# Examples,
# Enter on the command line of your terminal, the line
#-----------------------------------
# python get_seq_following_seq_from_multiFASTA.py seq.fa record_id pattern amount_to_get > extracted_seq.fa
#-----------------------------------
#
# Issue `get_seq_following_seq_from_multiFASTA.py -h` for details.
# 
#
#
# To use this after importing/pasting or loading into a cell in a Jupyter 
# notebook, specify at least the results file (or results as a string) in the 
# call to the main function similar to below:
# following_seq = get_seq_following_seq_from_multiFASTA("seq.fa", "Skluv", "GAAATTTCCCCCAAAATGT" , 200)
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
following_seq = get_seq_following_seq_from_multiFASTA("seq.fa", "Skluv", "GAAATTTCCCCCAAAATGT" , 200)
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



def get_start_n_ends_for_match_to_pattern(pattern_obj,a_string):
    '''
    Takes a compiled pattern object and returns a list of the start 
    and end indices in tuple form of the occurences of the pattern in
    the provided string.
    
    Returns a list of tuples where first value is the start and
    second is the end of the span of the pattern match.
    '''
    start_end_tuples = []
    for m in pattern_obj.finditer(a_string.lower()):
        start_end_tuples.append(m.span()) # from https://stackoverflow.com/a/250306/8508004
    return start_end_tuples

###--------------------------END OF HELPER FUNCTIONS-------------------------###
###--------------------------END OF HELPER FUNCTIONS-------------------------###



#*******************************************************************************
###------------------------'main' function of script--------------------------##

def get_seq_following_seq_from_multiFASTA(
    sequence_file, record_id, seq_to_find, amount_to_get, filter_dashes = True):
    '''
    Main function of script.
    Takes a sequence string, a sequence file (FASTA-format), and a record id and 
    extracts a sequence of specified size following the sequence string. 
    (The FASTA-formatted sequence file is assumed by default to be a 
    multi-FASTA, i.e., multiple sequences in the provided file, although it 
    definitely doesn't have to be. In case it is only a single sequence, the 
    record id becomes moot and users can provide anything for this parameter. 

    A sequence string of the specified length will be returned.
    '''
    # make a records list because need to adjust record_id because don't 
    # xpect it to matter in such a case of only one entry in sequence file & 
    # give feedback of that sort.
    records = []  #intialize a list of records for the present FASTA file
    for record in SeqIO.parse(sequence_file, "fasta"):
        records.append(record)
    single_record = False
    if len(records) == 1:
        record_id = records[0].id
        single_record = True
        # feedback
        sys.stderr.write("Single sequence with id of {} provided in the "
            "sequence file. It will be "
            "used to search for the provided sequence pattern and provide the "
            "{} residues after it.".format(
            record_id, amount_to_get))

    # Prepare the pattern to be searched for.
    if filter_dashes:
        seq_to_find = seq_to_find.replace("-","")
    # make a pattern object since it could be used often if I change things to 
    # search the file for multiple records instead of just one at a time
    pat_obj = re.compile(seq_to_find.lower())
    # luckily it seemed it didn't (and wouldn't work) to use 
    # https://stackoverflow.com/a/12989308/8508004 to bring regular expression
    # search term from shell. Seems like according to what I found at while 
    # working on `collapse_large_unknown_blocks_in_DNA_sequence.py` using 
    # https://stackoverflow.com/a/5466478/8508004 you need to double up 
    # brackets on command line to use those search terms. I put some examples in
    # the demo notebook, `demo get_seq_following_seq_from_multiFASTA.ipynb`.
    
    # If more than one, go through parsed records and collect the record to act 
    # on.
    record_to_search  = None
    if single_record:
        record_to_search = records[0]
    else:
        for record in records:
            if record.id == record_id:
                record_to_search = record

    # Now search the record and extract the sequence after the pattern match.
    match_locations = get_start_n_ends_for_match_to_pattern(
                pat_obj,str(record_to_search.seq))
    if len(match_locations) > 1:
        sys.stderr.write("{} matches to the sequence found in the specified "
            "sequence. The sequence\nthat follows the one that occurs first "
            "has been returned.".format(len(match_locations)))
    if match_locations:
        start = match_locations[0][1]
        end = match_locations[0][1]+amount_to_get
        if len(record_to_search.seq[start:end]) < amount_to_get:
            sys.stderr.write("Note that the end of the sequence was "
            "encountered, and so less than the\nspecified {} residues were "
            "returned.".format(amount_to_get))
        return str(record_to_search.seq[start:end])
    else:
        sys.stderr.write("***NO MATCHES FOUND. RETURNING EMPTY STRING.***** "
            "   **** ERROR???**")
        return("")

    



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
    kwargs['filter_dashes'] = filter_dashes
    result = get_seq_following_seq_from_multiFASTA(
        sequence_file, record_id, seq_to_find, amount_to_get,**kwargs)
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
        'get_seq_following_seq_from_multiFASTA.py',
        description="get_seq_following_seq_from_multiFASTA.py \
        takes a sequence pattern string, a sequence file (FASTA-format), and a \
        record id and extracts a sequence of specified size following the \
        sequence pattern. Importantly, the regular expression search term \
        syntax is acceptable in the provided sequence pattern, although \
        anything dealing with case will be ignored. (The FASTA-formatted \
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
        contains one sequence, that sequence to be mined and what is provided \
        for this parameter will be ignored. In other words, if the sequence \
        file is not a multi-FASTA file, you don't need to determine the \
        identifier and can instead just enter `blahblah` or any other \
        nonsensical string in this spot.", metavar="RECORD_ID")
    parser.add_argument("pattern", help="Sequence or sequence pattern to use \
        to locate site after which to get the sequence. Regular expressions \
        are accepted here; however any information about case will be ignored \
        as the provided sequence pattern and sequence will both be converted \
        to lower case to check for a match.", metavar="PATTERN")
    parser.add_argument("amount_to_get", type=int, help="Number (integer) of \
        residues \
        to retrieve following the match to the sequence. The length of this \
        sequence is to be given in common terms, where the first item is \
        referenced as `1` ,and so a provided argument of \
        `1` would a single residue following the match would be returned.", 
        metavar="NUMBER_TO_GET")
    parser.add_argument('-ld', '--leave_dashes', help="Add this flag when \
        calling the script in \
        order to be able to use gaps (represented as dashes) in the pattern \
        required to match. I.E., for matching with an aligned FASTA file. \
        (***ATYPICAL.***)", action="store_true")





    #I would also like trigger help to display if no arguments provided because 
    # need at least one for url
    if len(sys.argv)==1:    #from http://stackoverflow.com/questions/4042452/display-help-message-with-python-argparse-when-script-is-called-without-any-argu
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    sequence_file = args.sequence_file
    record_id = args.record_id
    seq_to_find = args.pattern
    amount_to_get = args.amount_to_get
    if args.leave_dashes:
        filter_dashes = False
    else:
        filter_dashes = True



    main()

#*******************************************************************************
###-***********************END MAIN PORTION OF SCRIPT***********************-###
#*******************************************************************************
