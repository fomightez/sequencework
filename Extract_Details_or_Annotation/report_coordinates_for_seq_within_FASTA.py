#!/usr/bin/env python
# report_coordinates_for_seq_within_FASTA.py
__author__ = "Wayne Decatur" #fomightez on GitHub
__license__ = "MIT"
__version__ = "0.1.0"


# report_coordinates_for_seq_within_FASTA.py by 
# Wayne Decatur
# ver 0.1
#
#*******************************************************************************
# Verified compatible with both Python 2.7 and Python 3.6; written initially in 
# Python 3. 
#
#
# PURPOSE: Takes a sequence pattern string, a sequence file (FASTA-format), & a 
# record id, and reports the start and end coordinates of that sequence within
# the specified FASTA record.
# Because inherently the sequence to search can be provided as a pattern, I am
# saying 'pattern'; however, the impetus for this script was for when you have
# a known sequence, that will only occour once in that record and you need the 
# coordinates.
# (The FASTA-formatted sequence file in which to search is assumed by default 
# to be a multi-FASTA, i.e., multiple sequences in the provided file, although it 
# definitely doesn't have to be. In case it is only a single sequence, the
# record id becomes moot, see below.) 
# Two values, the start and the end will be returned, separated by a tab. 
# Redirect the output to a file if using command line version and a file is 
# needed.
#
# Importantly, the coordinates numbering is in 'common' terms where the position 
# numbered one corresponds to the first position. (In other words, the 
# coordinates returned are in the more conventional form & are not zero-indexed 
# even though behind-the-scenes Python/Biopython is indeed zero-indexed.)
#
#
# The provided sequence pattern will be matched regardless of case, as both the
# input sequence and pattern to search will be converted behind-the-scenes to 
# lowercase for the comparison. Beyond being insensitive of the case, REGULAR 
# EXPRESSION SEARCH TERM SYNTAX IS ACCEPTABLE in the provided sequence pattern.
#
# Note that if there is only one record in the specified sequence file, the 
# record id is moot and you can instead provide any string for that parameter 
# as it will be ignored. This makes the script more flexible in cases where 
# sequence files aren't complex as the user doesn't need to provide an actual 
# record id.
# 
# It is designed to handle/filter gaps ('dashes') in the provided sequence 
# patterns. The idea being that the known sequence ends may be manually 
# extracted from sequence alignments. This way the user is not wasting time 
# removing the gap indications / dashes from the collected text lines.
#
# Note that why some aspects of this script may seem redundant with my scripts 
# `find_sequence_element_occurrences_in_sequence.py` and `blast_to_df.py` that 
# include coordinates in the output tables, those were meant to mainly be used 
# when looking for a table of information about matches or exploring matches 
# that may occur multiples times in many sequences in the supplied file(s). This 
# script is targeted at the case where you know record id and expect only one 
# match to that record with that identifier.
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
# (based on `get_seq_following_seq_from_multiFASTA.py`)
# If you are Wayne, see 
# `Gathering XXXXXXX from aligned XXXXXXXXX genomes collected from the XXXXX genomes project.ipynb` 
# (and maybe just after that?) for impetus behind this script.
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
# python report_coordinates_for_seq_within_FASTA.py seq.fa record_id pattern
#-----------------------------------
#
# Issue `report_coordinates_for_seq_within_FASTA.py -h` for details.
# 
#
#
# To use this after importing/pasting or loading into a cell in a Jupyter 
# notebook, call the main function similar to below:
# following_seq = report_coordinates_for_seq_within_FASTA("seq.fa", "Skluv", "GAAATTTCCCCCAAAATGT")
#
# 
#
'''
CURRENT ACTUAL CODE FOR RUNNING/TESTING IN A NOTEBOOK WHEN IMPORTED/LOADED OR 
PASTED IN ANOTHER CELL:
following_seq = report_coordinates_for_seq_within_FASTA("seq.fa", "Skluv", "GAAATTTCCCCCAAAATGT")
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
    second is the end of the span of the pattern match.<--actually it is the 
    index of the first item not in the match, like slice notation.
    '''
    start_end_tuples = []
    for m in pattern_obj.finditer(a_string.lower()):
        start_end_tuples.append(m.span()) # from https://stackoverflow.com/a/250306/8508004
    return start_end_tuples

###--------------------------END OF HELPER FUNCTIONS-------------------------###
###--------------------------END OF HELPER FUNCTIONS-------------------------###



#*******************************************************************************
###------------------------'main' function of script--------------------------##

def report_coordinates_for_seq_within_FASTA(
    sequence_file, record_id, seq_to_find, filter_dashes = True):
    '''
    Main function of script.
    Takes a sequence string, a sequence file (FASTA-format), and a record id and 
    and reports the start and end coordinates of that sequence within the 
    specified FASTA record.
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
            "used to search for the provided sequence pattern\nand provide the "
            "coordinates spanned.\n\n".format(record_id))

    # Prepare the pattern to be searched for.
    if filter_dashes:
        seq_to_find = seq_to_find.replace("-","")
    # make a pattern object since it could be used often if I change things to 
    # search the file for multiple records instead of just one at a time
    pat_obj = re.compile(seq_to_find.lower())
    # It seemed it didn't (and wouldn't work) to use 
    # https://stackoverflow.com/a/12989308/8508004 to bring regular expression
    # search term from shell. Seems that similar to what I found at while 
    # working on `collapse_large_unknown_blocks_in_DNA_sequence.py` using 
    # https://stackoverflow.com/a/5466478/8508004 you need to double up 
    # brackets on command line to use those search terms. I put several 
    # examples in the demo notebook, 
    # `demo report_coordinates_for_seq_within_FASTA.ipynb`.
    
    # If more than one, go through parsed records and collect the record to act 
    # on.
    record_to_search  = None
    if single_record:
        record_to_search = records[0]
    else:
        for record in records:
            if record.id == record_id:
                record_to_search = record

    # Now search the record and return the coordinates of the match.
    match_locations = get_start_n_ends_for_match_to_pattern(
                pat_obj,str(record_to_search.seq))
    if len(match_locations) > 1:
        sys.stderr.write("{} matches to the sequence found in the specified "
            "sequence. The coordinates\nof the match encountered first "
            "have been returned.".format(len(match_locations)))
    if match_locations:
        start = match_locations[0][0]+1 #`+1` hides zero-index behind-the-scenes
        end = match_locations[0][1] #`+1`is not necessary here to hide 
        # zero-index that behind-the-scenes b/c second value `span()` seems to
        # return seems to be the 'stop' position index like Python slices where
        # it is up to but not including the `stop` index despite it just saying
        # start and end index are returned at
        # https://stackoverflow.com/a/250306/8508004 . You can see it in the
        # example there, the 'm' in 'message' is at the fourth position in the
        # string analyzed but the 'e' at the end is at the tenth position and 
        # not the 11th, but `span` returns (4,11). So `+1` already built in.
        return "{}\t{}".format(start, end)
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
    result = report_coordinates_for_seq_within_FASTA(
        sequence_file, record_id, seq_to_find,**kwargs)
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
        'report_coordinates_for_seq_within_FASTA.py',
        description="report_coordinates_for_seq_within_FASTA.py \
        takes a sequence pattern string, a sequence file (FASTA-format), and a \
        record id, and reports the start and end coordinates of that sequence \
        within the specified FASTA record. Importantly, the coordinates \
        numbering is in 'common' terms where the position numbered one \
        corresponds to the first position.\
        **** Script by Wayne Decatur   \
        (fomightez @ github) ***")

    parser.add_argument("sequence_file", help="Name of sequence file to \
        use as input. Must be FASTA format. Can be a \
        multi-FASTA file, i.e., multiple sequences in FASTA format in one \
        file.", metavar="SEQUENCE_FILE")
    parser.add_argument("record_id", help="Specific identifier of sequence \
        entry in sequence file to search. If the provided sequence file only \
        contains one sequence, that sequence will be searchedand what is \
        provided for this parameter will be ignored. In other words, if the \
        sequence file is not a multi-FASTA file, you don't need to determine \
        the identifier and can instead just enter `blahblah` or any other \
        nonsensical string in this spot.", metavar="RECORD_ID")
    parser.add_argument("pattern", help="Sequence or sequence pattern to use \
        to determine the coordinates it spans. Regular expressions \
        are accepted here; however any information about case will be ignored \
        as the provided sequence pattern and sequence will both be converted \
        to lower case to check for a match.", metavar="PATTERN")






    #I would also like trigger help to display if no arguments provided because 
    # need at least one for url
    if len(sys.argv)==1:    #from http://stackoverflow.com/questions/4042452/display-help-message-with-python-argparse-when-script-is-called-without-any-argu
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    sequence_file = args.sequence_file
    record_id = args.record_id
    seq_to_find = args.pattern



    main()

#*******************************************************************************
###-***********************END MAIN PORTION OF SCRIPT***********************-###
#*******************************************************************************
