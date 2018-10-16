#!/usr/bin/env python
# delete_seq_following_pattern_within_multiFASTA.py
__author__ = "Wayne Decatur" #fomightez on GitHub
__license__ = "MIT"
__version__ = "0.1.0"


# delete_seq_following_pattern_within_multiFASTA.py by 
# Wayne Decatur
# ver 0.1
#
#*******************************************************************************
# Verified compatible with both Python 2.7 and Python 3.6; written initially in 
# Python 3. 
#
#
# PURPOSE: Takes a sequence pattern string, a sequence file (FASTA-format), & a 
# record id, and deletes any sequence following the sequence pattern. In other
# words it trims the specified sequence, to make the first match to the pattern
# the new end. (The FASTA-formatted sequence file is assumed by default to be a 
# multi-FASTA, i.e., multiple sequences in the provided file, although it 
# definitely doesn't have to be. In case it is only a single sequence, the
# record id becomes moot, see below.) 
# Nothing will be returned but a copy of the FASTA sequence file with the 
# truncated sequence will be produced.
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
# This script is meant to be used after you have performed a large alignment, 
# say of an entire chromosome, in order to have individual occurrences of 
# related segments fall linearly with where they match up along the span of the 
# sequence. Often due to large (seeming-to-be) arbitratrily-sized blocks of 
# repeated unknown nucleotides (which are often good to 'collapse', see 
# `collapse_large_unknown_blocks_in_DNA_sequence.py`) the 'ends' of regions 
# often fail to get extracted exactly right and you can end up with some 
# sequences that trail on for longer than they should.
#
# It is designed to handle/filter gaps ('dashes') in the provided sequence 
# patterns. The idea being that the known sequence ends may be manually 
# extracted from sequence alignments. This way the user is not wasting time 
# removing the gap indications / dashes from the collected text lines. The 
# default handling of removing the gaps to ignore them can be overriden. 
# The idea is that maybe you'll have a multiple sequence alignment 
# file saved as FASTA with dashes, i.e., aligned FASTA file format and may want 
# to use this script. 
#
#
# (based on `get_seq_following_seq_from_multiFASTA.py`)
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
# python delete_seq_following_pattern_within_multiFASTA.py seq.fa record_id pattern > seq_clipped.fa
#-----------------------------------
#
# Issue `delete_seq_following_pattern_within_multiFASTA.py -h` for details.
# 
#
#
# To use this after importing/pasting or loading into a cell in a Jupyter 
# notebook, specify at least the results file (or results as a string) in the 
# call to the main function similar to below:
# delete_seq_following_pattern_within_multiFASTA("seq.fa", "Skluv", "GAAATTTCCCCCAAAATGT")
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
delete_seq_following_pattern_within_multiFASTA("seq.fa", "Skluv", "GAAATTTCCCCCAAAATGT" , 200)
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

suffix_for_saving = "_clipped" #to be used for naming the output automatically
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
        ("sequence.fa", "_clipped")
    returns
        "sequence_clipped.fa"
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
    for m in pattern_obj.finditer(a_string.lower()):
        start_end_tuples.append(m.span()) # from https://stackoverflow.com/a/250306/8508004
    return start_end_tuples

###--------------------------END OF HELPER FUNCTIONS-------------------------###
###--------------------------END OF HELPER FUNCTIONS-------------------------###



#*******************************************************************************
###------------------------'main' function of script--------------------------##

def delete_seq_following_pattern_within_multiFASTA(
    sequence_file, record_id, seq_to_find, filter_dashes = True, 
    suffix_for_saving = suffix_for_saving):
    '''
    Main function of script.
    Takes a sequence string, a sequence file (FASTA-format), and a record id and 
    deletes all the seqeunce following the match to the sequence pattern string
    in the specified record. 
    The FASTA-formatted sequence file is assumed by default to be a 
    multi-FASTA, i.e., multiple sequences in the provided file, although it 
    definitely doesn't have to be. In case it is only a single sequence, the 
    record id becomes moot and users can provide anything for this parameter. 

    A modified sequence file will be made.
    '''
    # make a records list because will most likely need to modify a record. 
    # Plus, knowing how many records there are will affect how records handled
    # because don't want record_id to matter in such a case & 
    # be sure to give feedback of that sort.
    records = []  #intialize a list of records for the FASTA file
    for record in SeqIO.parse(sequence_file, "fasta"):
        records.append(record)
    single_record = False
    if len(records) == 1:
        record_id = records[0].id
        single_record = True
        # feedback
        sys.stderr.write("Single sequence with id of {} provided in the "
            "sequence file. It will be "
            "used to search for the provided sequence pattern and delete the "
            "residues after it.".format(record_id))

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
    # examples of using regular expression from command line and supplying to a 
    # function in the demo notebook, 
    # `demo get_seq_following_seq_from_multiFASTA.ipynb`.
    
    # If more than one, go through parsed records and collect the record to act 
    # on.
    record_to_search  = None
    if single_record:
        record_to_search = records[0]
        index_of_record = 0
        record_description = records[0].description
    else:
        for indx,record in enumerate(records):
            if record.id == record_id:
                record_to_search = record
                index_of_record = indx
                record_description = record.description

    # Now search the record and delete the sequence after the pattern match.
    match_locations = get_start_n_ends_for_match_to_pattern(
                pat_obj,str(record_to_search.seq))
    if len(match_locations) > 1:
        sys.stderr.write("{} matches to the sequence found in the specified "
            "sequence. The sequence\nthat follows the match encountered first "
            "has been deleted.".format(len(match_locations)))
    if match_locations:
        end = match_locations[0][1]
        modf_seq = str(record_to_search.seq[:end])
        # replace the appropriate record with the modified sequence, tagging
        # the description line
        records[index_of_record] = SeqRecord(
            Seq(modf_seq, generic_dna), 
            id=record_id, description=record_description+" CLIPPED")#based
        # on https://www.biostars.org/p/48797/ and `.ungap()` method, see
        # https://github.com/biopython/biopython/issues/1511 , and `description`
        # from what I've seen for `id` plus https://biopython.org/wiki/SeqIO
        #print (records[indx]) # ONLY FOR DEBUGGING
    else:
        modf_seq = None
        


    # Replace the FASTA file with the modified records
    output_file_name = generate_output_file_name(sequence_file,suffix_for_saving)
    SeqIO.write(records,output_file_name, "fasta");
    # Feedback
    sys.stderr.write("\n\n*****************DONE**************************\n")
    if modf_seq:
        sys.stderr.write("Sequence after the match to the provided pattern "
            "\ndeleted from '{}', and saved within a modified version \nof "
            "'{}' as the output file '{}'.".format(
            record_id, sequence_file, output_file_name))
    else:
        sys.stderr.write("***NO MATCHES FOUND. NO CHANGES MADE.***** "
            "   **** ERROR?!?!?**")
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
    kwargs['filter_dashes'] = filter_dashes
    kwargs['suffix_for_saving'] = suffix_for_saving
    delete_seq_following_pattern_within_multiFASTA(
        sequence_file, record_id, seq_to_find, **kwargs)
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
        'delete_seq_following_pattern_within_multiFASTA.py',
        description="delete_seq_following_pattern_within_multiFASTA.py \
        takes a sequence pattern string, a sequence file (FASTA-format), and a \
        record id, and deletes any sequence following the sequence pattern. In \
        other words it trims the specified sequence, to make the first match \
        to the pattern the new end. (The FASTA-formatted sequence file is \
        assumed by default to be a multi-FASTA, i.e., multiple sequences in \
        the provided file, although it definitely doesn't have to be. In case \
        it is only a single sequence, the record id becomes moot, see below.) \
        Nothing will be returned; however a copy of the FASTA sequence file \
        with the truncated sequence specified will be produced. \
        **** Script by Wayne Decatur   \
        (fomightez @ github) ***")

    parser.add_argument("sequence_file", help="Name of sequence file to \
        use as input. Must be FASTA format. Can be a \
        multi-FASTA file, i.e., multiple sequences in FASTA format in one \
        file.", metavar="SEQUENCE_FILE")
    parser.add_argument("record_id", help="Specific identifier of sequence \
        entry in sequence file to modify. If the provided sequence file only \
        contains one sequence, then that sequence will be altered, and \
        whatever is provided \
        for this parameter will be ignored. In other words, if the sequence \
        file is not a multi-FASTA file, you don't need to determine the \
        identifier and can instead just enter `blahblah` or any other \
        nonsensical string in this spot.", metavar="RECORD_ID")
    parser.add_argument("pattern", help="Sequence or sequence pattern to use \
        to locate site after which to delete in the specified sequence. \
        Regular expressions \
        are accepted here; however any information about case will be ignored \
        as the provided sequence pattern and sequence will both be converted \
        to lower case to check for a match.", metavar="PATTERN")
    parser.add_argument('-ld', '--leave_dashes', help="Add this flag when \
        calling the script in \
        order to be able to use gaps (represented as dashes) in the pattern \
        required to match. I.E., for matching with an aligned FASTA file. \
        (***ATYPICAL.***)", action="store_true")

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
    sequence_file = args.sequence_file
    record_id = args.record_id
    seq_to_find = args.pattern
    suffix_for_saving = args.output_suffix
    if args.leave_dashes:
        filter_dashes = False
    else:
        filter_dashes = True



    main()

#*******************************************************************************
###-***********************END MAIN PORTION OF SCRIPT***********************-###
#*******************************************************************************
