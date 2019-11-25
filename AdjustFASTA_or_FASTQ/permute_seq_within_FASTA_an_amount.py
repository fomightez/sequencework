#!/usr/bin/env python
# permute_seq_within_FASTA_an_amount.py
__author__ = "Wayne Decatur" #fomightez on GitHub
__license__ = "MIT"
__version__ = "0.1.0"


# permute_seq_within_FASTA_an_amount.py by 
# Wayne Decatur
# ver 0.1
#
#*******************************************************************************
# Verified compatible with both Python 2.7 and Python 3.6; written initially in 
# Python 3. 
#
#
# PURPOSE: Takes a sequence file (FASTA-format) & a record id (unless single 
# sequence in file) of a circular sequence such as a plasmid or mitochondrial 
# genome, and moves the specified distance in bps the start/end (breakppoint) of
# that sequence within the specified FASTA record. In other words, the so the 
# position that is the specified amount from the original 'start' becomes the 
# new 'start' position.
# (The FASTA-formatted sequence file in which to act is assumed by default 
# to be a multi-FASTA, i.e., multiple sequences in the provided file, although it 
# definitely doesn't have to be. In case it is only a single sequence, the
# record id becomes moot, see below.) 
# Nothing is returned when using this on the command line. A file will be saved
# of the permuted sequence. In the case of importing and using the main funcion,
# by default the name of the saved file is returned. Optionally, when calling 
# the function, returning the name of the permuted file can be disabled.
#
#
# Note that if there is only one record in the specified sequence file, the 
# record id is moot and you can instead provide any string for that parameter 
# as it will be ignored. This makes the script more flexible in cases where 
# sequence files aren't complex as the user doesn't need to provide an actual 
# record id.
#
# A good reference for the use of the term 'breakpoint' for the arbitrary and 
# artificial location of the 'start'/'end' of circular sequence is at 
# https://software.broadinstitute.org/gatk/blog?id=23598 .
#
# Typical use cases for this script include:
# 1) Permuting a circular sequence to match the published or 'standard' sequence 
# for ease in comparison. Or vice versa. For example if you have 
# experimentally-determined sequence that spans the artifical start/end 
# breakpoint and you want to compare.
# 2) For a circular sequence, if you don't find an expected sequence, you may 
# wish to permute to see if possible the sequence happens to fall across the
# artificial start/end (end/start) breakpoint.
#
# Note that this processing is related to processing I do where first find match
# to 'start' sequence. See `Permuting PacBio yeast mito chrs to match SGD.ipynb`
# and related notebooks. The difference is that here only adjusting a specified
# amount. As of writing this script, I haven't scripted the adjustment based on
# adjusting to a specific location, but this script should serve as a basis
# where would just need to determing the 'new' start location.
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
# `Making_dataframes_of_UNPERMUTED_details_of_XXXXXXs_identified_by_XXXXXXXXXXXX_approach_for_PacBio_and_1011.ipynb` 
# (for impetus behind this script.
# Backbone for script was `report_coordinates_for_seq_within_FASTA.py` only b/c
# I had happened to edit that script recently and knew it allowed for handling
# multi-fasta and ids where you may often only have one sequence and the id 
# would be moot. 
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
# python permute_seq_within_FASTA_an_amount.py seq.fa record_id distance
#-----------------------------------
#
# Issue `permute_seq_within_FASTA_an_amount.py -h` for details.
# 
#
#
# To use this after importing/pasting or loading into a cell in a Jupyter 
# notebook, call the main function similar to below:
# permuted_seq = permute_seq_within_FASTA_an_amount("seq.fa", "Skluv", 1001)
#
# 
#
'''
CURRENT ACTUAL CODE FOR RUNNING/TESTING IN A NOTEBOOK WHEN IMPORTED/LOADED OR 
PASTED IN ANOTHER CELL:
permuted_seq = permute_seq_within_FASTA_an_amount("seq.fa", "Skluv", 1001)
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
suffix_for_output_file = "_permuted.fa"


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


def generate_output_file_name(file_name, suffix_for_output_file):
    '''
    Takes a file name as an argument and returns string for the name of the
    output file. The generated name is tagged to indicate difference.


    Specific example
    =================
    Calling function with
        ("sequences.fa", "_permuted.fa")
    returns
        "sequences_permuted.fa"
    '''
    main_part_of_name, file_extension = os.path.splitext(
        file_name) #from http://stackoverflow.com/questions/541390/extracting-extension-from-filename-in-python
    if '.' in file_name:  #I don't know if this is needed with the os.path.splitext method but I had it before so left it
        return main_part_of_name + suffix_for_output_file
    else:
        return file_name + text_to_add_to_altered_file + suffix_for_output_file

###--------------------------END OF HELPER FUNCTIONS-------------------------###
###--------------------------END OF HELPER FUNCTIONS-------------------------###



#*******************************************************************************
###------------------------'main' function of script--------------------------##

def permute_seq_within_FASTA_an_amount(
    sequence_file, record_id, distance_to_permute = 1001,
    return_file_produced=True):
    '''
    Main function of script.
    Takes a sequence file (FASTA-format), and a record id and 
    and permutes that specified circular sequence so the position that is the
    specified amount from the original 'start' becomes the new 'start' position.
    The FASTA-formatted sequence file is assumed by default to be a 
    multi-FASTA, i.e., multiple sequences in the provided file, although it 
    definitely doesn't have to be. In case it is only a single sequence, the 
    record id becomes moot and users can provide anything for this parameter. 

    When calling from the command line, nothing will be returned and the name of 
    the sequence file made with the permuted sequence will be noted in stderr.
    When calling the function the setting for `return_file_produced` will 
    default to true and the name of the generated output file will be returned.
    This can be set to false to not return that.
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
        record_description= records[0].description
        single_record = True
        # feedback
        sys.stderr.write("Single sequence with id of '{}' provided in the "
            "sequence file.\nIt will be permuted.\n\n".format(record_id))

    # If more than one, go through parsed records and collect the record to act 
    # on.
    record_to_permute  = None
    if single_record:
        record_to_permute = records[0]
    else:
        for record in records:
            if record.id == record_id:
                record_to_permute = record
                record_description= record.description

    # Now permute the record. (Below permuting step based on `Permuting PacBio yeast mito chrs to match SGD.ipynb` and related notebooks.)
    # fix, based on where it says "i.e. shift the starting point on this plasmid," @ 
    #http://biopython.org/DIST/docs/api/Bio.SeqRecord.SeqRecord-class.html .
    # Specify the amount currently at start to move to end with 'left'.
    left = record_to_permute.seq[:distance_to_permute-1] # use one less than 
    # position desired to be new start position in order to shift to zero-index
    # behind the scenes.  And want users to refer to the standard way of the
    # firstbasepair being numbered one and oneward. For example if for the sake of 
    # simplicity  a 'sequence' matchiing the word 'serenade' was being permuted
    # so the third letter was becoming the first, you'd need to move first two 
    # letters to end and slice for that would be `"serenade"[:3-1]`.
    # Next, specify the part to become the new start with 'right' because it 
    # currently includes what will be the new start & what is currently to 
    # 'right' of it.
    right = record_to_permute.seq[distance_to_permute-1:] # Use `-1` because 
    # want user to refer to position that becomes new start in 'typical' way 
    # where first basepair is numbered as one. So this is just really shifting 
    # to zero-indexing behind-the-scenes.
    permuted_seq = right + left #where 'right' and 'left' refer to what was 
    # relative the position that is now the new start BEFORE permuted
    permuted_seq_record = SeqRecord(permuted_seq,  
        id= record_id, description=record_description.split(record_id,1)[1].strip()+"[PERMUTED]")

    # write result after altered (Writing based on sequence file writing of `
    # fix_lsu_rRNA_annotation_in_gff_resulting_from_mfannot.py`.)
    output_file_name = generate_output_file_name(
        sequence_file, suffix_for_output_file)
    SeqIO.write(
        permuted_seq_record, output_file_name, "fasta");
    if return_file_produced:
        return output_file_name
    else:
        sys.stderr.write("The file {} with the sequence start permuted to "
            "position {} has been generated.".format(
            output_file_name,distance_to_permute))

    



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
    kwargs['return_file_produced'] = False # when calling on command line, nothing
    # will be returned and sequence file generated will be noted in stderr 
    result = permute_seq_within_FASTA_an_amount(
        sequence_file, record_id, distance_to_permute,**kwargs)
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
        'permute_seq_within_FASTA_an_amount.py',
        description="permute_seq_within_FASTA_an_amount.py \
        takes a sequence file (FASTA-format), and a \
        record id, and permutes that circular sequence so the position that is \
        the specified amount from the original 'start' becomes the new 'start' \
        position..\
        **** Script by Wayne Decatur   \
        (fomightez @ github) ***")

    parser.add_argument("sequence_file", help="Name of sequence file to \
        use as input. Must be FASTA format. Can be a \
        multi-FASTA file, i.e., multiple sequences in FASTA format in one \
        file.", metavar="SEQUENCE_FILE")
    parser.add_argument("record_id", help="Specific identifier of sequence \
        entry in sequence file to permute. If the provided sequence file only \
        contains one sequence, that sequence will be permuted and what is \
        provided for this parameter will be ignored. In other words, if the \
        sequence file is not a multi-FASTA file, you don't need to determine \
        the identifier and can instead just enter `blahblah` or any other \
        nonsensical string in this spot.", metavar="RECORD_ID")
    parser.add_argument('distance', action="store", type=int, help="Amount of \
        base pairs to permute the sequence; the new start will be this amount \
        after the original start point",metavar="DISTANCE")
    






    #I would also like trigger help to display if no arguments provided because 
    # need at least one for url
    if len(sys.argv)==1:    #from http://stackoverflow.com/questions/4042452/display-help-message-with-python-argparse-when-script-is-called-without-any-argu
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    sequence_file = args.sequence_file
    record_id = args.record_id
    distance_to_permute = args.distance



    main()

#*******************************************************************************
###-***********************END MAIN PORTION OF SCRIPT***********************-###
#*******************************************************************************