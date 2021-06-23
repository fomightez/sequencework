#!/usr/bin/env python
# get_seq_from_multiFASTA_with_match_in_description.py
__author__ = "Wayne Decatur" #fomightez on GitHub
__license__ = "MIT"
__version__ = "0.1.0"


# get_seq_from_multiFASTA_with_match_in_description.py by 
# Wayne Decatur
# ver 0.1
#
#*******************************************************************************
# Verified compatible with both Python 2.7 and Python 3.6; written initially in 
# Python 3. 
#
#
# PURPOSE: Takes any sequences in FASTA format and gets the first sequence 
# with a description line containing a match to provided text string. For 
# example, if provided a multi-sequence FASTA file and a gene identifier, such
# as `YDL140C`, it will pull out the first sequence matching that anywhere in
# the description line.
# Defaults to ignoring case.
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
# pyfaidx
#
#
# VERSION HISTORY:
# v.0.1. basic working version.
#
#
# To do:
# - use what I did in `get_seq_following_seq_from_FASTA.py` to add ability
# to use regex in provided search text
# - Would I have need for multiple matches?
#
#
#
#
# TO RUN:
# Examples,
# Enter on the command line of your terminal, the line
#-----------------------------------
# python get_seq_from_multiFASTA_with_match_in_description.py sequence.fa text_to_match
#-----------------------------------
#
# Issue `get_seq_from_multiFASTA_with_match_in_description.py -h` for details.
# 
#
#
# To use this after importing/pasting or loading into a cell in a Jupyter 
# notebook, specify at least the sequence file (or results as a string) in the 
# call to the main function similar to below:
# get_seq_from_multiFASTA_with_match_in_description("sequence.fa", "YDL140C")
# 
#
# 
#
'''
CURRENT ACTUAL CODE FOR RUNNING/TESTING IN A NOTEBOOK WHEN IMPORTED/LOADED OR 
PASTED IN ANOTHER CELL:
get_seq_from_multiFASTA_with_match_in_description("sequence.fa", "YDL140C")
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


suffix_for_saving = "_" #to be used for naming the output automatically
# when running script from command line to act on an input file

case_sensitive = False # So defaults to insenstive when called from command line

#
#*******************************************************************************
#**********************END USER ADJUSTABLE VARIABLES****************************






















#*******************************************************************************
#*******************************************************************************
###DO NOT EDIT BELOW HERE - ENTER VALUES ABOVE###

import sys
import os
#from Bio import SeqIO
#from Bio.Seq import Seq 
#from Bio.SeqRecord import SeqRecord 
#from Bio.Alphabet import generic_dna
from pyfaidx import Fasta



###---------------------------HELPER FUNCTIONS-------------------------------###


def generate_output_file_name(file_name,text_to_match, suffix_for_saving):
    '''
    Takes a file name as an argument and returns string for the name of the
    output file. The generated name is based on the text to match.


    Specific example
    =================
    Calling function with
        ("sequence.fa", "YDL140C", "_")
    returns
        "seq_YDL140C.fa"
    '''
    return "seq" + suffix_for_saving + text_to_match + ".fa"

def generate_seq_chunks(seq_string, chunk_size = 70):
    '''
    This takes a sequence as a string and breaks it up into list of strings of
    set length with graceful handling of the last line that will most likely
    not be full length. 
    The list of sequence strings gets returned

    `chunk_size = 70` sets residues per line to have in FASTA; note this 
    #chunking to multiple lines is the opposite of what PatMatch's 
    # `unjustify_fasta.pl` does. 

    Note `chunk_size` defaults to 70 here but optionally a different one
    can be provided.

    I believe my chunking code is based on 
    https://stackoverflow.com/a/13673133/8508004 or 
    https://stackoverflow.com/a/9475354/8508004 , see my gist 
    https://gist.github.com/fomightez/ef7583919dde51f3569731ca1c5247ba for some 
    notes on it and more related.
    There's a related code in my script 
    `similarities_in_proteinprotein_interactions.py` inside my repo 
    structurework/pdbsum-utilities/  and in 
    `Searching for homologs among deduced proteins from a genome using BLAST and Python.ipynb`
    & `Searching for coding sequences in genomes using BLAST and Python.ipynb`
    & `notebooks/GSD/GSD Rpb1_orthologs_in_1011_genomes.ipynb` inside my repo 
    blast-binder.
    '''
    return [seq_string[i:i+chunk_size] for i in range(
        0, len(seq_string),chunk_size)]


    prot_seq_fa = ">" + prot_descr + "\n"+ "\n".join(prot_seq_chunks)


###--------------------------END OF HELPER FUNCTIONS-------------------------###
###--------------------------END OF HELPER FUNCTIONS-------------------------###



#*******************************************************************************
###------------------------'main' function of script--------------------------##

def get_seq_from_multiFASTA_with_match_in_description(sequence, text_to_match,
    case_sensitive = False,
    return_record_as_string = False,
    suffix_for_saving = suffix_for_saving):
    '''
    Main function of script. 
    Takes any sequences in FASTA format and gets the first sequence with a 
    description line containing a match to provided text string. For
    example, if provided a multi-sequence FASTA file and a gene identifier, such 
    as `YDL140C`, it will pull out the first sequence matching that.

    Has `return_record_as_string` option so it can be imported and used in
    IPython or Jupyter notebook and the result used directly without need to
    access generated file. Plus, this then can be used for more flexible output 
    FASTA file naming or making combined multi-sequene FASTA files, see the demo
    in my cl_sq_demo-binder.
    '''

    # get fasta records using pyfaidx
    records = Fasta(sequence)
    #for record in SeqIO.parse(sequence, "fasta"):
    #   records.append(record)

    # go through records and collect necessary sequence
    seq_fa = "NOT_ANY_FOUND"
    for record in records:
        if case_sensitive:
            if text_to_match in record.long_name:
                seq_fa = ">" + record.long_name + "\n"+"\n".join(
                    generate_seq_chunks(str(record)))
                break
        else:
            if text_to_match.lower() in record.long_name.lower():
                seq_fa = ">" + record.long_name + "\n"+"\n".join(
                    generate_seq_chunks(str(record)))
                break
    if seq_fa == "NOT_ANY_FOUND":
        sys.stderr.write("**ERROR:No match to provided text found in "
            "description line "
            "for ANY sequence record.  ***ERROR*** \nEXITING.\n")
        sys.exit(1)





    # Save the extracted FASTA record
    output_file_name = generate_output_file_name(
        sequence,text_to_match,suffix_for_saving)
    with open(output_file_name, 'w') as output:
        output.write(seq_fa)
    # add new line to the end of the saved file so when concatanating sequence 
    # it will start on next line
    with open(output_file_name, 'a') as output: 
        output.write('\n')
    # Feedback
    sys.stderr.write("\n\n*****************DONE**************************\n")
    sys.stderr.write("Extracted sequence saved in FASTA format as '{}'.".format(
            output_file_name))
    sys.stderr.write("\n*****************DONE**************************\n")


    # Return the sequence record as a string if requested
    if return_record_as_string:
        return seq_fa



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
    kwargs['case_sensitive'] = case_sensitive
    get_seq_from_multiFASTA_with_match_in_description(
        sequence,text_to_match,**kwargs)
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
        'get_seq_from_multiFASTA_with_match_in_description.py',
        description="get_seq_from_multiFASTA_with_match_in_description.py \
        takes any sequences in FASTA format and gets the first sequence  with \
        a description line containing a match to provided text string. For \
        example, if provided a multi-sequence FASTA file and a gene \
        identifier, such as `YDL140C`, it will pull out the first sequence \
        matching that anywhere in description line. Defaults to ignoring case. \
        **** Script by Wayne Decatur   \
        (fomightez @ github) ***")

    parser.add_argument("sequence_file", 
        help="Name of sequence file to search.", metavar="SEQUENCE_FILE")

    parser.add_argument("text_to_match", 
        help="Text to match.", metavar="TEXT_TO_MATCH")

    parser.add_argument('-cs', '--case_sensitive', action='store_true', 
    default= case_sensitive,
    help="Add this flag if you want to force matching to be case-sensitive.")
    
    #parser.add_argument('-os', '--output_suffix', action='store', type=str, 
    #default= suffix_for_saving, help="OPTIONAL: Set a suffix for including in \
    #file name of output. \
    #If none provided, '{}' will be used.".format(suffix_for_saving))




    #I would also like trigger help to display if no arguments provided because 
    # need at least one for url
    if len(sys.argv)==1:    #from http://stackoverflow.com/questions/4042452/display-help-message-with-python-argparse-when-script-is-called-without-any-argu
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    sequence = args.sequence_file
    text_to_match = args.text_to_match
    case_sensitive = args.case_sensitive
    #suffix_for_saving = args.output_suffix


    main()

#*******************************************************************************
###-***********************END MAIN PORTION OF SCRIPT***********************-###
#*******************************************************************************
