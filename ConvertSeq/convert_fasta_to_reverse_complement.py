#!/usr/bin/env python
# convert_fasta_to_reverse_complement.py
__author__ = "Wayne Decatur" #fomightez on GitHub
__license__ = "MIT"
__version__ = "0.1.0"


# convert_fasta_to_reverse_complement.py by 
# Wayne Decatur
# ver 0.1
#
#*******************************************************************************
# Verified compatible with both Python 2.7 and Python 3.6; written initially in 
# Python 3. 
#
#
# PURPOSE: Converts a provided sequence (or sequences) in FASTA format to the 
# reverse complement. The sequence remains in FASTA format.
# Also works if the provided FASTA file is a multi-FASTA, i.e., contains 
# multiple sequences in FASTA format in the one file. All sequences in the file
# will be converted to the reverse complement sequence.
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
# python convert_fasta_to_reverse_complement.py sequence.fa
#-----------------------------------
#
# Issue `convert_fasta_to_reverse_complement.py -h` for details.
# 
#
#
# To use this after importing/pasting or loading into a cell in a Jupyter 
# notebook, specify at least the sequence file (or results as a string) in the 
# call to the main function similar to below:
# convert_fasta_to_reverse_complement("sequence.fa")
#-or-
# To specify OPTIONAL suffix to add to the names of generated files, specify 
# `suffix_for_saving` like on next line when calling main function:
# extract_regions_from_clustal_alignment("sequence.fa",suffix_for_saving="rev_compl")
# 
#
# 
#
'''
CURRENT ACTUAL CODE FOR RUNNING/TESTING IN A NOTEBOOK WHEN IMPORTED/LOADED OR 
PASTED IN ANOTHER CELL:
convert_fasta_to_reverse_complement("sequence.fa")
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

suffix_for_saving = "_rc" #to be used for naming the output automatically
# when running script from command line to act on an input file

#
#*******************************************************************************
#**********************END USER ADJUSTABLE VARIABLES****************************






















#*******************************************************************************
#*******************************************************************************
###DO NOT EDIT BELOW HERE - ENTER VALUES ABOVE###

import sys
import os
from Bio import SeqIO
from Bio.Seq import Seq 




###---------------------------HELPER FUNCTIONS---------------------------------###


def generate_output_file_name(file_name,suffix_for_saving):
    '''
    Takes a file name as an argument and returns string for the name of the
    output file. The generated name is based on the original file
    name.


    Specific example
    =================
    Calling function with
        ("sequence.fa", "_rc")
    returns
        "sequence_rc.fa"
    '''
    main_part_of_name, file_extension = os.path.splitext(
        file_name) #from 
    #http://stackoverflow.com/questions/541390/extracting-extension-from-filename-in-python
    if '.' in file_name:  #I don't know if this is needed with the os.path.splitext method but I had it before so left it
        return main_part_of_name + suffix_for_saving  + file_extension
    else:
        return file_name + suffix_for_saving + ".fa"

###--------------------------END OF HELPER FUNCTIONS---------------------------###
###--------------------------END OF HELPER FUNCTIONS---------------------------###



#*******************************************************************************
###------------------------'main' function of script---------------------------##

def convert_fasta_to_reverse_complement(sequence,
    suffix_for_saving = suffix_for_saving):
    '''
    Main function of script. 
    Converts a provided sequence in FASTA format to the reverse complement.
    The sequence remains in FASTA format.
    '''
    records = []
    records_rc = []
    for record in SeqIO.parse(sequence, "fasta"):
        #records.append(record)
        records_rc.append(record.reverse_complement(id=True,description=True)) #
        # see https://biopython.org/DIST/docs/api/Bio.SeqRecord.SeqRecord-class.html#reverse_complement 
        # for arguments. They were necessary to pass information in FASTA 
        # description line from original into new. Otherwise got 
        # `><unknown id> <unknown description> ` on first line.
        
    # Save reverse complement FASTA
    output_file_name = generate_output_file_name(sequence,suffix_for_saving)
    SeqIO.write(records_rc,output_file_name, "fasta");
    # based on https://www.biostars.org/p/48797/#48801 "If that's the case then 
    # note SeqIO.write() can take a list or a generator of SeqRecords so you 
    # should pass one of those"
    # Feedback
    sys.stderr.write("\n\n*****************DONE**************************\n"
        "Sequences in FASTA file '{}' converted to reverse "
        "complement\nand saved as '{}'.\n"
        "*****************DONE**************************\n".format(
        sequence,output_file_name))



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
    kwargs['suffix_for_saving'] = suffix_for_saving
    convert_fasta_to_reverse_complement(sequence,**kwargs)
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
        'convert_fasta_to_reverse_complement.py',
        description="convert_fasta_to_reverse_complement.py \
        Converts provided sequence(s) in FASTA format to the reverse \
        complement sequence.\
        **** Script by Wayne Decatur   \
        (fomightez @ github) ***")

    parser.add_argument("sequence_file", help="Name of sequence file to \
        use as input for converting to reverse complement. Can be a \
        multi-FASTA file, i.e., multiple sequences in FASTA format in one file.\
        All included will be converted.", metavar="SEQUENCE_FILE")

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
    suffix_for_saving = args.output_suffix


    main()

#*******************************************************************************
###-***********************END MAIN PORTION OF SCRIPT***********************-###
#*******************************************************************************
