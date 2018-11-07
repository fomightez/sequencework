#!/usr/bin/env python
# reverse_complement_of_clustal_alignment.py
__author__ = "Wayne Decatur" #fomightez on GitHub
__license__ = "MIT"
__version__ = "0.1.0"


# reverse_complement_of_clustal_alignment.py by 
# Wayne Decatur
# ver 0.1
#
#*******************************************************************************
# Verified compatible with both Python 2.7 and Python 3.6; written initially in 
# Python 3. 
#
#
# PURPOSE: Takes a text document of an alignment in CLUSTAL format and generates
# that same alignment but with composing sequences as reverse complements.
# 
# Intended to be used after `extract_regions_from_clustal_alignment.py` for 
# sequences on the other (a.k.a. reverse) strand because some might be on other 
# strand when extract from entire chromosome and sometimes still it would make 
# more sense to view the extracted sequence in the other direction. For 
# example, if typically used to looking at the regions in one direction, or that
# make 'biological sense', etc.
# 
# Because alignment files contain no 'description' information, the script has 
# an option to provide the file of the ungapped, FASTA-formatted sequences 
# that are output as `*_extracted_ungapped.fa` files from the 
# `extract_regions_from_clustal_alignment.py` in order
# to have the description lines added to the ungpapped, reverse complement FASTA
# -formatted sequences this script includes as output.
#
#
# 
#
# Written to run from command line or imported into/pasted/loaded inside a 
# Jupyter notebook cell. When doing in Jupyter (IPython, I believe) you can skip
# the file save intermediate, see https://git.io/vh8Mi for similar advanced 
# examples.
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
# If you are Wayne, see `Collecting yeast XXXXX XXXXXX (XXXXXX) May 2018.md` for 
# impetus behind this script.
#
#
#
# Dependencies beyond the mostly standard libraries/modules:
#
#
#
# VERSION HISTORY:
# v.0.1. basic working version

#
# To do:
# - any way to make it work with a string?
#
#
#
#
# TO RUN:
# Examples,
# Enter on the command line of your terminal, the line
#-----------------------------------
# python reverse_complement_of_clustal_alignment.py ALIGNMENT_TEXT_FILE
#-----------------------------------
#
# Issue `reverse_complement_of_clustal_alignment.py -h` for details.
# 
#
#
# To use this after importing/pasting or loading into a cell in a Jupyter 
# notebook, specify at least the results file (or results as a string) in the 
# call to the main function similar to below:
# reverse_complement_of_clustal_alignment("test.clustal")
# 
#-or-
# To specify OPTIONAL source of description line for the description lines made
# by this script you can use 'descr_source' to indicate the file of ungapped 
# extracted FASTA-formatted seqeuence outout by 
# `extract_regions_from_clustal_alignment.py` similar to below:
# reverse_complement_of_clustal_alignment("test.clustal",descr_source="BlockA_extracted_ungapped.fa")
#
# 
#
'''
CURRENT ACTUAL CODE FOR RUNNING/TESTING IN A NOTEBOOK WHEN IMPORTED/LOADED OR 
PASTED IN ANOTHER CELL:
reverse_complement_of_clustal_alignment("test.clustal")
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


suffix_for_saving = "_rev_compl" #to be used for naming the output automatically
# when running script from command line to act on an input file

gap_indicator = "-" #change to "." if using dots. Used by 
# `get_aln_index_and_real_pos()` and `.ungap()`.



#
#*******************************************************************************
#**********************END USER ADJUSTABLE VARIABLES****************************













#*******************************************************************************
#*******************************************************************************
###DO NOT EDIT BELOW HERE - ENTER VALUES ABOVE###

import sys
import os
from collections import Counter
from Bio import SeqIO
from Bio import AlignIO
from Bio.Alphabet import generic_dna 
from Bio.Seq import Seq 
from Bio.SeqRecord import SeqRecord 
from Bio.Align import MultipleSeqAlignment 




###---------------------------HELPER FUNCTIONS---------------------------------###

def get_aln_index_and_real_pos(sequence):
    '''
    a generator that takes a sequence and then returns next position in 
    alignment (equivalent to index, i.e., zero-indexed) and the actual position 
    in the contiguous sequence to which that corresponds. So the second-value
    is the ungapped position in common terms.
    For example, if sequence is `a-b`,the first values returned are `0,1`.
    The second values returned are `1,1`. And third returned values are `2,2`.
    '''
    indx = 0
    while True:
        yield indx, len(sequence[:indx+1].replace(gap_indicator,"")) #because second 
        # value after colon means up to but not including, the `+1` there allows
        # getting first character when index is set to zero
        indx+=1

def generate_output_file_name(file_name,suffix_for_saving, fa=False, tb = False):
    '''
    Takes a file name as an argument and returns string for the name of the
    output file. The generated name is based on the original file
    name.

    Added `fa` parameter so can say if FASTA and then make file extension `.fa`.
    Added `tb` parameter so can say if tab-delimited table and then make file 
    extension `.tsv`.

    Specific example
    =================
    Calling function with
        ("alignment.clustal")
    returns
        "alignment_ADJUSTED.clustal"
    '''
    main_part_of_name, file_extension = os.path.splitext(
        file_name) #from 
    #http://stackoverflow.com/questions/541390/extracting-extension-from-filename-in-python
    if '.' in file_name:  #I don't know if this is needed with the os.path.splitext method but I had it before so left it
        if fa:
            return main_part_of_name + suffix_for_saving  + ".fa"
        elif tb:
            return main_part_of_name + suffix_for_saving  + ".tsv"
        else:
            return main_part_of_name + suffix_for_saving  + file_extension
    else:
        if fa:
            return file_name + suffix_for_saving + ".fa"
        elif tb:
            return main_part_of_name + suffix_for_saving  + ".tsv"
        else:
            return file_name + suffix_for_saving + ".clustal"


def any_first_words_occur_four_times(repeated_words_list):
    '''
    Return True if any words in the list occur four times
    '''
    most_common,num_most_common = Counter(
        repeated_words_list).most_common(1)[0] # based on
        # https://stackoverflow.com/a/6987358/8508004
    return num_most_common >= 4

###--------------------------END OF HELPER FUNCTIONS---------------------------###
###--------------------------END OF HELPER FUNCTIONS---------------------------###



#*******************************************************************************
###------------------------'main' function of script---------------------------##

def reverse_complement_of_clustal_alignment(
    alignment, descr_source = None,
    suffix_for_saving = suffix_for_saving, name_basis="alignment.clustal"):
    '''
    Main function of script. 
    It will take an alignment text file (clustal form) (or the alignment as a 
    Python string BUT THIS DOESN'T SEEM TO WORK BECAUSE GET OSERROR ABOUT File
    NAME EVEN THOUGH EXCEPT `OSError` already there.) 
    and put sequences within the alignment in reverse complement.

    The option to provide the results as a string is to handle where sending the
    data directly from shell script to Python without a typical file 
    intermediate, see the advanced notebook at https://git.io/vpr7i for 
    examples. The obvious use case for that is when working in the Jupyter 
    # notebook environment.
    '''
    # feedback about important option
    sys.stderr.write("\n**NOTE: gap indicator in this script is currently set "
        "to '{}'. If\nthat does not match what provided alignment uses to "
        "indicate gaps,\nplease change the setting within the script code "
        "under\n'USER ADJUSTABLE VALUES' around line 120 "
        "(give or take a few).**\n".format(gap_indicator))


    # Bring in the necessary alignment data:
    #---------------------------------------------------------------------------

    try:
        file_name = alignment
        alignment = AlignIO.read(alignment, "clustal") 
    except (TypeError,OSError,IOError) as e:
        file_name = name_basis
        try:
            from StringIO import StringIO
        except ImportError:
            from io import StringIO
        alignment = AlignIO.read(StringIO(alignment), "clustal") 
        # Note "FileNotFoundError is a subclass of OSError"
        # (https://stackoverflow.com/a/28633573/8508004)

    # feedback
    sys.stderr.write("Alignment file read...")



    # Go through multiple sequence alignment making reverse complement for each,
    # plus grab id for each for use later when define SeqRecord
    #---------------------------------------------------------------------------
    seq_n_id_list = []
    for record in alignment:
        seq_n_id_list.append((record.seq.reverse_complement(),record.id))


        
    





    # Save aligned sequence region in clustal format
    # based on http://biopython.org/DIST/docs/api/Bio.Align-pysrc.html
    # where it says "You would normally load a MSA from a file using 
    # Bio.AlignIO, but you can do this from a list of SeqRecord objects too:"
    # Note that each `record` in `record in aligmnet` is a SeqRecord but I 
    # kept getting errors when trying to use that SeqRecord to try and feed
    # it to `align = MultipleSeqAlignment()` directly because
    # was a probelm with 'New seq object is not a seq' or something like that 
    # (issue with alphabet?), it seemed. Lucily, I knew with text had
    # worked in `extract_regions_from_clustal_alignment.py`
    records = []
    for seq_n_id in seq_n_id_list:
        records.append(
            SeqRecord(Seq(str(seq_n_id[0]), generic_dna), seq_n_id[1])) # based
        # on https://www.biostars.org/p/48797/ and worked in 
        # `extract_regions_from_clustal_alignment.py`
    align = MultipleSeqAlignment(records)
    out_alignment_name = generate_output_file_name(file_name,suffix_for_saving)
    AlignIO.write(align, out_alignment_name, "clustal") # based on 
    # https://biopython.org/wiki/AlignIO

    # Feedback
    sys.stderr.write("\n\nAlignment block converted to reverse complement and "
        "saved as '{}'.".format(out_alignment_name))



    # Extract descriptions from provided file, if provided
    if descr_source:
        ids_n_descriptions = {}
        # read in the fasta-formatted file
        with open(descr_source, 'r') as input_handler:
            for line in input_handler:
                if line.startswith('>'):
                    line_parts = line.strip().split(' ',1) # maxsplit for
                    # just first word at index zero
                    ids_n_descriptions[line_parts[0][1:]] = line_parts[1]#`[1:]`
                    # after `line_parts[0]` so that the first character, which
                    # is `>` since used `line.startswith('>')`, gets dropped.
                    # so ids match and don't have `>identifier` at instead.
        # use that in description for new ungapped versions


    # Make raw, ungapped fasta of each aligned sequence.
    # Seems  need `.ungap()` method seq.io includes
    ungapped_records = []
    for seq_n_id in seq_n_id_list:
        if descr_source:
            # use `ids_n_descriptions` to get appropriate description text
            corrspnd_description=ids_n_descriptions[seq_n_id[1]]
            ungapped_records.append(
                SeqRecord(
                Seq(str(seq_n_id[0]), generic_dna).ungap(
                gap_indicator), seq_n_id[1], description=corrspnd_description))# 
            # basedon https://www.biostars.org/p/48797/ and `.ungap()` method, 
            # see https://github.com/biopython/biopython/issues/1511 , and 
            # `description` from from what I seen for `id` plus 
            # https://biopython.org/wiki/SeqIO
        else:
            ungapped_records.append(
                SeqRecord(
                Seq(str(seq_n_id[0]), generic_dna).ungap(
                gap_indicator), seq_n_id[1], description=""))# based
            # on https://www.biostars.org/p/48797/ and `.ungap()` method, see
            # https://github.com/biopython/biopython/issues/1511 , and 
            # `description` from from what I seen for `id` plus 
            # https://biopython.org/wiki/SeqIO

    # save records as one multi-fasta file
    ungapped_file_name = generate_output_file_name(
        file_name,suffix_for_saving+"_ungapped",fa=True)
    SeqIO.write(ungapped_records,ungapped_file_name, "fasta");
    # Feedback
    sys.stderr.write("\n\nAlignment block converted to reverse complement and "
        "saved in\nFASTA-formatted, ungapped form "
        "in '{}'.".format(ungapped_file_name))





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
    kwargs['descr_source'] = descr_source
    kwargs['suffix_for_saving'] = suffix_for_saving
    reverse_complement_of_clustal_alignment(alignment,**kwargs)
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
    ###-----------------for parsing command line arguments-----------------------###
    import argparse
    parser = argparse.ArgumentParser(prog=
        'reverse_complement_of_clustal_alignment.py',
        description="reverse_complement_of_clustal_alignment.py \
        takes a text document of an alignment in CLUSTAL format and outputs \
        that same alignment, but with composing sequences as reverse \
        complements.\
        **** Script by Wayne Decatur   \
        (fomightez @ github) ***")

    parser.add_argument("align_file", help="Name of file of alignment text \
        file.\
        ", metavar="ALIGNMENT_FILE")


    parser.add_argument('-ds', '--descr_source', action='store', type=str, 
    default= None, help="OPTIONAL: Provide FASTA file with same \
    ids to use as source of alignment to use as source of descriptions. The \
    typical file is output as `*_extracted_ungapped.fa` files from the \
    `extract_regions_from_clustal_alignment.py` script\
    .")

    parser.add_argument('-os', '--output_suffix', action='store', type=str, 
    default= suffix_for_saving, help="OPTIONAL: Set a suffix for including in file \
    name of output. \
    If none provided, '{}' will be used.".format(suffix_for_saving))




    #I would also like trigger help to display if no arguments provided because 
    # need at least one for url
    if len(sys.argv)==1:    #from http://stackoverflow.com/questions/4042452/display-help-message-with-python-argparse-when-script-is-called-without-any-argu
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    alignment = args.align_file
    descr_source = args.descr_source
    suffix_for_saving = args.output_suffix


    main()

#*******************************************************************************
###-***********************END MAIN PORTION OF SCRIPT***********************-###
#*******************************************************************************
