#!/usr/bin/env python
# score_columns_in_clustal_nucleic.py
__author__ = "Wayne Decatur" #fomightez on GitHub
__license__ = "MIT"
__version__ = "0.1.0"


# score_columns_in_clustal_nucleic.py by 
# Wayne Decatur
# ver 0.1
#
#*******************************************************************************
# Verified compatible with both Python 2.7 and Python 3.6; written initially in 
# Python 3. 
#
#
# PURPOSE: Takes a text document of an alignment in CLUSTAL format and generates
# a couple of quick-n-dirty 'scores' based on percent overall conservation in an 
# alignment. Specifically, determines conservation of residues in alignment 
# columns both including and excluding counting gaps in the column. This is done
# in order to assign some rough scoring metric(s) to an alignment. The metric(s)
# can then be used when comparing separate alignments of related, but distinct 
# repetitive genetic elements in order to objectively assess which of the 
# related elements are more highly conserved overall given what a particular
# 'score' emphasizes. Mainly meant to be a rough guide to replace (or confirm) 
# visual inspection.
#
# Excluding gaps in the column is meant to allow insertions in a few sequences 
# not to be penalized as much and instead allow greater impact from segments
# that may match well.
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
# If you are Wayne, see `scoring XXXXXXX XXXXXX among strains in multiple 
# sequence alignments.md` for specific impetus behind this script.
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
# - 
#
#
#
#
# TO RUN:
# Examples,
# Enter on the command line of your terminal, the line
#-----------------------------------
# python score_columns_in_clustal_nucleic.py -ALIGNMENT_TEXT_FILE
#-----------------------------------
#
# Issue `score_columns_in_clustal_nucleic.py -h` for details.
# 
#
#
# To use this after importing/pasting or loading into a cell in a Jupyter 
# notebook, specify at least the results file (or results as a string) in the 
# call to the main function similar to below:
# score_columns_in_clustal_nucleic("test.clustal")
# 
#-or-
# To specify OPTIONAL ....
#
# 
#
'''
CURRENT ACTUAL CODE FOR RUNNING/TESTING IN A NOTEBOOK WHEN IMPORTED/LOADED OR 
PASTED IN ANOTHER CELL:
score_columns_in_clustal_nucleic("test.clustal")
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


suffix_for_saving = "_qnd_scores"#to be used for naming the output automatically
# when running script from command line to act on an input file

gap_indicator = "-" #change to "." if using dots. Used by 
# `get_aln_index_and_real_pos()` and `.ungap()`.

percent_threshold = 1.0 # percent that matches need to match in column to be
# counted as conserved



#
#*******************************************************************************
#**********************END USER ADJUSTABLE VARIABLES****************************













#*******************************************************************************
#*******************************************************************************
###DO NOT EDIT BELOW HERE - ENTER VALUES ABOVE###

import sys
import os
from collections import Counter
from Bio import AlignIO





###---------------------------HELPER FUNCTIONS-------------------------------###

def increment_warranted(column_seq, count_gaps=True):
    '''
    Takes a column from a multiple sequence alignment and returns a `True` if
    the percent of matches of residues meets the `percent_threshold`.

    If `count_gaps=False`, then rows with gaps are not considered in trying
    to determine if conserved. This setting would lessen impact from an insert
    in the sequence.
    No matter what `count_gaps` and `percent_threshold` are, gaps can never be
    considered conserved
    '''
    if not count_gaps:
        column_seq = column_seq.replace(gap_indicator,"")
    total_in_column = len(column_seq)
    match_percent = (
        Counter(column_seq.replace(gap_indicator,"")).most_common(
        1)[0][1])/float(total_in_column)
    # note that because of the inclusion of 
    # `column_seq.replace(gap_indicator,"")` gaps can never be considered as 
    # most abundant. That may becaome important if ever lower threshold from 
    # 100% might because with `Counter(column_seq)` nine gap characters and one 
    # nucleotide in a column would be deemed conserved.
    if match_percent >= percent_threshold:
        return True
    else:
        return False


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




###--------------------------END OF HELPER FUNCTIONS---------------------------###
###--------------------------END OF HELPER FUNCTIONS---------------------------###



#*******************************************************************************
###------------------------'main' function of script---------------------------##

def score_columns_in_clustal_nucleic(
    alignment, 
    suffix_for_saving = suffix_for_saving, name_basis="alignment.clustal"):
    '''
    Main function of script. 
    It will take an alignment text file (clustal form) (or the alignment as a 
    Python string BUT THIS DOESN'T SEEM TO WORK BECAUSE GET OSERROR ABOUT File
    NAME EVEN THOUGH EXCEPT `OSError` already there.) 
    and return some quick-n-dirty scoring assessments that reflect conservation.

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



    # Extract length of alignment
    #---------------------------------------------------------------------------
    # Sanity check first:
    # Verify that alignments are equal in length
    lengths_of_alignments = [len(record) for record in alignment]
    assert lengths_of_alignments.count(
        lengths_of_alignments[0]) == len(lengths_of_alignments), ("The length "
        "of all parsed alignments should be the same.") # see where 
        # `conservation_line` being made in 
        # `pretty_msa_maker_from_clustal_nucleic.py ` for mechanism used in 
        # assertion test involving Ivo van der Wijk's solution from 
        # https://stackoverflow.com/a/3844948/8508004

    len_alignment = lengths_of_alignments[0] #-or- use `len(alignment[0])`

    

    # Go though column by column of alignment and increment counter if 
    # appropriate
    #---------------------------------------------------------------------------
    matches_with_gaps = 0
    matches_without_gaps = 0
    for column_no in range(len_alignment):
        column_seq = alignment[:, column_no]
        if increment_warranted(column_seq, count_gaps=True):
            matches_with_gaps += 1
        if increment_warranted(column_seq, count_gaps=False):
            matches_without_gaps += 1



    # Normalize the collected scores based on alignment length
    #---------------------------------------------------------------------------
    score_with_gaps = matches_with_gaps/len_alignment
    score_not_counting_gaps = matches_without_gaps/len_alignment

        
    





    # Feedback
    #---------------------------------------------------------------------------
    # report scores
    sys.stderr.write("\n\nSimilarity percent including gaps is:\n"
        "'{}'.".format(score_with_gaps))
    sys.stderr.write("\nSimilarity percent not counting gaps is:\n"
        "'{}'.".format(score_not_counting_gaps))
    score_with_gaps
    score_not_counting_gaps

    # Save a table of the scores
    import pandas as pd
    scores_dict = {
        'alignment': file_name,
        'score_including_gaps': [score_with_gaps],
        'score_without_considering_gaps': [score_not_counting_gaps],
        }
    # Even though using single items for scores above in making dict, I have to 
    # make them a list or otherwise when try to make a dataframe with them in 
    # nextline get error: `ValueError: If using all scalar values, you must 
    # pass an index`. See https://stackoverflow.com/a/17840195/8508004 .
    scores_df = pd.DataFrame(
        scores_dict,columns = ['alignment',
        'score_including_gaps','score_without_considering_gaps']) 
    scores_file_name = generate_output_file_name(
        file_name,suffix_for_saving,tb=True)
    scores_df.to_csv(scores_file_name, sep='\t',index = False)
    # Feedback
    sys.stderr.write("\n\nScores for provided alignment saved as\n"
        "a table '{}'.".format(scores_file_name))





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
    #kwargs['descr_source'] = descr_source
    kwargs['suffix_for_saving'] = suffix_for_saving
    score_columns_in_clustal_nucleic(alignment,**kwargs)
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
        'score_columns_in_clustal_nucleic.py',
        description="score_columns_in_clustal_nucleic.py \
        takes a text document of an alignment in CLUSTAL format  and generates \
        a couple of quick-n-dirty 'scores' based on percent overall \
        conservation in an alignment. Specifically, determines conservation of \
        residues in alignment columns both including and excluding counting \
        gaps in the column. Intended to be used when comparing separate \
        alignments of related but distinct repetitive genetic elements in \
        order to objectively assess which of the related elements are more \
        highly conserved overall. \
        **** Script by Wayne Decatur   \
        (fomightez @ github) ***")

    parser.add_argument("align_file", help="Name of file of alignmnet text \
        file.\
        ", metavar="ALIGNMENT_FILE")


    '''parser.add_argument('-ds', '--descr_source', action='store', type=str, 
    default= None, help="OPTIONAL: Provide FASTA file with same \
    ids to use as source of alignment to use as source of descriptions. The \
    typical file is output as `*_extracted_ungapped.fa` files from the \
    `extract_regions_from_clustal_alignment.py` script\
    .")
    '''
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
    suffix_for_saving = args.output_suffix


    main()

#*******************************************************************************
###-***********************END MAIN PORTION OF SCRIPT***********************-###
#*******************************************************************************
