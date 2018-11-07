#!/usr/bin/env python
# score_sequences_in_clustal_msa_favoring_top_line.py 
__author__ = "Wayne Decatur" #fomightez on GitHub
__license__ = "MIT"
__version__ = "0.1.0"


# score_sequences_in_clustal_msa_favoring_top_line.py  by 
# Wayne Decatur
# ver 0.1
#
#*******************************************************************************
# Verified compatible with both Python 2.7 and Python 3.6; written initially in 
# Python 3. 
#
#
# PURPOSE: Takes a multiple sequence alignment (MSA) in clustal format and 
# generates a quick-n-dirty 'score' of percent match to the top line sequence 
# for each individual sequence in the alignment. Mainly meant to be a rough 
# guide to replace (or confirm) visual inspection. Returns a dataframe with the
# id and score for each sequence in the alignment.
# It is not exactly relative the top line because gaps in the top line are not
# counted as conserved and thus even the top line won't be 100% if in the 
# alignment the top line has gaps; this was a compromise to get those sequences
# that are all gaps to score as `0.0`, which makes more sense then them 
# accumulating score where the gaps happen to match with gaps in top line. 
# Hence `favoring` in name and not `relative`.
# This is only meant to be run when `score_sequences_in_clustal_msa.py` fails 
# and you have additional knowledge that the sequence on the top line of the 
# multiple sequence alignment is a good guide, i.e., it is a reference
# seqeunce. The tell-tale sign that `score_sequences_in_clustal_msa.py` failed 
# is when all sequences have scores of 0.0. `score_sequences_in_clustal_msa.py`
# scores purely on a match to the majority; however, if the majority of 
# sequences included in an alignment are unrelated to several sequences, the 
# subset of 'related' sequences will be scored poorly even though they really 
# are the subject of interest. This script allows placing the subject 
# of interest on the first line of the mutliple sequence alignment to indicate 
# it is to be considered as 'the reference'. (This can be done, for example, by 
# using `input` as the order in Clustal Omega by selecting that under 
# `More options` at https://www.ebi.ac.uk/Tools/msa/clustalo/ when submitting 
# sequences for alignment.)
#
#
# 
#
# Written to run from command line or imported into/pasted/loaded inside a 
# Jupyter notebook cell.  
#
#
#
# If you are Wayne, see `Gathering XXXXXXs from aligned XXXXXXXXXXXXX genomes 
# collected from the 1011 genomes project.ipynb` for specific impetus behind 
# this script.
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
#
#
#
# VERSION HISTORY:
# v.0.1. basic working version

#
# to do:
# - since this is based on `score_sequences_in_clustal_msa.py` it inherits
# a lot of the 'to do's of that,
# - make so works with alignments where length is less than (2 * number 
    # characters in aligned block per line) + 1. I.E. two blocks plus one more 
    #bp.  The reason for this now is that it 
    # uses first words occuring at least two times to distiguish indentifiers from 
    # any header line and plus it needs to hit third occurence of first id in 
    # order to trigger this. Hence the `+ 1` part of  `(2 times + 1)`. Ideas
    # to get this would be to rework how this is done to only need one complete  
    # block and identify the first_id first as the first firs_word that occurs 
    # twice. Then reset to reread the file from the start to collect 
    # first_words between the first occurence of first_id and the second 
    # occurence of it (exlcude symbols normally used in consensus line). Then 
    # add in a way to set ids by hand in the script for case one block.
    # When I do this in the comments for the modified version that this was 
    # adapted from the original approach in 
    # `extract_regions_from_clustal_alignment.py` that relied on more repeats. 
    # There I was making the assumption I'd have lots of repeats because it 
    # would probably be extracting a good size alignment from a larger context 
    # (like an alignment of entire region or chromosome) because if smaller 
    # could do by hand and so I could just continue on processing several 
    # blocks to clearly delineate the ids. This is a modified version of that 
    # approach that tries to default to not needing hard coding of the ids even 
    # if alignment is rather short, but provides a way to hardcode them if 
    # the alignment is indeed too short to make clear which are ids and which 
    # are extraneous text.
    # REMOVE CAVEAT WHEN DONE. 
# - add `name_basis` use?
#   
#
#
#
#
# TO RUN:
# Examples,
# Enter on the command line of your terminal, the line
#-----------------------------------
# python score_sequences_in_clustal_msa_favoring_top_line.py msa.txt 
#-----------------------------------
#
# Issue `score_sequences_in_clustal_msa_favoring_top_line.py  -h` for details.
# 
#
#
# To use this after importing/pasting or loading into a cell in a Jupyter 
# notebook, specify at least the results file (or results as a string) in the 
# call to the main function similar to below:
# score_sequences_in_clustal_msa_favoring_top_line("msa.clustal")
# 
#-or-
# To specify OPTIONAL 'file name'-like string to base names of generated files 
# on, specify `name_basis` like on next line when providing alignment as string:
# score_sequences_in_clustal_msa_favoring_top_line(test.out,name_basis="rt_arm_chrV_sensu.clustal")
# **`name_basis` ignored when a real file is supplied; presently, 
# "alignment_plusCONS.clustal" is default.** 
#
# 
#
'''
CURRENT ACTUAL CODE FOR RUNNING/TESTING IN A NOTEBOOK WHEN IMPORTED/LOADED OR 
PASTED IN ANOTHER CELL:
from score_sequences_in_clustal_msa_favoring_top_line import score_sequences_in_clustal_msa_favoring_top_line
score_sequences_in_clustal_msa_favoring_top_line("msa.clustal")
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


suffix_for_saving = "_tl_match_scores" #to be used for naming the output 
# automatically when running script from command line to act on an input file



#
#*******************************************************************************
#**********************END USER ADJUSTABLE VARIABLES****************************













#*******************************************************************************
#*******************************************************************************
###DO NOT EDIT BELOW HERE - ENTER VALUES ABOVE###

import sys
import os
from collections import Counter




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
        yield indx, len(sequence[:indx+1].replace("-","")) #because second 
        # value after colon means up to but not including, the `+1` there allows
        # getting first character when index is set to zero
        indx+=1

def generate_output_file_name(file_name,suffix_for_saving):
    '''
    Takes a file name as an argument and returns string for the name of the
    output file. The generated name is based on the original file
    name.


    Specific example
    =================
    Calling function with
        ("alignment.clustal")
    returns
        "alignment_tl_match_scores.tsv"
    '''
    main_part_of_name, file_extension = os.path.splitext(
        file_name) #from 
    #http://stackoverflow.com/questions/541390/extracting-extension-from-filename-in-python
    if '.' in file_name:  #I don't know if this is needed with the os.path.splitext method but I had it before so left it
        return main_part_of_name + suffix_for_saving  + ".tsv"
    else:
        return file_name + suffix_for_saving + ".tsv"


def any_first_words_occur_three_times(repeated_words_list):
    '''
    Return True if any words in the list occur three times
    '''
    most_common,num_most_common = Counter(
        repeated_words_list).most_common(1)[0] # based on
        # https://stackoverflow.com/a/6987358/8508004
    return num_most_common >= 3

###--------------------------END OF HELPER FUNCTIONS---------------------------###
###--------------------------END OF HELPER FUNCTIONS---------------------------###






#*******************************************************************************
###------------------------'main' function of script---------------------------##

def score_sequences_in_clustal_msa_favoring_top_line(
    alignment, 
    suffix_for_saving = suffix_for_saving, name_basis="alignment.clustal",
    return_df = True):
    '''
    Main function of script. 
    It will take an alignment text file (or the alignment as a Python string(?)) 
    and add position indicators to the top line that is the actual numbering of 
    the residues in contiguous sequence.
    Optionally also returns a dataframe of the results data. Meant for use in a 
    Jupyter notebook.

    DOES THIS WORK???--->The option to provide the results as a string is to 
    handle where sending the
    data directly from shell script to Python without a typical file 
    intermediate, see the advanced notebook at https://git.io/vpr7i for 
    examples. The obvious use case for that is when working in the Jupyter 
    notebook environment. <---DOES THIS WORK???
    '''
    # Bring in the necessary alignment data:
    #---------------------------------------------------------------------------

    try:
        with open(alignment, 'r') as the_file:
            file_name = alignment
            alignment = the_file.read()
    except (TypeError,OSError,IOError) as e:
        file_name = name_basis
        pass # pass because instead just want to use results as a string because
        # probably provided as a string directly piped from PatMatch to Python
        # without file intermediate. The results variable will already be a
        # string in that case and so ready to go and don't need to read file.
        # The `except` above the pass has three errors it excepts because the
        # first was a TypeError seen when I tried to pass a file-like string to
        # `with open` because it seems `with open` method is incompatible with
        # use of StringIO, I think, which I usually use to try to pass things 
        # associated with file methods string. (I qualifed it with 'I think' b/c 
        # questions on stackoverlflow seemed to agree but I didn't try every 
        # possibility because realized this would probably be a better way to 
        # handle anyway.) That TypeError except got me to the next issue which
        # was trying the string as a file name and getting it was too long, and 
        # so I added the `OSError` catch and that seemed to make passing a 
        # string into the function work. `IOError` seemed to handle that same
        # thing in Python 2.7.
        # Note "FileNotFoundError is a subclass of OSError"
        # (https://stackoverflow.com/a/28633573/8508004)

    # feedback
    sys.stderr.write("Alignment file read...")



    # Go through multiple sequence alignment file parsing it to collect the 
    # identifiers. Also, identify the first identifier in each alignment block.
    #---------------------------------------------------------------------------
    first_words = []
    first_id_needed = True
    for line in alignment.split("\n"):
        if line:
            first_word = line.split()[0]
            if first_word in first_words:
                if first_id_needed:
                    first_identifier = first_word
                    first_id_needed = False
                first_words.append(first_word)
                if any_first_words_occur_three_times(first_words):
                    # want to collect those that occur at least two times 
                    # because if an identifier has occured three times, than 
                    # others should number as two and don't want anything 
                    # occurring less since usually there is a header line in
                    # Clustal alignments describing source (and/or version)
                    the_count = Counter(first_words)
                    aln_ids = [k for k, v in the_count.items() if v > 1] # based
                    # on https://stackoverflow.com/a/26773120/8508004 and
                    # https://stackoverflow.com/a/30418498/8508004, to work in 
                    # 2.7 and 3
                    last_identifier = first_words[-2]
                    break # because have gone far enough to collect the identifiers
            else:
                first_words.append(first_word)
    # feedback
    sys.stderr.write(
        "top line identifier determined as '{}'...".format(first_identifier))
    # NOTE WILL GET `UnboundLocalError: local variable 'first_identifier' 
    # referenced before assignment` if no alignment/file present!!!!
    sys.stderr.write(
        "bottom line identifier determined as '{}'."
        "..".format(last_identifier))




    # Go through multiple sequence alignment file parsing it as needed to 
    # accumulate full alignment for each sequence identifier code.
    # Since collected identifiers above, can use them.
    #---------------------------------------------------------------------------
    alignment_dict = {}
    num_chars_in_first_MSA_block = 0  #best to collect since differs in MSAs
    for line in alignment.split("\n"):
        # determing if it is one of the alignment lines
        for id_ in aln_ids:
            if line.strip().startswith(id_):
                if id_ in alignment_dict:
                    alignment_dict[id_] += line[len(
                        id_):].strip().split()[0].strip()
                else:
                    alignment_dict[id_] =  line[len(
                        id_):].strip().split()[0].strip()
                    if num_chars_in_first_MSA_block == 0:
                        num_chars_in_first_MSA_block = len(alignment_dict[id_])

    # Because: "Note the website should have an option about showing gaps as 
    # periods (dots) or dashes, we've shown dashes above." - SOURCE:
    # http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc75
    # and I am adjusting sequence to allow for that by changing from '.' to
    # dashes.
    # First I had check if any have periods, but maybe just faster/easier to 
    # run `.replace()` on all anyway? Rather then check and then process.
    '''
    if any([True if '.' in v else False for k,v in alignment_dict.items()]):
        alignment_dict= {k: v.replace(".","-") for k,v in alignment_dict.items()}
    '''
    alignment_dict= {k: v.replace(".","-") for k,v in alignment_dict.items()}

    # feedback
    sys.stderr.write("\nindividual lines for each sequence identifier parsed...")





    # sanity check
    # Verify that alignments are equal in length
    lengths_of_alignments = [len(v) for v in alignment_dict.values()]
    assert lengths_of_alignments.count(
        lengths_of_alignments[0]) == len(lengths_of_alignments), ("The length "
        "of all parsed alignments should be the same.") # see where 
        # `conservation_line` being made in 
        # `pretty_msa_maker_from_clustal_nucleic.py ` for mechanism used in 
        # assertion test involving Ivo van der Wijk's solution from 
        # https://stackoverflow.com/a/3844948/8508004

    len_alignment = lengths_of_alignments[0] #-or- use `len(alignment[0])`



    # Now with full alignments in hand, calculate conservation.
    #---------------------------------------------------------------------------
    # Here top line, excluding gapped regions will be used in the scoring!
    sys.stderr.write("determining what matches top line at each position of "
        "aligned sequences, exlcuding gaps in top line...")
    majority = ''
    for indx, nt in enumerate(alignment_dict[aln_ids[0]]):
        if alignment_dict[first_identifier][indx] == "-":
            majority += " "
        else:
            majority += alignment_dict[first_identifier][indx]


    # Now with the 'majority' line in hand, calculate scores for each sequence
    #---------------------------------------------------------------------------
    scores_dict = {}
    for id_ in alignment_dict:
        seq = alignment_dict[id_]
        # Go though letter by letter of aligned sequence and increment counter 
        # if appropriate
        matches= 0
        for indx, nt in enumerate(seq):
            if majority[indx] == nt:
                matches += 1
        scores_dict[id_] = matches/len_alignment


    # Save a table of the scores
    #---------------------------------------------------------------------------
    import pandas as pd
    scores_df = pd.DataFrame(list(scores_dict.items()),columns = ['id','score']) 
    output_file_name = generate_output_file_name(
        file_name,suffix_for_saving)
    scores_df.to_csv(output_file_name, sep='\t',index = False)
    # Feedback
    sys.stderr.write("\n\nScores (%) for how well each sequence matched to "
        "the top line saved as\n"
        "a table '{}'.".format(output_file_name))


    # Return dataframe (optional)
    #---------------------------------------------------------------------------
    if return_df:
        sys.stderr.write( "\n\nReturning a dataframe with the scores "
                "as well.")
        return scores_df

    sys.stderr.write("\nFinished.\n")







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
    kwargs['return_df'] = False #probably don't want dataframe returned if 
    # calling script from command line
    score_sequences_in_clustal_msa_favoring_top_line(alignment,**kwargs)
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
        'score_sequences_in_clustal_msa_favoring_top_line.py',
        description="score_sequences_in_clustal_msa_favoring_top_line.py  \
        takes a multiple sequence alignment (MSA) in clustal format and \
        generates a quick-n-dirty 'score' of percent match to the top \
        line sequence in the alignment. Gaps are counted as penalties from the \
        top possible score of 100%, even for the top line sequence. Mainly \
        meant to be a \
        rough guide to replace (or confirm) visual inspection. Returns a \
        dataframe with the id and score for each sequence in the alignment. \
        **** Script by Wayne Decatur   \
        (fomightez @ github) ***")

    parser.add_argument("MSA", help="Name of file containing multiple sequence \
        alignment in CLUSTAL format. REQUIRED. \
        ", metavar="MSA_FILE")

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
    alignment = args.MSA
    suffix_for_saving = args.output_suffix


    main()

#*******************************************************************************
###-***********************END MAIN PORTION OF SCRIPT***********************-###
#*******************************************************************************
