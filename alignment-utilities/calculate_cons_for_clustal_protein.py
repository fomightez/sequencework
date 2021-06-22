#!/usr/bin/env python
# calculate_cons_for_clustal_protein.py 
__author__ = "Wayne Decatur" #fomightez on GitHub
__license__ = "MIT"
__version__ = "0.1.0"


# calculate_cons_for_clustal_protein.py  by 
# Wayne Decatur
# ver 0.1
#
#*******************************************************************************
# Verified compatible with both Python 2.7 and Python 3.6; written initially in 
# Python 3. 
#
#
# PURPOSE: Takes a multiple sequence alignment (MSA) in Clustal format and adds 
# symbols indicating conservation/consensus below each set of sequence 
# alignment blocks.  Different symbols indicate identical residues and conserved 
# substitutions in a output alignment file.
# Good for replacing absent information in alignments. Or good
# for replacing outdated /lost information in alignments upon hand-editing with 
# Seqotron on a Mac. Better than `cons` in EMBOSS suite because makes format 
# that matches t-coffee more so no hand editing, plus adds conserved 
# substitutions as colons and periods depending if weak or strongly similar, see
# https://www.ebi.ac.uk/seqdb/confluence/display/JDSAT/Bioinformatics+Tools+FAQ#BioinformaticsToolsFAQ-WhatdoconsensussymbolsrepresentinaMultipleSequenceAlignment?
# 
# Also useful in conjunction with Kalign (https://www.ebi.ac.uk/Tools/msa/kalign/) 
# tool as that doesn't add asterisks indicating conservation. (Several others
# that generate Clustal formatted alignments also seem to save with consensus 
# symbol lines.)
# CAVIEAT: right now in order to correctly identify the sequence identifiers at
# the start of alignment lines, they need to occur at least three times so there
# has to be at least two full blocks of sequence + at least a partial block. The
# last one can be as small as a single bp. Since most alignment blocks default
# to 60, it means this won't work at present if less than 121 bps long sequence.
#
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
# If you are Wayne, see
# `publication-quality pretty version of alignment using Reportlab.md`
# and `edited t-coffee alignment for true xxxxxxxxxxxxx xxxxxxs.md`
# for info about this script and development.
# Wrote my own when didn't see easy way to customize or even add the result of
# Biopython's related method shown at 
# http://biopython.org/DIST/docs/api/Bio.Align.AlignInfo.SummaryInfo-class.html
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
#
#
#
# VERSION HISTORY:
# v.0.1. basic working version

#
# to do:
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
# python calculate_cons_for_clustal_protein.py msa.txt 
#-----------------------------------
#
# Issue `calculate_cons_for_clustal_protein.py  -h` for details.
# 
#
#
# To use this after importing/pasting or loading into a cell in a Jupyter 
# notebook, specify at least the results file (or results as a string) in the 
# call to the main function similar to below:
# calculate_cons_for_clustal_protein("msa.clustal")
# 
#-or-
# To specify OPTIONAL 'file name'-like string to base names of generated files 
# on, specify `name_basis` like on next line when providing alignment as string:
# calculate_cons_for_clustal_protein(test.out,name_basis="rt_arm_chrV_sensu.clustal")
# **`name_basis` ignored when a real file is supplied; presently, 
# "alignment_plusCONS.clustal" is default.** 
#
# 
#
'''
CURRENT ACTUAL CODE FOR RUNNING/TESTING IN A NOTEBOOK WHEN IMPORTED/LOADED OR 
PASTED IN ANOTHER CELL:
calculate_cons_for_clustal_protein("msa.clustal")
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


suffix_for_saving = "_plusCONS" #to be used for naming the output automatically
# when running script from command line to act on an input file

strongly_similar_aa_tuples = [("S","T","A"), ("N","E","Q","K"), 
("N","H","Q","K"), ("N","D","E","Q"), ("Q","H","R","K"), ("M","I","L","V"),
("M","I","L","F"),("H","Y"),("F","Y","W")] #used to decide if
# conservative substitution (strongly similar) in case where not all identical.
# Based on
# https://www.ebi.ac.uk/seqdb/confluence/display/JDSAT/Bioinformatics+Tools+FAQ#BioinformaticsToolsFAQ-WhatdoconsensussymbolsrepresentinaMultipleSequenceAlignment?

weakly_similar_aa_tuples = [("C","C","A"), ("A","T","V"), ("S","A","G"), 
("S","T","N","K"), ("S","T","P","A"), ("S","G","N","D"), 
("S", "N", "D", "E","Q","K"), ("N", "D", "E", "Q","H","K"),
("N", "E", "H", "Q","R","K"), ("F", "V", "I", "L","M"), ("H","F","Y")] #used to 
# decide if conservative substitution (weakly similar) in case where not identical.
# Based on
# https://www.ebi.ac.uk/seqdb/confluence/display/JDSAT/Bioinformatics+Tools+FAQ#BioinformaticsToolsFAQ-WhatdoconsensussymbolsrepresentinaMultipleSequenceAlignment?


#
#*******************************************************************************
#**********************END USER ADJUSTABLE VARIABLES****************************













#*******************************************************************************
#*******************************************************************************
###DO NOT EDIT BELOW HERE - ENTER VALUES ABOVE###

import sys
import os
from collections import Counter




###-------------------------HELPER FUNCTIONS---------------------------------###

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
        "alignment_plusCONS.clustal"
    '''
    main_part_of_name, file_extension = os.path.splitext(
        file_name) #from 
    #http://stackoverflow.com/questions/541390/extracting-extension-from-filename-in-python
    if '.' in file_name:  #I don't know if this is needed with the os.path.splitext method but I had it before so left it
        return main_part_of_name + suffix_for_saving  + file_extension
    else:
        return file_name + suffix_for_saving + ".clustal"


def any_first_words_occur_three_times(repeated_words_list):
    '''
    Return True if any words in the list occur three times
    '''
    most_common,num_most_common = Counter(
        repeated_words_list).most_common(1)[0] # based on
        # https://stackoverflow.com/a/6987358/8508004
    return num_most_common >= 3

###--------------------------END OF HELPER FUNCTIONS-------------------------###
###--------------------------END OF HELPER FUNCTIONS-------------------------###






#*******************************************************************************
###------------------------'main' function of script--------------------------##

def calculate_cons_for_clustal_protein(
    alignment, 
    suffix_for_saving = suffix_for_saving, name_basis="alignment.clustal"):
    '''
    Main function of script. 
    It will take an alignment text file (or the alignment as a Python string(?)) 
    and add  symbols indicating identical residues and conserved substitutions, 
    and then output a new alignment file.

    The option to provide the results as a string is to handle where sending the
    data directly from shell script to Python without a typical file 
    intermediate, see the advanced notebook at https://git.io/vpr7i for 
    examples. The obvious use case for that is when working in the Jupyter 
    # notebook environment. BUt I don't think it works here??? maybe?
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
    # sys.stderr.write("Alignment file read...") # move to below so doesn't show
    # if file not really read



    # Go through multiple sequence alignment file parsing it to collect the 
    # identifiers. Also, identify the first identifier in each alignment block.
    #---------------------------------------------------------------------------
    first_words = []
    first_id_needed = True
    for line in alignment.split("\n"):
        if line.strip():
            first_word = line.split()[0]
            if first_word in first_words:
                if first_id_needed:
                    first_identifier = first_word
                    sys.stderr.write("Alignment file read...") 
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



    # Now with full alignments in hand, calculate conservation.
    #---------------------------------------------------------------------------
    # `conservation_line` will be used for logging if identical in all 
    # alignments at that position or if conservatively substituted (either 
        # strongly or weakly). 
    # `majority` will simply be what is the majority at that position and later 
    # be used in another script, i.e., my pretty msa maker one, for determining 
    # if nucleotides represented as bold. 
    sys.stderr.write("determining conservation of aligned sequences...")
    conservation_line = ""
    majority = ""
    for indx, nt in enumerate(alignment_dict[aln_ids[0]]):
        # make one by one a list of each nucleotide (or residue) at same position
        # in all alignments.
        current_pos_list = [nt]
        for x in range(1,len(aln_ids)):
            current_pos_list.append(alignment_dict[aln_ids[x]][indx])
        # If "-" present all then it isn't conserved in any way, and write space
        # to conservation line
        # if all are same, i.e. `x.count(x[0]) == len(x)`, where x is the list for
        # that position , then write asterisk to conservation line. (Uses
        # Ivo van der Wijk's solution from
        # https://stackoverflow.com/questions/3844801/check-if-all-elements-in-a-list-are-identical )
        # If all conserved substitution, write "." to conservation line.
        if '-' in current_pos_list:
            conservation_line += " "
        elif current_pos_list.count(current_pos_list[0]) == len(current_pos_list):
            conservation_line += "*"
        elif any([set(current_pos_list).issubset(x) for x in strongly_similar_aa_tuples]):
            conservation_line+= ":"
        elif any([set(current_pos_list).issubset(x) for x in weakly_similar_aa_tuples]):
            conservation_line+= "."
        else:
            conservation_line += " "
        if Counter(current_pos_list).most_common(1)[0][0] == "-":
            majority += " "
        else:
            most_common,num_most_common = Counter(
                current_pos_list).most_common(1)[0] # based on
            # https://stackoverflow.com/questions/6987285/python-find-the-item-with-maximum-occurrences-in-a-list
            if num_most_common > len(current_pos_list)/2:
                majority += most_common
            else:
                majority += " "





    # Add the conservation/consensus indicator symbols line to the MSA file below 
    # each block of sequences block and save the modified file.
    #---------------------------------------------------------------------------
    n = num_chars_in_first_MSA_block  
    chunks_for_conservation = [conservation_line[i:i+n] for i in range(
        0, len(conservation_line), n)]
    # above line based on satomacoto's answer at 
    # https://stackoverflow.com/questions/9475241/split-python-string-every-nth-character
    #print (list2text(chunks_for_conservation))
    chunk = iter(chunks_for_conservation)

    out_alignment_name = generate_output_file_name(file_name,suffix_for_saving)
    # prepare output file for saving so it will be open and ready
    with open(out_alignment_name, 'w') as output_file:

        # read in the alignment file line by line
        for line in alignment.split("\n"):
            if line.strip().startswith(last_identifier):
                # After line with last_identifier written to new file, add the 
                # next chunk of conservation to insert conservation indicator 
                # lines into MSA in output.
                output_file.write(line.strip()+"\n")
                # determine index of sequence
                seq = line.strip().split()[1].strip()
                seq_index = line.strip().index(seq)
                output_file.write((" "* seq_index)+next(chunk)+"\n")
            else:
                output_file.write(line.strip()+"\n")


    # Feedback
    sys.stderr.write("\n\nAlignment with conservation indication symbols added "
        "saved as '{}'.".format(out_alignment_name))
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
    calculate_cons_for_clustal_protein(alignment,**kwargs)
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
        'calculate_cons_for_clustal_protein.py ',
        description="calculate_cons_for_clustal_protein.py  \
        takes a multiple seuence alignment (MSA) in clustal format and adds \
        symbols indicating conservation/consensus below each set of sequence \
        alignment blocks . \
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
