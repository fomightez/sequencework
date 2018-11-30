#!/usr/bin/env python
# MSA_to_corresponding_residue_numbers.py
__author__ = "Wayne Decatur" #fomightez on GitHub
__license__ = "MIT"
__version__ = "0.1.0"


# MSA_to_corresponding_residue_numbers.py by 
# Wayne Decatur
# ver 0.1
#
#*******************************************************************************
# Verified compatible with both Python 2.7 and Python 3.6; written initially in 
# Python 3. 
#
#
# PURPOSE: Takes a multiple sequence alignment in Clustal format and determines
# the corresponding residue numbers for aligned residues for a specified 
# reference sequence in the alignment and the pairing of it with every other 
# aligned sequence. The reference sequence is specified by the identifier at 
# the start of each line in the sequence blocks.
# The sequences provided do not have to begin at reside one, but in that
# case the starting position represented for each sequence in the alignment must
# be provided. Even if all the others are simply `1`. The ends don't have to 
# represent full sequences either, and in that case nothing needs to be done if 
# they are truncated. They just won't be addressed. For the numberings extracted 
# to be correct, the only requirements is that the sequence between the first 
# and last residue be contiguous discounting gaps. The gaps don't matter but the
# sequence has to be intact once it starts are the specified residue number, 
# default of which is one. In other words, once the chain starts, no internal 
# stretches of residues or even a single residue can be absent. This can be 
# verified by using my script `check_seq_frag_in_MSAclustal_intact_viaFASTA.py`
# to compare a sequence in a multiple sequence alignment to that in a FASTA and
# validate it contains a continuous fragment. Optionally, 
# `check_seq_frag_in_MSAclustal_intact_viaFASTA.py` can return the value that
# corresponds to the first position of the fragment represented in the multiple
# sequence alignment, and so the information it yeilds can be very useful when 
# writing commands to call the script or main function when using fragments of 
# sequences in the multiple sequence alignment.
#
# Meant to be used to determine the specific corresponding pairs of residues of
# sequences represented in the alignment relative the sequence of a chain 
# available in a structure model. With these pairs of residues macromolecular 
# visualization commands, such as those that fit/dock/compare or highlight 
# corresponding residues via distinguishing colors and representations, can be 
# written. (For the latter, you may wish to also see 
# `categorize_residues_based_on_conservation_relative_consensus_line.py` as a 
# related & distinct alternative.)
# There a notebook demonstrating of use of this script to automagically write
# fit/dock/compare commands available via 
# https://github.com/fomightez/cl_demo-binder
# Note that all numbers returned by the script will be in common terms where
# the fist residue in a sequence is numbered ONE. (In other words, the 
# zero-indexing behind-the-scenes will not be visible to user.)
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
# If you are Wayne, see 
# `alignment of yeast XXX1 vs XXX2 vs XXX12 vs archael XXXX from Swiss Model.clw` 
# for specific impetus behind this script.
#
#
# Dependencies beyond the mostly standard libraries/modules:
# biopython
#
#
# VERSION HISTORY:
# v.0.1. basic working version

#
# To do:
#
#
#
# TO RUN:
# Examples,
# Enter on the command line of your terminal, the line
#-----------------------------------
# python MSA_to_corresponding_residue_numbers.py ALIGNMENT_TEXT_FILE reference_seq_id
#-----------------------------------
#
# Issue `MSA_to_corresponding_residue_numbers.py -h` for details.
# 
#
#
# To use this after importing/pasting or loading into a cell in a Jupyter 
# notebook, specify at least the alignment file, id of sequence to analyze, and 
# ungapped FASTA file of the sequence in the call to the main function similar 
# to below:
# ref_id, dfs_by_id = MSA_to_corresponding_residue_numbers("alignment.clw", "VPH1", return_dfs = True)
# -OR-
# If any of the sequences begin with a residue number other than ONE, provide
# a list of the starting positions for ALL sequences, even if all others are 
# ONE, as a python list. For example, if four sequences and the second one from 
# top one begins at residue 37, provide `supplied_start_pos = [1,37,1,1]` in 
# call to main function of script, like so:
# ref_id, dfs_by_id  = MSA_to_corresponding_residue_numbers("alignment.clw", "VPH1", return_dfs = True, supplied_start_pos = [5,1,9,8])
# 
#
# 
#
'''
CURRENT ACTUAL CODE FOR RUNNING/TESTING IN A NOTEBOOK WHEN IMPORTED/LOADED OR 
PASTED IN ANOTHER CELL:
ref_id, dfs_by_id  = MSA_to_corresponding_residue_numbers("alignment.clw", "VPH1", return_dfs = True)
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

gap_indicator = "-" #change to "." if using dots. Used by 
# `.ungap()`.



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
import pandas as pd





###---------------------------HELPER FUNCTIONS-------------------------------###


def get_equivalents(ref_seq_to_here, ref_offset, query_to_here, query_offset):
    '''
    Check to make sure both of these are valid residues and not gaps.
    If they are valid, it returns a tuple of the equivalent residues positions,
    adjusted by the appropriate offset. Returned numbers are not zero-indexed.

    Return `None` if one of them matches gap indicator.
    '''
    current_chars = [ref_seq_to_here[-1],query_to_here[-1]]
    if any(x == gap_indicator for x in current_chars):
        return None
    else:
        return (len(ref_seq_to_here.ungap(gap_indicator)) + ref_offset, 
            len(query_to_here.ungap(gap_indicator)) + query_offset)


def end_of_match_block(current_set,next_set):
    '''
    Takes current matching residues and the next set of matching residues 
    and determines if there is  a 'break' next because element for either chain is 
    not the next consecutive position for each relative current.

    Return True if that match block is the end of that accumulating block of 
    positions for both chains or false if accumulation of block of matches is to 
    continue.
    '''
    # signal end if any element is not next increment
    if current_set[0] + 1 == next_set[0] and current_set[1] + 1 == next_set[1]:
        return False
    else:
        return True

def range_extract(lst):
    'Yield 2-tuple ranges or 1-tuple single elements from list of increasing ints'
    'modified from  https://www.rosettacode.org/wiki/Range_extraction#Python'
    lenlst = len(lst)
    i = 0
    while i< lenlst:
        low = lst[i]
        while i <lenlst-1 and lst[i]+1 == lst[i+1]: i +=1
        hi = lst[i]
        if hi - low >= 1:    #<---MAIN DIFFERENCE from https://www.rosettacode.org/wiki/Range_extraction#Python
            yield (low, hi)
        else:
            yield (low,)
        i += 1

def make_intervals(lizt_o_lists):
    '''
    Takes a list of list and works through each making ranges. Reason it has to 
    be a list of lists is so respects breaks in matched sets correspond to 
    breaks in sequence of either chain that occurs in the reference-query 
    matching.
    The final result is a list of intervals or residues that is returned
    '''
    intervals = []
    for l in lizt_o_lists:
        current_intervals_and_residues = list(range_extract(l))
        intervals.extend(current_intervals_and_residues)
    return intervals
    # flatten the separate lists of tuples to a single list; turns out it was 
    # easier than the example at https://stackoverflow.com/a/952952/8508004
    flat_list = [item for item in intervals]
    return flat_list
    

###--------------------------END OF HELPER FUNCTIONS---------------------------###
###--------------------------END OF HELPER FUNCTIONS---------------------------###



#*******************************************************************************
###------------------------'main' function of script---------------------------##

def MSA_to_corresponding_residue_numbers(
    alignment, ref_id, supplied_start_pos = None, return_dfs = False):
    '''
    Main function of script. 
    It will take an alignment text file (CLUSTAL form) determine the determines
    the corresponding residue numbers for aligned residues for a specified 
    reference sequence in the alignment and the pairing of it with every other 
    aligned sequence. The reference sequence is specified by `ref_id`.
    Presumably the reference sequence is one that has been solved structurally 
    and you want to know the exact residues that match so you can use the 
    backbone (a.k.a., 'spine') of the chain as a guide to fit or highlight or 
    otherwise process the other chains in the alignment in some manner.

    The default is to assume that all sequences in the multiple sequence 
    alignmentbeing with residue number one. `supplied_start_pos` allows for a 
    list to be supplied with the first position for each sequence if at least 
    one differs from starting at 1. The list must be the same length as the 
    number of sequences in the multiple sequence alignment and so supply `1` for
    each of those that begin at the first residue.

    While the sequences don't need to start at residue number one if 
    `supplied_start_pos` is used, they do need to be consecutive once they 
    start in order for numbers extracted to make any sense. In other words,
    no residues can be skipped internally. You can use my script 
    `check_seq_frag_in_MSAclustal_intact_viaFASTA.py` to aid in validating this.

    Sequences at the end of the sequence need not be present in the multiple
    sequence alignment; howerver, if that is the case, no equivalent residues
    can be supplied for the missing end. This may be reasonable if that part is 
    of no concern for downstream applications.

    Note that all numbers returned by the script will be in common terms where
    the fist residue in a sequence is numbered ONE. (In other words, the 
    zero-indexing behind-the-scenes will not be visible to user.)

    If `return_dfs` is supplied as `True`, it will return the `ref_id` and a
    dictionary of the dataframes of the matched residue blocks (start and end 
    for each). This is meant to be used when calling from a Jupyter notebook or 
    IPython.
    The script saves a tab-separated file of the matched residue blocks for each 
    reference sequence-query pairing.
    '''
    # feedback about important option
    sys.stderr.write("\n**NOTE: gap indicator in this script is currently set "
        "to '{}'. If\nthat does not match what provided alignment uses to "
        "indicate gaps,\nplease change the setting within the script code "
        "under\n'USER ADJUSTABLE VALUES' around line 130 "
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



    # Get reference sequence user specified from alignment and leave the others
    # as a those to compare, in turn
    #---------------------------------------------------------------------------
    msa_ids = []
    query_seqs = {} # make a dictionary for collecting all the sequences that 
    # aren't the reference; ids as keys. Sequences are values.
    order_of_nonref_records = {}
    for indx,record in enumerate(alignment):
        msa_ids.append(record.id)
        if record.id == ref_id:
            ref_seq = record.seq
            sys.stderr.write(
                "{} reference sequence collected from alignment...".format(ref_id))
            # track position of reference sequence among the records so if 
            # any start position information will get associated.
            ref_index_among_records = indx
        else:
            query_seqs[record.id] = record.seq
            # track position for each non-reference id so can be associated with
            # any provided start position information later
            order_of_nonref_records[record.id] = indx
    try: 
        if ref_seq != None:
            pass
    except UnboundLocalError:
        sys.stderr.write("\n***    ERROR    ***   ERROR   ***\nCAUSE OF ERROR:")
        sys.stderr.write("\nNo matches to provided identifier were found in "
            "the multiple sequence alignment.\nUser specified id of "
            "reference: {}\n"
            "Identifiers in the provided alignment: {}.\n".format(id_, repr(msa_ids)))
        sys.stderr.write("***  END OF ERROR DESCRIPTION ***\n")
        sys.exit(1)



    # Deal with start positions
    #---------------------------------------------------------------------------
    # build in assert to check `len(start_pos)` equals number of `msa_ids`
    if supplied_start_pos:
        start_pos = supplied_start_pos
    else:
        start_pos = [1] * len(msa_ids)
    assert len(start_pos) == len(msa_ids), ("The number of supplied start "
        "positions must equal the number of sequences in the provided "
        "alignment.\nThere were {} supplied start posistions and {} sequences "
        "in the alignment.\nNote: if you provide even one start position, you "
        "must provide one for each sequence in the alignment."
        .format(len(start_pos),len(msa_ids)))
    # associate reference sequence start position with reference sequence
    offset_for_reference = start_pos[ref_index_among_records] - 1 # Using ` - 1`
    # so `offset_for_reference` results in zero if there is no offset
    # associate non-reference sequence start positions with related
    offsets_by_id = {}
    for id_ in order_of_nonref_records:
        offsets_by_id[id_] = start_pos[order_of_nonref_records[id_]] - 1



    

    # Go through each query set and compare to reference resulting in a pair of 
    # matched residues for each matching position of each reference-query 
    # pairing. Reference is always first.
    #---------------------------------------------------------------------------
    query_matches_by_id = {} # a dictionary where id of query alignment is the 
    # key and values are a lit of tuples the matched residues with reference 
    # sequence position listed first and corresponding query position listed second.
    for id_ in query_seqs:
        # Using the squence associated with the current id as the current 
        # query, step through each character listed in reference sequence and 
        # compare to corresponding character for current query sequence. If 
        # both aren't gap indicators, then count as matched.
        matched_positions = []
        for indx,charct in enumerate(ref_seq):
            query_matches = get_equivalents(ref_seq[:indx + 1], 
                offset_for_reference, 
                query_seqs[id_][:indx + 1], offsets_by_id[id_])
            if query_matches != None:
                matched_positions.append(query_matches)

        query_matches_by_id[id_] = matched_positions
        # feedback
        sys.stderr.write("query compared...")



    # For each id, go through the matched individual residues and make 
    # 2 item-tuple ranges or 1 item-tuple elements for both reference and query 
    # chains. (NOTE: 1 item-tuple element == single matched residue occurences) 
    #---------------------------------------------------------------------------
    # For each pairing, store these as a list with reference chain list first 
    # and query chain list second.
    # What makes this trickier than the solution I used from 
    # https://www.rosettacode.org/wiki/Range_extraction#Python
    # in the notebook entitled 
    # 'Using Biopython's PDB module to list resolved residues and construct fit commands.ipynb' 
    # in 'cl_demo-binder' @ github is that here I cannot just make one list and 
    # package consecutive residues. It has to be done without respect to both 
    # sequences so any condensed ranges (intervals) match up. For example, if 
    # there is a gap between say the 
    # query sequence while the say reference sequence is contiguous and I make 
    # a range using it, they won't come out with matching numbers of ranges. The 
    # reference sequence will be represented as one range and the query sequence 
    # will be represented as two and you cannot make a fit command with 
    # unequal numbers of backbone atoms.
    # I think the solution is to make multiple lists for each sequence using 
    # the other sequence to guide where the lists should end and then 
    # process all those separate into ranges and then stitch them together for 
    # the appropriate chain at the end. So that they match up with equal numbers 
    # of backbone atoms represented.

    # To start, go through for each id nd put the reference and query residues 
    # in lists, with one match set in each list. Lists divide into new lists 
    # at each non-consecutive break in either chain. List of lists for 
    # reference chain is first in tuple and reference is second.
    # Then convert these lists of lists to range intervals where possible.
    lists_of_ref_n_query_residues_pairings_by_id = {}
    lists_of_ref_n_query_residues_block_pairings_by_id = {}
    for id_, match_pairings in query_matches_by_id.items():
        ref_n_query_lists = ([],[])
        ref_res_list = []
        query_res_list = []
        for indx,match_pair in enumerate(match_pairings):
            ref_res, query_res = match_pair
            ref_res_list.append(ref_res)
            query_res_list.append(query_res)
            if indx == (len(match_pairings) - 1):
                ref_n_query_lists[0].append(ref_res_list)
                ref_n_query_lists[1].append(query_res_list)
            elif end_of_match_block(match_pair,match_pairings[indx+1]):
                ref_n_query_lists[0].append(ref_res_list)
                ref_n_query_lists[1].append(query_res_list)
                ref_res_list = []
                query_res_list = []


        # add the 2 item - tuple of lists of lists for each sequence to the 
        # dictionary. First item in tuple will be list of lists for reference
        # sequence and the second itme in tuple will be list of lists for the
        # query corresponding to that id (key)
        lists_of_ref_n_query_residues_pairings_by_id[id_] = ref_n_query_lists

        # Collapse these to blocks of start and end intervals where possible, 
        # i.e., contiguous
        ref_collapsed = make_intervals(ref_n_query_lists[0])
        query_collapsed = make_intervals(ref_n_query_lists[1])

        lists_of_ref_n_query_residues_block_pairings_by_id[id_] = (
            ref_collapsed, query_collapsed)

    sys.stderr.write("\nmade match pairings into intervals of start and end "
        "where possible...\n")



    # Convert the values `lists_of_ref_n_query_residues_block_pairings_by_id` 
    # to dataframe for easy options of how to proceed. Having it as that would 
    # let me pivot any number of ways. For example, I can easily store as a 
    # tabular text file or return in dataframe form for further use
    #---------------------------------------------------------------------------
    residue_block_pairing_dfs_by_id = {}
    for id_, l_o_paired_tuples in (
        lists_of_ref_n_query_residues_block_pairings_by_id.items()):
        # make dataframe from `l_o_paired_tuples` (meaning 
        # 'list of paired tuples'). 
        # Example of a `l_o_paired_tuples`:
        # ([(1, 76), (82, 84), (85, 99), (100, 144), (145, 163), (164, 178), 
        # (179, 214), (217, 259), (263, 317), (320, 652), (653, 667), 
        # (668, 698), (699, 756), (757, 815), (822, 825)], [(1, 76), (77, 79), 
        # (84, 98), (110, 154), (160, 178), (182, 196), (233, 268), (269, 311), 
        # (312, 366), (367, 699), (709, 723), (728, 758), (760, 817), 
        # (821, 879), (880, 883)])
        # Want to cast first item of 2 item-tuple element to the list for 
        # reference sequence and the other to the query (a.k.a.,current id)
        # Each tuple in both lists will be start and end and start and end of a 
        # row so that the matched pairs are kept. In other words the columns 
        # will be `'ref_start', 'ref_end', 'id_start', 'id_end'`.
        # For 1 item-tuple elements (i.e, single matched residue occurences) 
        # want to cast the same number to both `start` and `end`
        rows_parsed_out = []
        ref_list = l_o_paired_tuples[0]
        query_list = l_o_paired_tuples[1]
        assert len(ref_list) == len(query_list), ("Matched pairings should "
            "mean lists are same length.")
        for indx,tup in enumerate(ref_list):
            if len(tup) == 1:
                rows_parsed_out.append((tup[0],tup[0],
                    query_list[indx][0],query_list[indx][0]))
            else:
                rows_parsed_out.append((tup[0],tup[1],
                    query_list[indx][0],query_list[indx][1]))
        '''labels = ['ref_seq_start',
                'ref_seq_end','query_seq_start','query_seq_end']'''
        labels = ['{}_start'.format(ref_id),'{}_end'.format(ref_id),
                    '{}_start'.format(id_),'{}_end'.format(id_)]
        df = pd.DataFrame.from_records(rows_parsed_out, columns=labels)
        residue_block_pairing_dfs_by_id[id_] = df
        # Resulting dataframe example using data from example above (before 
        # column ids adjusted to be more specific):
        '''
                    ref_seq_start  ref_seq_end  query_seq_start  query_seq_end
        0               1           76                1             76
        1              82           84               77             79
        2              85           99               84             98
        3             100          144              110            154
        4             145          163              160            178
        5             164          178              182            196
        6             179          214              233            268
        7             217          259              269            311
        8             263          317              312            366
        9             320          652              367            699
        10            653          667              709            723
        11            668          698              728            758
        12            699          756              760            817
        13            757          815              821            879
        14            822          825              880            883
        '''
    if len(residue_block_pairing_dfs_by_id) == 1:
        sys.stderr.write("\nmade dataframe of matched start and ends...")
    else:
        sys.stderr.write("\nmade dataframes of matched start and ends...")
    sys.stderr.write("DONE.\n")






    # Return reference id and a dictionary of the lists of the dataframes with 
    # ids as keys if called with `return_dfs = True`. Otherwise, consider 
    # called from command line and output the dataframes as tab separated values 
    # with a start and end for each query and reference chain. For 1 item-tuple 
    # elements (i.e., sing residue pairs), there will be the same number for 
    # start and end columns. For each id, there will be a tab-separated values 
    # file.
    #---------------------------------------------------------------------------
    if return_dfs:
        return ref_id, residue_block_pairing_dfs_by_id
        sys.stderr.write("Reference identifier and a dictionary of dataframes "
            "have been returned.\n")
    else:
        sys.stderr.write("\n")
        for id_, df in residue_block_pairing_dfs_by_id.items():
            tsv_nom = '{}_residues_matched_to_{}.tsv'.format(ref_id,id_)
            df.to_csv(tsv_nom, sep='\t',index = False)
            sys.stderr.write("File '{}' saved.\n".format(tsv_nom))





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
    #kwargs['suffix_for_saving'] = suffix_for_saving
    kwargs['supplied_start_pos'] = supplied_start_pos
    MSA_to_corresponding_residue_numbers(alignment, ref_id, **kwargs)
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
    # For 'nargs', see https://stackoverflow.com/a/15753721/8508004 for where
    # use of `nargs='*'` allows specifying zero or more in the list.
    ###-----------------for parsing command line arguments-----------------------###
    import argparse
    parser = argparse.ArgumentParser(prog=
        'MSA_to_corresponding_residue_numbers.py',
        description="MSA_to_corresponding_residue_numbers.py \
        takes akes a multiple sequence alignment in Clustal format and \
        determines the corresponding residue numbers for aligned residues for \
        a specified reference sequence in the alignment and the pairing of it \
        with every other aligned sequence in the multiple sequence alignment. \
        The reference sequence is specified by the identifier at the start of \
        each line in the sequence \
        blocks.\
        **** Script by Wayne Decatur   \
        (fomightez @ github) ***")

    parser.add_argument("align_file", help="Name of file of alignment text \
        file (CLUSTAL format).\
        ", metavar="ALIGNMENT_FILE")

    parser.add_argument("ref_id", help="Identifier that matches the sequence \
        in the MSA to use as a reference to get matches to other sequences.\
        Presumably, this sequence has a corresponding structural model \
        available.\
        ", metavar="REF_ID")

    parser.add_argument('integers', metavar='N', type=int, nargs='*',
        help="OPTIONAL IN MANY CASES: Provide an integer \
        specifying the start position for each sequence in the provided \
        alignment, if at \
        least one of them isn't beginning at the first residue in the \
        sequence. Seperate each integer with a space, such as `5 1 9 8` \
        without tickmarks. IF ALL BEGIN AT ONE, \
        NO REASON TO INCLUDE SUCH A LIST.")

    '''parser.add_argument("fasta_file", help="Name of sequence file to check \
        against (FASTA format).\
        ", metavar="FASTA_FILE")


    parser.add_argument('-ds', '--descr_source', action='store', type=str, 
    default= None, help="OPTIONAL: Provide FASTA file with same \
    ids to use as source of alignment to use as source of descriptions. The \
    typical file is output as `*_extracted_ungapped.fa` files from the \
    `extract_regions_from_clustal_alignment.py` script\
    .")
    '''





    #I would also like trigger help to display if no arguments provided because 
    # need at least one for url
    if len(sys.argv)==1:    #from http://stackoverflow.com/questions/4042452/display-help-message-with-python-argparse-when-script-is-called-without-any-argu
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    alignment = args.align_file
    ref_id = args.ref_id
    if args.integers:
        supplied_start_pos = args.integers
    else:
        supplied_start_pos = None


    main()

#*******************************************************************************
###-***********************END MAIN PORTION OF SCRIPT***********************-###
#*******************************************************************************
