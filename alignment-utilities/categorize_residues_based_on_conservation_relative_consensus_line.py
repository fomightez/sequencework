#!/usr/bin/env python
# categorize_residues_based_on_conservation_relative_consensus_line.py
__author__ = "Wayne Decatur" #fomightez on GitHub
__license__ = "MIT"
__version__ = "0.1.0"


# categorize_residues_based_on_conservation_relative_consensus_line.py by 
# Wayne Decatur
# ver 0.1
#
#*******************************************************************************
# Verified compatible with both Python 2.7 and Python 3.6; written initially in 
# Python 3. 
#
#
# PURPOSE: Takes an multiple sequence alignment (in CLUSTAL format) that has a 
# consensus line, say from MUSCLE at 
# https://www.ebi.ac.uk/Tools/msa/muscle/, and for a specific sequence in the 
# alignment categorizes the residues that are identical, strongly, similar, or 
# weakly similar in the alignment. See 'What do the consensus symbols mean in 
# the alignment?' at 
# https://www.ebi.ac.uk/Tools/msa/clustalw2/help/faq.html#23 for indicators. 
# Plus also categorizes unconserved, while at it.
# Importantly, residue positions in the results are in common terms where the 
# first residue is number one. (I.E., zero-indexing of Python is hidden 
# behind-the-scenes.)
#
#
# Assumes the CLUSTAL alignment is provided with a header line as it would 
# come from EMBl-EBI's Muscle. Doesn't check it matches anything. Just assumes
# nothing in the first line is of interest.
# 
# Note, because white space is critical for the consensus symbols line, it is 
# best to save the alignment file directly from EMBL-EBI and use that file as 
# input for this script, rather than doing copy-paste from the site. For 
# example, via copy-paste it may be easy to miss the spaces on the last line of 
# the consensus symbols line in the case of two sequences that mismatch for the 
# span of the entire last row of an alignment.
#
#
# Written to run from command line or with the main function imported into a 
# Jupyter notebook. When using main function imported into a notebook it can
# pass back dataframe(s) with the results. See the example commands and demo.
#
#
#
#
#
# Dependencies beyond the mostly standard libraries/modules:
#
#
# VERSION HISTORY:
# v.0.1. very basic working version
#
#
# To do:
# - 
#
#
#
# TO RUN:
# Examples,
# Enter on the command line of your terminal, the line
#-----------------------------------
# python categorize_residues_based_on_conservation_relative_consensus_line.py alignment.txt seq_id
#-----------------------------------
#
# Issue `categorize_residues_based_on_conservation_relative_consensus_line.py -h` for details.
# 
#
#
# To use this after importing/pasting or loading into a cell in a Jupyter 
# notebook, specify at least the alignment file and id of sequence to analyze 
# in the call to the main function similar to below:
# categorize_residues_based_on_conservation_relative_consensus_line("alignment.clw","VPH1")
# 
#-or-
# If `return_panel_data` is `True` and `output_separate` is `True` then it will
# return a dictionary of the dataframes with the categories as keys, one dataframe
# with the positions for each category. For that all to the main function 
# similar to below:
# categorize_residues_based_on_conservation_relative_consensus_line("alignment.clw","VPH1",output_separate = True,return_panel_data=True)
#
# There are other options, see the demo for illustration.
# 
#
'''
CURRENT ACTUAL CODE FOR RUNNING/TESTING IN A NOTEBOOK WHEN IMPORTED/LOADED OR 
PASTED IN ANOTHER CELL:
categorize_residues_based_on_conservation_relative_consensus_line("alignment.clw","VPH1",output_separate=True,return_panel_data=True)
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

gap_indicator = "-" # Possible tip: Since '.' is the other popular one I know   
# of, `category_symbols_meaning` (below) may need changing also.

# from https://www.ebi.ac.uk/Tools/msa/clustalw2/help/faq.html#23
category_symbols_meaning = {
    '*':'identical',
    ':':'strongly_similar',
    '.':'weakly_similar',
}

#
#*******************************************************************************
#**********************END USER ADJUSTABLE VARIABLES****************************


















#*******************************************************************************
#*******************************************************************************
###DO NOT EDIT BELOW HERE - ENTER VALUES ABOVE###

import sys
import os
from collections import defaultdict
import pandas as pd



###--------------------------HELPER FUNCTIONS--------------------------------###



def update_caetgorization(
    position, consensus_symbol, categorization_dict):
    '''
    Takes a current position (integer) and places that in the correct category 
    based on the provided consensus symbol. 

    Residue positions in the returned dictionary are in common terms where first
    residue is number one.

    Returns an updated dictionary.
    '''
    if consensus_symbol in category_symbols_meaning.keys():
        categorization_dict[category_symbols_meaning[consensus_symbol]].append(
            position+1) #`+1` so numbers will match common numbering system 
            # where first residue is numbered one.
    else:
        categorization_dict['not_conserved'].append(position+1) #`+1` so numbers
            # will match common numbering system where first residue is 
            # numbered one.
    return categorization_dict

def to_comma_sep_string (list_of_ints):
    '''
    Takes a list of integers and returns a string of those values separated by 
    commas.
    '''
    return ",".join([str(i) for i in list_of_ints])


###-------------------------END OF HELPER FUNCTIONS--------------------------###
###-------------------------END OF HELPER FUNCTIONS--------------------------###



#*******************************************************************************
###-----------------------'main' function of script---------------------------##

def categorize_residues_based_on_conservation_relative_consensus_line(
    alignment, id_, output_single = True, output_separate = False, 
    output_table_text = True, output_pickled_panel_data = False, 
    return_panel_data = False):
    '''
    Main function of script. 
    It will take a multiple sequence alignment (in CLUSTAL format) that has a 
    consensus line with symbols indicating conservation, say from MUSCLE and 
    for the user specified sequence it categorizes each residue as identical, 
    strongly similar, weakly similar, or unconserved.

    Assumes the CLUSTAL alignment is provided with a header line as it would 
    come from EMBl-EBI's Muscle. Doesn't check it matches anything. Just assumes
    nothing in the first line is of interest.

    Returns a dataframe, or dictionary of dataframes with positions for each 
    category with conservation categories as keys, and/or makes a tab-separate 
    tabular text file of the results.

    Residue positions in the resulting lists are in common terms where first
    residue is number one.
    '''
    # feedback about important option
    sys.stderr.write("\n**NOTE: gap indicator in this script is currently set "
        "to '{}'. If\nthat does not match what provided alignment uses to "
        "indicate gaps,\nplease change the setting within the script code "
        "under\n'USER ADJUSTABLE VALUES' around line 110 "
        "(give or take a few).**\n".format(gap_indicator))



    # Bring in the necessary data:
    #---------------------------------------------------------------------------
    # Reading in of the alignment file is adapted from 
    # `add_actual_position_numbering_to_sequence_in_top_line_of_clustal_output.py`
    # and `calculate_cons_for_clustal_nucleic.py`
    # and `score_sequences_in_clustal_msa_favoring_top_line.py` and 
    # `score_sequences_in_clustal_msa.py`
    # because Biopython seems to not put the consensus line in any form I 
    # could access. (I looked around quite a bit and didn't find anything, but 
    # that doesn't mean it doesn't actually exist in Bipython, just that it 
    # isn't obvious.) Thankfully I had previously parsed alignments manual via those
    # script and could build on it. (Note this is as opposed to using Biopython
    # `AlignIO.read()`which worked for sequences in 
    # `extract_regions_from_clustal_alignment.py`.)

    try:
        with open(alignment, 'r') as the_file:
            alignment = the_file.read()
    except (TypeError,OSError,IOError) as e:
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
        # Note "FileNotFoundError is a subclass of OSError"(https://stackoverflow.com/a/28633573/8508004)

    # feedback
    sys.stderr.write("Alignment read...")

    

    # Process
    #---------------------------------------------------------------------------

    # Identify the first identifier in each alignment block
    first_words = []
    for line in alignment.split("\n"):
        if line:
            first_word = line.split()[0]
            if first_word in first_words:
                first_id = first_word
                break
            else:
                first_words.append(first_word)
    # feedback
    sys.stderr.write(
        "top line identifier determined as '{}'...\n".format(first_id))

    # Now that have first identifier, get column where sequences begin and end.
    for line in alignment.split("\n"):
         if line.startswith(first_id):
            seq = line.split(first_id)[1].strip()
            start_col = line.index(seq)
            end_col = start_col + len(seq)
            #print (start_col) # FOR DEBUGGING ONLY
            #print (end_col)   # FOR DEBUGGING ONLY
            break


    identifier_list = first_words[1:-1] # slice [1:-:1] because don't want the 
    # header line (first word of firstline) or the consensus line text 
    # (last first_word) among the indentifiers
    #print (identifier_list )  # FOR DEBUGGING ONLY

    assert id_ in identifier_list, (
        "The user specified sequence identifier, '{}', is not among the "
        "extracted\nidentifiers. Recognized ids:'{}'.".format(
        id_,", ".join(identifier_list)))

    # Go through each line of the sub-blocks, assuming a blank line between each
    # sub-block. The blank line would be right before the first_id. Collect 
    # each gapped string to grow into one string based on any previous one.
    # In other words, make a string for each identifier and the consensus.
    # keep growing strings in a dictionary with identifiers as consensus for 
    # those with identifiers and `consensus` for the consensus line.
    growing_seqs_and_consensus = defaultdict(str)
    current_position = 0
    in_subblock = False
    last_identifier_insublock_encountered = False
    for line in alignment.split("\n"):
        if len(line) > 2: # removed `.strip` & made bigger than 0 so if sequence
        #ends in a gap it will still have correnctly sized consensus
            if line.startswith(first_id) and not in_subblock:
                #sub-block begins here
                in_subblock = True
            if in_subblock:
                try:
                    first_wd = line.split()[0]
                except IndexError: # This is because if ends in gaps wasn't 
                #getting last line until I removed strip above but then this 
                #failed
                    first_wd = " "
                if first_wd in identifier_list:
                    #print (line[start_col:end_col])       # FOR DEBUGGING ONLY
                    #print(type(line[start_col:end_col]))  # FOR DEBUGGING ONLY
                    growing_seqs_and_consensus[first_wd] += line[start_col:end_col]
                if line.startswith(identifier_list[-1]):
                    last_identifier_insublock_encountered = True
                    continue
                elif last_identifier_insublock_encountered:
                    growing_seqs_and_consensus['consensus'] += line[start_col:end_col]
                    # reset for moving to the next subblock
                    in_subblock = False
                    last_identifier_insublock_encountered = False
    #print(growing_seqs_and_consensus)  # FOR DEBUGGING ONLY

    # Now go through the sequence corresponding to the designated id and remove 
    # the gaps from the sequence and delete the corresponding column 
    # in the consensus line. REmoving gaps so can count from first residue so
    # cetegorized residues match position in chain. 
    # Initially, I mistakenly thought I could remove gaps and categorize at same 
    # time, however...
    # Cannot combine because found very efficient to count actual residues 
    # when sequences get really large and so easier to just remove gaps first. 
    # Plus need to remove first so iterate over it and adjustments and 
    # assignments will be easier.
    sys.stderr.write(
        "...Processing information about sequence with id "
        "of '{}'.\n".format(id_))
    gapped_seq_o_interest = growing_seqs_and_consensus[id_]
    seq_o_interest = ''
    gapped_consensus = growing_seqs_and_consensus['consensus']
    consensus = ''
    assert len(gapped_seq_o_interest) == len(gapped_consensus), (
    "The sequence with the id specified and the consensus must be the " 
    "same size.")
    for indx,char in enumerate(gapped_seq_o_interest):
        if char != gap_indicator:
            seq_o_interest += char

            # Also only collect corresponding columns of consensus symbols to 
            # keep in register relative the sequence.
            consensus += gapped_consensus[indx]

    assert len(seq_o_interest) == len(consensus), (
    "The sequence and the consensus, both with the gaps removed, must be the " 
    "same size.")

    # Then iterate over the ungapped sequence and categorize conservation.
    categorization = defaultdict(list) # make lists for each category for 
    # making into dataframe at end
    for indx,char in enumerate(seq_o_interest):
        categorization = update_caetgorization(
            indx, consensus[indx], categorization)
    #print(categorization)                      # FOR DEBUGGING ONLY
    if output_table_text:
        # Want positions separated by comma for possible tabular text output
        same_cat_dict_for_tbl = {k: ",".join(
            [str(i) for i in v]) for (k, v) in categorization.items()}
        # print(same_cat_dict_for_tbl)                   # FOR DEBUGGING ONLY

    # Make the results into a dataframe
    categorized_residue_positions_df = pd.DataFrame(
        list(categorization.items()),columns = ['category','residue_positions'])
    # That works but I want categories with most conserved as top line and 
    # `not_conserved` on bottom
    # Because I think the dictionary will have these as arbitrary orders I
    # cannot simply base order on what I saw in development. More robust would
    # be to extract what `new_index` order should be
    #print(categorized_residue_positions_df)     # FOR DEBUGGING ONLY
    default_indx = {}
    for i, row in categorized_residue_positions_df.iterrows():
        default_indx[row.category] = i
    new_index = ([default_indx['identical'],
        default_indx['strongly_similar'],
        default_indx['weakly_similar'],
        default_indx['not_conserved'],])
    categorized_residue_positions_df = categorized_residue_positions_df.reindex(
        new_index)# based on https://stackoverflow.com/a/30010004/8508004
    categorized_residue_positions_df = categorized_residue_positions_df.reset_index(
        drop=True)
    #print(categorized_residue_positions_df)              # FOR DEBUGGING ONLY

    # Make results into dataframes if need multiple
    if output_separate:
        sep_dfs = {}
        for category, pos_list in categorization.items():
            df = pd.DataFrame({category: pos_list})
            df = df.rename(columns={category:'{}_residue_pos'.format(category)})
            sep_dfs[category] = df
        #make a list of dataframes ordered most to least conserved for 
        # potentially returning
        sep_dfs_ordered = ([sep_dfs['identical'],
            sep_dfs['strongly_similar'],sep_dfs['weakly_similar'],
            sep_dfs['not_conserved']])
        #print(sep_dfs)                               # FOR DEBUGGING ONLY
        #print(sep_dfs_ordered)                       # FOR DEBUGGING ONLY


    
    # Saving and reporting or handling returning results
    #---------------------------------------------------------------------------


    if output_table_text and output_single:
        # Want positions separated by comma for possible tabular text output
        same_cat_df_for_tbl_df = categorized_residue_positions_df.copy() #first copy df with `deep=True` so subsequent changes are 
        # not reflected in original
        # Now change so `residue_positions` column is the listed values now as 
        #string separated by commas
        same_cat_df_for_tbl_df['residue_positions'] = (
            same_cat_df_for_tbl_df['residue_positions'].apply(
            to_comma_sep_string))
        tsv_file_name = "categorized_consv_{}_residues.tsv".format(id_)
        same_cat_df_for_tbl_df.to_csv(tsv_file_name, sep='\t',index = False)
        sys.stderr.write(
            "\nData a saved in tabular text form (tab-separated form)"
            " as '{}'.\n".format(tsv_file_name))
        #print(same_cat_df_for_tbl_df)                      # FOR DEBUGGING ONLY

    if output_table_text and output_separate:
        sys.stderr.write("\nData saved in tabular text form "
            "(tab-separated form) as:")
        for category, df in sep_dfs.items():
            tsv_file_name = "{}_{}_residues.tsv".format(category,id_)
            df.to_csv(tsv_file_name, sep='\t',index = False)
            sys.stderr.write(
                "\n       '{}'.\n".format(tsv_file_name))


    

    if output_pickled_panel_data and output_single:
        single_pickled_file_name="categorized_consv_{}_residues.pkl".format(id_)
        categorized_residue_positions_df.to_pickle(single_pickled_file_name)
        sys.stderr.write(
            "\nData saved in pickled dataframe form"
            " as '{}'.\n".format(single_pickled_file_name))

    if output_pickled_panel_data and output_separate:
        sys.stderr.write("\nData saved in pickled dataframe form as:")
        for category, df in sep_dfs.items():
            pickled_file_name = "{}_{}_residues.pkl".format(category,id_)
            df.to_pickle(pickled_file_name)
            sys.stderr.write(
                "\n       '{}'.\n".format(pickled_file_name))



    #Handle potential returns last so calling function allows specifying 
    # multiple forms of output, plus returning dataframe(s).
    # (Because `output_separate` defaults to false, & output_single defaults to
    # True, put output_separate first so will return before testing single when
    # separate is specified.)
    if return_panel_data:
        if output_separate:
            #return sep_dfs_ordered # was considering returning a dataframe
            # for each category ordered by conservation from highest to lowest
            # but makes more sense to return dictionary with that info in keys.
            # Leaving that here in case want to switch at some point.

            sys.stderr.write("\nReturning a a dictionary of {} dataframes "
                ", one dataframe for each category. The dictionary keys\nare "
                "the categories:{}.".format(
                    len(sep_dfs),",".join([repr(x) for x in sep_dfs.keys()])))
            return sep_dfs # the dictionary of dataframes with categories as
            # keys. There is a dataframe with the residue positions for each
            # category as rows. The `repr()` trick to get the quotes for each
            # to stay is from https://stackoverflow.com/a/13207713/8508004.
        if output_single:
            sys.stderr.write("\nReturning a single dataframe of residue "
                "positions, one list per category.")
            return categorized_residue_positions_df


    

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
    kwargs['output_single'] = output_single
    kwargs['output_separate'] = output_separate
    kwargs['output_pickled_panel_data'] = output_pickled_panel_data
    kwargs['output_table_text'] = output_table_text
    kwargs['output_pickled_panel_data'] = output_pickled_panel_data
    categorize_residues_based_on_conservation_relative_consensus_line(
        alignment,id_,**kwargs)
    # using https://www.saltycrane.com/blog/2008/01/how-to-use-args-and-kwargs-in-python/#calling-a-function
    # to build keyword arguments to pass to the function above
    # (see https://stackoverflow.com/a/28986876/8508004 and
    # https://stackoverflow.com/a/1496355/8508004 
    # (maybe https://stackoverflow.com/a/7437238/8508004 might help too) for 
    # related help). Makes it easy to add more as work out which are needed.





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
        'categorize_residues_based_on_conservation_relative_consensus_line.py',
        description="categorize_residues_based_on_conservation_relative_consensus_line.py \
         takes an multiple sequence alignment (in CLUSTAL format) that has a \
         consensus line, say from MUSCLE at \
         https://www.ebi.ac.uk/Tools/msa/muscle/, and for a specific sequence \
         in the alignment categorizes the residues that are identical, strongly\
         , similar, weakly similar, or unconserved in the alignment.\
        **** Script by Wayne Decatur   \
        (fomightez @ github) ***")

    parser.add_argument("align_file", help="Name of file of alignmnet text \
        file to use to categorize residues.\
        ", metavar="ALIGNMENT_FILE")

    parser.add_argument("id", help="Identifier that has the residues to \
        categorize relative consensus symbols.\
        ", metavar="ID")
    parser.add_argument("-og", "--output_grouping", type=str, 
        default= "single", choices=["single", "separate"],
        help="OPTIONAL: Specify grouping of output with this option. Choose \
        `-og single` for one table or dataframe for all categories. Or choose \
        `-og separate` for a separate table or dataframe for each category. \
        If this option is not specified, {} will be used.".format("single"))
    parser.add_argument("-ot", "--output_type", type=str, 
        default= "tabular_text", choices=["tabular_text", "panel_data"],
        help="OPTIONAL: Specify output type with this option. Choose \
        `-ot tabular_text` for tabular data in a tab-separated text file. \
        Choose `-ot panel_data` for dataframe or dataframes. If this option \
        is not specified, {} will be used. For either, whether output files \
        singular or multiple depends on `og` option setting.".format(
        "tabular_text"))



    #I would also like trigger help to display if no arguments provided because 
    # need at least one for url
    if len(sys.argv)==1:    #from http://stackoverflow.com/questions/4042452/display-help-message-with-python-argparse-when-script-is-called-without-any-argu
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    alignment = args.align_file
    id_ = args.id


    #grouping
    output_single = True
    output_separate = False
    if args.output_grouping == "single":
        output_single = True
    else:
        output_separate = True
        output_single = False

    #text or pickled dataframe(s)
    output_table_text = True
    output_pickled_panel_data = False
    if args.output_type == "tabular_text":
        output_table_text = True
    else:
        output_pickled_panel_data = True
        output_table_text = False



    main()

#*******************************************************************************
###-***********************END MAIN PORTION OF SCRIPT***********************-###
#*******************************************************************************
