#!/usr/bin/env python
# summarize_all_nts_even_ambiguous_present_in_FASTA.py
__author__ = "Wayne Decatur" #fomightez on GitHub
__license__ = "MIT"
__version__ = "0.1.0"


# summarize_all_nts_even_ambiguous_present_in_FASTA.py by 
# Wayne Decatur
# ver 0.1
#
#*******************************************************************************
# Verified compatible with both Python 2.7 and Python 3.6; written initially in 
# Python 3. 
#
#
# PURPOSE: Takes a sequence file (FASTA-format) and summarizes the nts in the 
# sequence(s). Assumes multi-FASTA, but single sequence entry is fine, too. 
#
# When running on the command line, it will print out a summary table of
# counts of nucleotides and other character in each sequence and totals. When 
# calling the main function it will, by default, return a dataframe with this 
# information.
#
# Only valid for DNA sequences; script has no step checking for data type, and 
# so you are responsible for verifying appropriate input.
#
# (incorporates quantification code adapted from my 
# 'Assessing_ambiguous_nts..' series of notebooks, which first was used to make the summarizing/accounting part of `replace_unusual_nts_within_FASTA.py`. However, 
# since may want the summary/accounting of letters/nts without invoking 
# making a file where anything is not A,T,G,C or N gets replaced, separated 
# that portion out to this script.)
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
# - verify compatible with Python 2.7
#
#
#
# TO RUN:
# Examples,
# Enter on the command line of your terminal, the line
#-----------------------------------
# python summarize_all_nts_even_ambiguous_present_in_FASTA.py seq.fa 
#-----------------------------------
#
# Issue `summarize_all_nts_even_ambiguous_present_in_FASTA.py -h` for details.
# 
#
#
# To use this after importing/pasting or loading into a cell in a Jupyter 
# notebook, call the main function similar to below:
# summarize_all_nts_even_ambiguous_present_in_FASTA("<name_of_sequence_file>.fa")
#
# 
#
'''
CURRENT ACTUAL CODE FOR RUNNING/TESTING IN A NOTEBOOK WHEN IMPORTED/LOADED OR 
PASTED IN ANOTHER CELL:
df = summarize_all_nts_even_ambiguous_present_in_FASTA("seq.fa")
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

#NONE YET; PLACEHOLDER FOR SECTION, FOR NOW.

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
import pandas as pd
import collections
try:
    from rich_dataframe import prettify
except ImportError:
    sys.stderr.write("Run `pip install rich-dataframe` on your command line "
        "or if\nusing in a Jupyter notebook, run `%pip install rich-dataframe`"
        ", without the tick marks.")




###---------------------------HELPER FUNCTIONS-------------------------------###
'''
from typing import Optional
def df_to_table(
        pandas_dataframe: pd.DataFrame,
        rich_table: Table,
        show_index: bool = True,
        index_name: Optional[str] = None,
    ) -> Table:
    """Convert a pandas.DataFrame obj into a rich.Table obj.
    Args:
        pandas_dataframe (DataFrame): A Pandas DataFrame to be converted to a rich Table.
        rich_table (Table): A rich Table that should be populated by the DataFrame values.
        show_index (bool): Add a column with a row count to the table. Defaults to True.
        index_name (str, optional): The column name to give to the index column. Defaults to None, showing no value.
    Returns:
        Table: The rich Table instance passed, populated with the DataFrame values.

    from: https://gist.github.com/neelabalan/33ab34cf65b43e305c3f12ec6db05938
    """

    if show_index:
        index_name = str(index_name) if index_name else ""
        rich_table.add_column(index_name)

    for column in pandas_dataframe.columns:
        rich_table.add_column(str(column))

    for index, value_list in enumerate(pandas_dataframe.values.tolist()):
        row = [str(index)] if show_index else []
        row += [str(x) for x in value_list]
        rich_table.add_row(*row)

    return rich_table
'''


def generate_output_file_name(file_name,suffix_for_saving):
    '''
    Takes a file name as an argument and returns string for the name of the
    output file. The generated name is based on the original file name.


    Specific example
    =================
    Calling function with
        ("sequence.fa", "_subbed)
    returns
        "sequence_subbed.fa"
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


def percent_calc(items):
    '''
    takes a list of two items and calculates percentage of first item
    within total (second item)
    '''
    return items[0]/items[1]

def make_count_of_nucleotides(seq):
    '''
    Take a sequence and return a dictionary of the counts of each letter in the
    sequence.
    '''
    return collections.Counter(seq)

def make_dataframe_accounting_of_nucleotides(dict_of_seq_letter_counts):
    '''
    Take a dictonary of counts of letters in the sequences and returns a 
    dataframe putting the data for each sequence as a row with nice summary
    features with respect totals and %N.
    Quantification code adapted from my 'Assessing_ambiguous_nts..' series of 
    notebooks.
    '''
    nt_count_df = pd.DataFrame.from_dict(
        dict_of_seq_letter_counts, orient='index').fillna(0)
    # Because high quality sequences won't necessarily have any N's, add a column
    # including that so below `%N` can be caculated without triggering  
    # `KeyError: "['N'] not in index"`.
    if 'N' not in nt_count_df.columns:
        nt_count_df['N'] = 0
    nt_count_df["Total_nts"] = nt_count_df.sum(1).astype(dtype='int64') 
    nt_count_df['% N'] = nt_count_df[['N','Total_nts']].apply(percent_calc, axis=1)
    nt_count_df.loc['TOTAL']= nt_count_df.sum().astype(dtype='int64') 
    return nt_count_df

def display_df_in_terminal(df):
    '''
    Display provided dataframe in command line/terminal nicely.

    Takes: a Pandas dataframe
    Returns: Nothing, just handles dsiplaying

    Moved to a helper function so can edit this one function and still get same 
    output other places that use core of this script such as in
    `replace_unusual_nts_within_FASTA.py`.
    '''
    prettify(df, row_limit=100,col_limit=100, 
       delay_time = 1, clear_console=False)

    #from rich import box
    #from rich.console import Console
    #from rich.table import Table

    #console = Console()
    #table = Table(show_header=True, header_style="bold magenta")

    # Modify the table instance to have the data from the DataFrame
    #table = df_to_table(df, table)

    # Update the style of the table
    #table.row_styles = ["none", "dim"]
    #table.box = box.SIMPLE_HEAD

    #console.print(table)

###--------------------------END OF HELPER FUNCTIONS-------------------------###
###--------------------------END OF HELPER FUNCTIONS-------------------------###



#*******************************************************************************
###------------------------'main' function of script--------------------------##

def summarize_all_nts_even_ambiguous_present_in_FASTA(
    sequence_file, return_df = True, display_text_of_dataframe_df = False):
    '''
    Main function of script.
    Takes a sequence string, a pattern, and a record id and 
    deletes all the seqeunce following the match to the sequence pattern string
    in the specified record. 
    The FASTA-formatted sequence file is assumed by default to be a 
    multi-FASTA, i.e., multiple sequences in the provided file, although it 
    definitely doesn't have to be. In case it is only a single sequence, the 
    record id becomes moot and users can provide anything for this parameter. 
    It also will count all the nucleotides in the data and produce a dataframe
    that can be printed in the terminal using rich-dataframe when called from
    the command line.

    A modified sequence file will be made.
    '''
    # ITERATE ON SEQUENCES IN SEQUENCE FILE:
    #--------------------------------------------------------------------------#
    # make a records list because will most likely need to modify a record. 
    # Plus, knowing how many records there are will affect how records handled
    # because don't want record_id to matter in such a case & 
    # be sure to give feedback of that sort.
    records = []  #intialize a list of records for the FASTA file
    unusual_nt_encountered = False # set this to False so have flag for knowing
    # if unusual encountered. Assume an unsual nt will get encountered 
    # (otherwise why would this script get selected for use?), & so in the 
    # interest being ready, run the substitution on all sequences as they are 
    # encounteted, but only outputfile with new ones if  
    # `unusual_nt_encountered` is True
    for record in SeqIO.parse(sequence_file, "fasta"):
        records.append(record)
    single_record = False
    if len(records) == 1:
        record_id = records[0].id
        single_record = True
        # feedback
        sys.stderr.write("Single sequence with id of '{}' provided in the "
            "sequence file.\nIt will be and alyzed & the unusual characters "
            "replaced with {}.\n\n".format(record_id, character_for_subbing))
    else:
         sys.stderr.write("{} sequences provided in the "
            "sequence file.\n\n".format(len(records)))

    

    # Go through each record and collect the counts on all the letters in the 
    # sequences. Also, if any sequences have anything other than A,T,G,C or N,
    # make new versions of the records where the characters that aren't 
    # A,T,G,C or N are substituted with the specifed character.
    dict_of_seq_letter_counts = {}
    for indx,record in enumerate(records):
        # allow for sequences with same record.id (I don't know if biopython 
        # allows that) If it happens appending the sequence number in all 
        # the entries will make it unique, and so if it gets through
        # biopython, at least it will be handled.
        if record.id in dict_of_seq_letter_counts:
            record_id = record.id+f"_{indx}"
        else:
            record_id = record.id
        dict_of_seq_letter_counts[record.id] = make_count_of_nucleotides(
            str(record.seq).upper())

    # Now that everything has been analyzed, make summary of letter counts as a 
    # dataframe
    df = make_dataframe_accounting_of_nucleotides(dict_of_seq_letter_counts)


    # display dataframe summary in terminal if using command line
    if display_text_of_dataframe_df:
        display_df_in_terminal(df)
    elif return_df:
        return df
    

    
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
    kwargs['return_df'] = False #probably don't want dataframe returned if 
    # calling script from command line
    kwargs['display_text_of_dataframe_df'] = True # if calling script from 
    # command line, a try to use rich-dataframe to print contents of count 
    # dataframe to terminal  
    summarize_all_nts_even_ambiguous_present_in_FASTA(
        sequence_file, **kwargs)
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
        'summarize_all_nts_even_ambiguous_present_in_FASTA.py',
        description="summarize_all_nts_even_ambiguous_present_in_FASTA.py \
        takes a sequence file (FASTA-format) and summarizes the nts in the \
        sequence(s). Assumes multi-FASTA, but single sequence entry is fine, \
        too. When running on the command line, it will print out a summary \
        table of counts of nucleotides and other character in each sequence \
        and totals. When calling the main function it will, by default, return \
        a dataframe with this information. Only valid for DNA sequences; \
        script has no step checking for data type, and so you are responsible \
        for verifying appropriate input. \
        **** Script by Wayne Decatur   \
        (fomightez @ github) ***")

    parser.add_argument("sequence_file", help="Name of sequence file to \
        use as input. Must be FASTA format. Can be a \
        multi-FASTA file, i.e., multiple sequences in FASTA format in one \
        file.", metavar="SEQUENCE_FILE")
    






    #I would also like trigger help to display if no arguments provided because 
    # need at least one for url
    if len(sys.argv)==1:    #from http://stackoverflow.com/questions/4042452/display-help-message-with-python-argparse-when-script-is-called-without-any-argu
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    sequence_file = args.sequence_file




    main()

#*******************************************************************************
###-***********************END MAIN PORTION OF SCRIPT***********************-###
#*******************************************************************************
