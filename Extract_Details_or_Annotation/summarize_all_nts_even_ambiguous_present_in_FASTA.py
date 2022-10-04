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
# Case-insensitive by default; however, both main routes to use the script /
# main function allow for setting it to 'case-sensitive' so that the numbers of
# lowercase & uppercase letters are reported separately.
#
# Dependent on the package 'rich-dataframe'
#
# (incorporates quantifying code adapted from my 
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
# - check if using case-sensitive setting necessary if have things like asterisk
# or other symbols in a sequence, meaning does `.upper()` cause any issue for
# when symbols are in sequence when run with default setting?
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
# Before built-in my customized version of rich-dataframe, this was how I was
# handling import. Search 'CUSTOMIZED RICH-DATAFRAME' in this document to find
# where added.
'''
try:
    from rich_dataframe import prettify
except ImportError:
    sys.stderr.write("Run `pip install rich-dataframe` on your command line "
        "or if\nusing in a Jupyter notebook, run `%pip install rich-dataframe`"
        ", without the tick marks.\n**EXITING.**.\n")
    sys.exit(1)
'''


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
    Take a sequence and return a dictionary of the counts of each character in 
    the sequence.
    '''
    return collections.Counter(seq)

def make_dataframe_accounting_of_nucleotides(dict_of_seq_letter_counts,
    case_sensitive=False):
    '''
    Take a dictonary of counts of letters in the sequences and returns a 
    dataframe putting the data for each sequence as a row with nice summary
    features with respect totals and %N.
    Quantification code adapted from my 'Assessing_ambiguous_nts..' series of 
    notebooks.
    By default it is case-insensitive; however, the lowercase & uppercase 
    letters will be tallied separately if `case_sensitive = True` in the call
    to this function.
    '''
    nt_count_df = pd.DataFrame.from_dict(
        dict_of_seq_letter_counts, orient='index').fillna(0)
    # Because high quality sequences won't necessarily have any N's, add a 
    # column including that so below `%N` can be caculated without triggering  
    # `KeyError: "['N'] not in index"`. Do same for lowercase `n` in case of
    # case-sensitivity setting.
    if 'N' not in nt_count_df.columns:
        nt_count_df['N'] = 0
    if case_sensitive and 'n' not in nt_count_df.columns:
        nt_count_df['n'] = 0
    nt_count_df["Total_nts"] = nt_count_df.sum(1).astype(dtype='int64') 
    nt_count_df['% N'] = nt_count_df[['N','Total_nts']].apply(
        percent_calc, axis=1)
    if case_sensitive:
        nt_count_df['% n'] = nt_count_df[['n','Total_nts']].apply(
            percent_calc, axis=1)
        extra_nt_count_df = nt_count_df[['N','n','Total_nts']].copy() # 
        # make a separate dataframe for combining counts of `N` & `n` without 
        # adding that column to main dataframe. The combined counts can then be
        # used to calculate total percent if consider both lower and upper case 
        # Ns as one. (Note that I added use of `.copy` to avoid 
        # `SettingWithCopyWarning` as per https://www.dataquest.io/blog/settingwithcopywarning/
        extra_nt_count_df['Nn'] = extra_nt_count_df[['N','n']].sum(
            1).astype(dtype='int64') 
        nt_count_df['% N&n'] = extra_nt_count_df[['Nn','Total_nts']].apply(
            percent_calc, axis=1)
    nt_count_df.loc['TOTAL']= nt_count_df.sum().astype(dtype='int64') 
    return nt_count_df



###--------------------------END OF HELPER FUNCTIONS-------------------------###
###--------------------------END OF HELPER FUNCTIONS-------------------------###



#*******************************************************************************
###------------------------'main' function of script--------------------------##

def summarize_all_nts_even_ambiguous_present_in_FASTA(
    sequence_file, case_sensitive=False, mono=False, return_df = True, 
    display_text_of_dataframe_df = False):
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
    the command line. The summary table displayed in the terminal will feature 
    colors as provided by rich-dataframe's features unless you set `mono` to 
    `True`.
    By default the counting of the letters will be case-insensitive; however,
    setting case_sensitive to `True` will cause the numbers of lowercase &
    uppercase letters to be tallied separately in the summary data table.

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

    

    # SUMMARIZE LETTERS PRESENT IN SEQUENCES:
    #--------------------------------------------------------------------------#
    # Go through each record and collect the counts on all the letters in the 
    # sequences. Also, if any sequences have anything other than A,T,G,C or N,
    # make new versions of the records where the characters that aren't 
    # A,T,G,C or N are substituted with the specifed character.
    dict_of_seq_letter_counts = {}
    for indx,record in enumerate(records):
        # allow for sequences with same record.id (I don't know if biopython 
        # allows that) If it happens, then appending the position number of that
        # record among all the entries will make it unique, and so if it gets 
        # through biopython, at least it will be handled.
        if record.id in dict_of_seq_letter_counts:
            record_id = record.id+f"_{indx}"
        else:
            record_id = record.id
        if case_sensitive: 
            dict_of_seq_letter_counts[record.id] = make_count_of_nucleotides(
            str(record.seq))
        else:
            dict_of_seq_letter_counts[record.id] = make_count_of_nucleotides(
            str(record.seq).upper()) 

    # Now that everything has been analyzed, make summary of letter counts as a 
    # dataframe
    df = make_dataframe_accounting_of_nucleotides(dict_of_seq_letter_counts, 
        case_sensitive)


    # HANDLE DISPLAY OR RETURN OF SUMMARY INFO:
    #--------------------------------------------------------------------------#
    # display dataframe summary in terminal if using command line
    if display_text_of_dataframe_df:
        prettify(df, row_limit=100,col_limit=100, 
           delay_time = 1, clear_console=False, mono=mono)
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
    kwargs['case_sensitive'] = case_sensitive
    kwargs['mono'] = mono
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
        for verifying appropriate input. By default, the summary is \
        case-insensitive; however, a flag can be added at the time of calling \
        the script so that it will tally the lowercase & uppercase letters \
        separately in the summary data table. You can also choose to have \
        the summary table display in the terminal in a simpler, color-less \
        form. \
        **** Script by Wayne Decatur   \
        (fomightez @ github) ***")

    parser.add_argument("sequence_file", help="Name of sequence file to \
        use as input. Must be FASTA format. Can be a \
        multi-FASTA file, i.e., multiple sequences in FASTA format in one \
        file.", metavar="SEQUENCE_FILE")

    parser.add_argument('-cs', '--case_sensitive', help="Add this flag when \
        calling the script if you want the numbers of lowercase & uppercase \
        letters reported separately in the summary. \
        ", action="store_true")
    parser.add_argument('-mc', '--mono', help="Add this flag when \
        calling the script if you want the display of the summary in the \
        terminal to be black-and-white. \
        ", action="store_true")
    






    #I would also like trigger help to display if no arguments provided because 
    # need at least one for url
    if len(sys.argv)==1:    #from http://stackoverflow.com/questions/4042452/display-help-message-with-python-argparse-when-script-is-called-without-any-argu
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    sequence_file = args.sequence_file
    case_sensitive = args.case_sensitive
    mono = args.mono


    ###--------------------CUSTOMIZED RICH-DATAFRAME-------------------------###
    # Based, almost directly, on https://github.com/khuyentran1401/rich-dataframe 
    # Building it in because it is small and hasn't been updated much in a while &
    # most importantly because I want to customize how it handles the captioning.
    # I put this in my fork of https://github.com/khuyentran1401/rich-dataframe now .
    # SPECIFIC CUSTOMIZATIONS:
    #-------------------------
    # I don't want the caption describing how many rows or cols shown unless number 
    # of rows or columns is larger than how many are shown.
    # Plus, adding note to install rich.
    # Plus, removing the animation to the 'beat' because it causes weird spacing
    # in Jupyter and I had already turned the speed up so high because I wasn't
    # interested in the animated aspect that the default rich-dataframe makes.
    # To do this and get output to show up even when using `%run` in Jupyter,
    # I had to define the color and style of the columns when first made instead
    # of updating after using the beat, `_add_random_color()`, & `_add_style()`,
    # AND DELETE `with Live()`.
    # Plus, removed `_change_width()` section since I'm not seeing a difference
    # without it when not running animation.
    # Plus, I added a setting so you could just dislay in terminal or equivalent 
    # in monochrome if you didn't like the color features Rich-dataframe adds
    # but like the table styling it allows. More like a plainer style but the 
    # table still is easy to read in the terminal, like 
    # https://stackoverflow.com/a/72747970/8508004 , but better because header 
    # handled by rich-dataframe.
    #-------------------------
    # Development/Trouble-shooting cycle that worked best for developing my
    # custom implementation:
    # Starting a MyBinder session from 
    # https://github.com/binder-examples/requirements and installing only `rich`
    # and then runningan edited version of `example.py` from 
    # https://github.com/fomightez/rich-dataframe allowed me to see I shouldn't 
    # need `rich-dataframe` installed to get run `%run example.py` to show the 
    # output in a JupyterLab cell. EDITS: `example.py` had the entire 
    # `rich_dataframe.py` contents placed in it. I also changed what data it uses
    # because it was annoying to get the large data it used and I wanted a 
    # dataframe I was familiar with. I used the iris dataset that's built into 
    # seaborn, see 
    # https://github.com/mwaskom/seaborn-data . `iris = pd.read_csv('https://raw.githubusercontent.com/mwaskom/seaborn-data/master/iris.csv')` or 
    # !curl -OL https://raw.githubusercontent.com/mwaskom/seaborn-data/master/iris.csv 
    # can be used to use it with pandas or get it, respecitivey.  
    # I didn't want to install seaborn to keep things simple in the session. Note that 
    # `rich-dataframe` didn't seem to disply `iris` dataframe with the settings 
    # in `example.py`  and so I further changed `example.py` to remove
    # the reference to `first_rows` and `first_cols=False`.
    # Then I could edit the top section from `rich_dataframe.py` to see what 
    # could be removed to allow display of a static dataframe.
    #-------------------------
    # Adding this in this section of my script so only need rich installed if
    # using command line. Use of the main function directly wouldn't necessarily
    # need the code handling printng the dataframe in the terminal and so best
    # rich not required in such cases so not asking user to install something
    # that isn't used.
    try:
        from rich import print
    except ImportError:
        sys.stderr.write("Run `pip install rich` on your command line "
            "or if\nusing in a Jupyter notebook, run `%pip install rich`"
            ", without the tick marks.\n**EXITING.**.\n")
        sys.exit(1)
    from rich.box import MINIMAL, SIMPLE, SIMPLE_HEAD, SQUARE
    from rich.columns import Columns
    from rich.console import Console
    from rich.live import Live
    from rich.measure import Measurement
    from rich.table import Table
    console = Console()
    COLORS = ["cyan", "magenta", "red", "green", "blue", "purple"]
    class DataFramePrettify:
        """Create animated and pretty Pandas DataFrame

        Parameters
        ----------
        df : pd.DataFrame
            The data you want to prettify
        row_limit : int, optional
            Number of rows to show, by default 20
        col_limit : int, optional
            Number of columns to show, by default 10
        first_rows : bool, optional
            Whether to show first n rows or last n rows, by default True. If this is set to False, show last n rows.
        first_cols : bool, optional
            Whether to show first n columns or last n columns, by default True. If this is set to False, show last n rows.
        delay_time : int, optional
            How fast is the animation, by default 5. Increase this to have slower animation.
        clear_console: bool, optional
             Clear the console before printing the table, by default True. If this is set to False the previous console input/output is maintained
        """
        def __init__(
            self,
            df: pd.DataFrame,
            row_limit: int = 20,
            col_limit: int = 10,
            first_rows: bool = True,
            first_cols: bool = True,
            delay_time: int = 5,
            clear_console: bool = True,
            mono: bool = False,
        ) -> None:
            self.df = df.reset_index().rename(columns={"index": ""})
            self.table = Table(show_footer=False)
            self.table_centered = Columns(
                (self.table,), align="center", expand=True
            )
            self.num_colors = len(COLORS)
            self.delay_time = delay_time
            self.row_limit = row_limit
            self.first_rows = first_rows
            self.col_limit = col_limit
            self.first_cols = first_cols
            self.clear_console = clear_console
            self.mono = mono
            if first_cols:
                self.columns = self.df.columns[:col_limit]
            else:
                self.columns = list(self.df.columns[-col_limit:])
                self.columns.insert(0, "index")

            if first_rows:
                self.rows = self.df.values[:row_limit]
            else:
                self.rows = self.df.values[-row_limit:]
            
            if self.clear_console:
                console.clear()

        def _add_columns(self):
            for i,col in enumerate(self.columns):
                if self.mono:
                    self.table.add_column(str(col))
                else:
                    color4col = COLORS[i % self.num_colors]# based on https://github.com/khuyentran1401/rich-dataframe/blob/fff92c5fb735babcec580b88ef94b9325b5b8558/rich_dataframe/rich_dataframe.py#L110
                    self.table.add_column(str(col),style="bold "+color4col, header_style="bold "+color4col,) # based on https://rich.readthedocs.io/en/stable/tables.html#tables 
                    # and https://github.com/khuyentran1401/rich-dataframe/blob/fff92c5fb735babcec580b88ef94b9325b5b8558/rich_dataframe/rich_dataframe.py#L110
        def _add_rows(self):
            for row in self.rows:
                if self.first_cols:
                    row = row[: self.col_limit]
                else:
                    row = row[-self.col_limit :]
                row = [str(item) for item in row]
                self.table.add_row(*list(row))
        def _move_text_to_right(self):
            for i in range(len(self.table.columns)):
                self.table.columns[i].justify = "right"
        def _adjust_box(self):
            for box in [SIMPLE_HEAD, SIMPLE, MINIMAL, SQUARE]:
                self.table.box = box
        def _dim_row(self):
            self.table.row_styles = ["none", "dim"]
        def _adjust_border_color(self):
            if not self.mono:
                self.table.border_style = "orange1" # normally `bright_yellow` in 
            # rich-dataframe
        def _add_caption(self):
            if self.first_rows:
                row_text = "first"
            else:
                row_text = "last"
            if self.first_cols:
                col_text = "first"
            else:
                col_text = "last"

            if (len(self.df) > self.row_limit) and (len(self.df.columns) > self.col_limit):
                self.table.caption = f"Only the [bold magenta not dim] {row_text} {self.row_limit} rows[/bold magenta not dim] and the [bold green not dim]{col_text} {self.col_limit} columns[/bold green not dim] are shown here."
            elif len(self.df) > self.row_limit:
                self.table.caption = f"Only the [bold magenta not dim] {row_text} {self.row_limit} rows[/bold magenta not dim] are shown here."
            elif len(self.df.columns) > self.col_limit:
                self.table.caption = f"Only the [bold green not dim]{col_text} {self.col_limit} columns[/bold green not dim] are shown here."

        def prettify(self):
            self._add_columns()
            self._add_rows()
            self._move_text_to_right()
            self._adjust_border_color()
            self._add_caption()
            console.print(self.table) # based on https://rich.readthedocs.io/en/stable/tables.html#tables and https://stackoverflow.com/a/72747970/8508004
            return self.table
    def prettify(
        df: pd.DataFrame,
        row_limit: int = 20,
        col_limit: int = 10,
        first_rows: bool = True,
        first_cols: bool = True,
        delay_time: int = 5,
        clear_console: bool = True,
        mono: bool = False,
    ):
        """Create animated and pretty Pandas DataFrame

        Parameters
        ----------
        df : pd.DataFrame
            The data you want to prettify
        row_limit : int, optional
            Number of rows to show, by default 20
        col_limit : int, optional
            Number of columns to show, by default 10
        first_rows : bool, optional
            Whether to show first n rows or last n rows, by default True. If this is set to False, show last n rows.
        first_cols : bool, optional
            Whether to show first n columns or last n columns, by default True. If this is set to False, show last n rows.
        delay_time : int, optional
            How fast is the animation, by default 5. Increase this to have slower animation.
        clear_console: bool, optional
            Clear the console before priting the table, by default True. If this is set to false the previous console input/output is maintained.
        mono: bool, optional
            Print a plain display without color, by default False. If this is set to true, then the table generated will be 'black-and-white', a.k.a. monochrome.
        """
        if isinstance(df, pd.DataFrame) or isinstance(df, pd.DataFrame):
            DataFramePrettify(
                df, row_limit, col_limit, first_rows, first_cols, delay_time,
                clear_console, mono
            ).prettify()

        else:
            # In case users accidentally pass a non-datafame input, use rich's 
            # print instead
            print(df)
    ###---------------END OF CUSTOMIZED RICH-DATAFRAME-----------------------###







    main()

#*******************************************************************************
###-***********************END MAIN PORTION OF SCRIPT***********************-###
#*******************************************************************************
