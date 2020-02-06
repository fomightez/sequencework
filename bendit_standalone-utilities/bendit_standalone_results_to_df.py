#!/usr/bin/env python
# bendit_standalone_results_to_df.py
__author__ = "Wayne Decatur" #fomightez on GitHub
__license__ = "MIT"
__version__ = "0.1.0"


# bendit_standalone_results_to_df.py by Wayne Decatur
# ver 0.1
#
#*******************************************************************************
# Verified compatible with both Python 2.7 and Python 3.6; written initially in 
# Python 3. 
#
#
# PURPOSE: Takes raw output from command line-based standalone version of bendIt 
# & brings it into Python as a Pandas dataframe and saves a file of that 
# dataframe for use elsewhere. 
# Optionally, it can also return that dataframe for use inside a Jupyter 
# notebook.
#
# This script is meant to be a utility script for working with the standalone 
# version of bendIt and Python, see a demonstration of use in
# https://git.io/Jvsad
# 
# If you are working with bend.it server results, see
# https://github.com/fomightez/sequencework/tree/master/bendit_server-utilities
# 
#
#
#
# Written to run from command line or imported into a Jupyter notebook. 
#
#
#
# 
#
#
#
#
# Developed by adapting backbone of `blast_to_df.py`.
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
# - verify works in Python 2 (low priority)
# - check which Python 3 version verified it works in and update comments above
# - if encounter any special way bendIt has to be run note above int he header.
# - update examples below as settle on syntax
# - when settle if making one or two dataframes, edit text in description of 
# function `def bendit_standalone_results_to_df()` to reflect if one or both?
# - update comment just below 'Read the bendIt results & process to get ready for going into dataframe:' if determine why positions listed twice
#
#
#
#
# TO RUN:
# Examples,
# Enter on the command line of your terminal, the line
#-----------------------------------
# python bendit_standalone_results_to_df.py RESULTS_FILE
#-----------------------------------
# Issue `bendit_standalone_results_to_df.py -h` for details.
# 
# More examples from running from the command line are at the links below: 
# https://git.io/vh8M9
#
#
#
# To use this after importing/pasting or loading into a cell in a Jupyter 
# notebook, specify at least the results file (or results as a string) in the 
# call to the main function similar to below:
# df = bendit_standalone_results_to_df("test.out")
# df
#
#
# A more in-depth series of examples of using this script within a notebook 
# without need to save file intermediates is found at:
# https://git.io/vh8M9
#
# 
#
'''
CURRENT ACTUAL CODE FOR RUNNING/TESTING IN A NOTEBOOK WHEN IMPORTED/LOADED OR 
PASTED IN ANOTHER CELL:
df = bendit_standalone_results_to_df("test.out")
df
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

## Settings and options for output
df_save_as_name = 'bendit_pickled_df.pkl' # name for saving pickled dataframe
text_save_as_name = 'bendit_results.tsv'


#
#*******************************************************************************
#**********************END USER ADJUSTABLE VARIABLES****************************













#*******************************************************************************
#*******************************************************************************
###DO NOT EDIT BELOW HERE - ENTER VALUES ABOVE###

import sys
import os



###---------------------------HELPER FUNCTIONS---------------------------------###

def write_string_to_file(s, fn):
    '''
    Takes a string and a name for a file and writes the string to the file.
    '''
    with open(fn, 'w') as output_file:
        output_file.write(s)

###--------------------------END OF HELPER FUNCTIONS---------------------------###
###--------------------------END OF HELPER FUNCTIONS---------------------------###

#*******************************************************************************
###------------------------'main' function of script---------------------------##

def bendit_standalone_results_to_df(results, also_reported, return_df = True, 
    pickle_df=True, sequence_string = "NONE_WAS_PROVIDED", 
    sequence_file = "NONE_WAS_PROVIDED"):
    '''
    Main function of script. 
    raw bendit standalone version results to a Pandas dataframe.

    It will take a text file of raw results from standalone verison of bendIt & 
    make a dataframe that will be more useful with Python or Juputer contexts.

    Optionally also returns a dataframe of the results data. Meant for use in 
    a Jupyter notebook.

    There are two OPTIONAL ways to provide a sequence as discussed in the next 
    two paragraphs. This will result in output more similar to the bendit server 
    produces that includes a column for the sequence. The sequence being 
    provided is optional; without the sequence provided, the 'Seqence' column
    will simply be absent in the text tables/dataframes produced.

    The option to provide the sequence as a string is to handle where using the 
    main function of the script imported into IPython or jupyter and therefore 
    can dispense with an intermediate file (FASTA format) to passs in the 
    sequence. Note this setting takes priority over any setting `sequence_file`
    may have. And so it is meant to really be using one or the other.

    The option to provide a sequence as a FASTA file is meant for when using 
    the script as script, i.e.m when calling with `%run ...` or `python ...`. 

    By default, the dataframe is pickled but this can be opted out of as 
    well.

    The astute user will have noted the raw output from the standalone version 
    of the bendIt software includes two sets of the same results for each 
    position. At present, I don't know why the results are repeated in the 
    output. The documentation for the standalone version of bendIt is very 
    sparse, and I haven't combed through the code to see if any reason is noted
    there. Whatever, the reason, this script verifies that both versions are 
    indeed the same and uses the first one.
    '''
    import pandas as pd
    if sequence_file == "NONE_WAS_PROVIDED" and (
        sequence_string == "NONE_WAS_PROVIDED"):
        col_names = ['Position', 'Predicted_curvature', also_reported]
    else:
        col_names = ['Position','Sequence','Predicted_curvature', also_reported]
    # If provided, bring in the sequence from the sequence file
    if sequence_string == "NONE_WAS_PROVIDED" and (sequence_file != "NONE_WAS_PROVIDED"):
        from pyfaidx import Fasta
        sequence_entries = Fasta(sequence_file, read_long_names=True)
        sequence_string = sequence_entries[0]
    # Read the bendIt results & process to get ready for going into dataframe:
    #---------------------------------------------------------------------------
    # read line by line so can skip past gnuplot instructions header, add in 
    # sequence data if needed, and split out what seems to be results twice(<====still deciphering why every position listed twice)
    # note it looks like gnuplot instructions header would end after tenth line
    # but I am not 100% sure that is 100% consistent and so best to skip header 
    # is to first skip anything that begins with `set` or `plot` and then 
    # start parsing next lines once hit `plot`
    # Because I am seeing each position listed twice in output from standalone 
    # bendIt, parsing each set to separate temporary file for comparison:
    # 1. `results_parsing_temp1.tsv`
    # 2. `results_parsing_temp2.tsv`
    gnuplot_instruction_header_passed = False
    gnuplot_header_line_starts = ("set","plot","unset","#", "replot", "splot") # 
    # I added more beyond `set` & `plot` that I saw in examples from running 
    # bendit based on https://people.duke.edu/~hpgavin/gnuplot.html ; it has to
    # be tuples and not a list for `.startswith()` method
    data_end_marker = "e"
    temp_files_for_passing_data = ["results_parsing_temp1.tsv",
                                    "results_parsing_temp2.tsv"]
    file_to_save = temp_files_for_passing_data[0]
    growing_string_to_save = "\t".join(col_names) +"\n"# add column names for easier
    # reading by pandas later
    with open(results, 'r') as input:
        # prepare to give feeback later, or MORE LIKELY, allow skipping to 
        # certain start point if ever determine header is ALWAYS ten lines. 
        # Right now being conservative & not just skipping by line number but 
        # screening contents
        lines_processed = 0
        for line in input:
            lines_processed += 1
            if line.startswith("plot"):
                gnuplot_instruction_header_passed = True
                print(f"Header ended on line {lines_processed}.") #temp. to check if ever deviates
                continue
            if line.startswith(gnuplot_header_line_starts) and not gnuplot_instruction_header_passed:
                continue
            if line.startswith(data_end_marker) and (
                line.strip() == data_end_marker):
                write_string_to_file(growing_string_to_save, file_to_save)
                # reset for next potential (block). If that had been last block,
                # no harm resetting these.
                file_to_save = temp_files_for_passing_data[1]
                growing_string_to_save = "\t".join(col_names)+"\n"
            else:
                if 'Sequence' in col_names:
                    # need to add in sequence alongside results if sequence 
                    # provided to get something closer to bendit server output
                    line_components = line.split()
                    line_components = [x.strip() for x in line_components]
                    corresponding_nt = sequence_string[int(line_components[0])]#
                    # <---output seems zero-indexed in output from bendit 
                    # server, & so output from bendit alone likely same and so 
                    # no need to adjust because Python also zero-indexed.
                    new_line = ("{}\t{}\t{}\t{}\n".format(
                        (line_components[0]),corresponding_nt,
                        line_components[1],line_components[2]))
                    growing_string_to_save += new_line
                else:
                    growing_string_to_save += line


    # Make collected results into dataframe:
    #---------------------------------------------------------------------------
    df = {}
    for i,f in enumerate(temp_files_for_passing_data):
        df[i] = pd.read_csv(f, sep='\t', header=0)

    # Test the two sets of data the standalone version of bendIt reports in the
    # output file are the same
    assert df[0].equals(df[1]), ("The two dataframes are different. NOT SEEN "
    "BEFORE. What is different?")
    # Seems to be the case that they are always the same so making it an assert
    # because it would be BIG if not! Since now assuming are AND ONLY PASSING
    # FIRST BACK.
    '''
    if not df[0].equals(df[1]):
        sys.stderr.write("***WARNING: The two sets of data in the results file "
            "are different. Maybe an issue?")
    '''
    i = 0 # use first
    # FROM NOW ON, THIS SCRIPT WILL USE THE FIRST SET OF RESULTS bendIt GIVES.



    # feedback
    sys.stderr.write("Provided results read and converted to a dataframe...")


    # Reporting and Saving
    #---------------------------------------------------------------------------
    #print(df)#originally for debugging during development,added..
    # Document the full set of data collected in the terminal or 
    # Jupyter notebook display in some manner. 
    # Using `df.to_string()` because more universal than `print(df)` 
    # or Jupyter's `display(df)`.
    #sys.stderr.write("\nFor documenting purposes, the following lists the "
    #    "parsed data:\n")
    #with pd.option_context('display.max_rows', None, 'display.max_columns', None):
    #    display(df)
    #sys.stderr.write(df.to_string())

    # Handle pickling the dataframe
    ''' # COMMENTED OUT THE PICKLING AND SAVING AS TEXT OF BOTH DATAFRAMES SINCE CHECKING A LOT SUGGEST THEY ARE ALWAYS THE SAME
    df_save_as_names = [df_save_as_name,"second_"+df_save_as_name]
    text_save_as_names = [text_save_as_name,"second_"+text_save_as_name]
    for i,f in enumerate(temp_files_for_passing_data):
        if pickle_df == False:
            sys.stderr.write("\n\nA dataframe of the data "
            "was not stored for use\nelsewhere "
            "because `no_pickling` was specified.\n")
        else:
            df[i].to_pickle(df_save_as_names[i])
            # Let user know
            sys.stderr.write("\n\nA dataframe of the data "
            "has been saved as a file\nin a manner where other "
            "Python programs can access it (pickled form).\n"
            "RESULTING DATAFRAME is stored as ==> '{}'".format(
                df_save_as_names[i] ))
            df[i].to_csv(text_save_as_names[i] , sep='\t',index = False) 
            # Let user know
            sys.stderr.write("\n\nA text table of the data "
            "has been saved as a file\nin tab-delimited form.\n"
            "RESULTING TAB-SEPARATED TEXT FILE is stored as ==> '{}'\n".format(
                df_save_as_names[i] ))
    '''
    if pickle_df == False:
        sys.stderr.write("\n\nA dataframe of the data "
        "was not stored for use\nelsewhere "
        "because `no_pickling` was specified.\n")
    else:
        df[i].to_pickle(df_save_as_name)
        # Let user know
        sys.stderr.write("\n\nA dataframe of the data "
        "has been saved as a file\nin a manner where other "
        "Python programs can access it (pickled form).\n"
        "RESULTING DATAFRAME is stored as ==> '{}'".format(
            df_save_as_name ))
        df[i].to_csv(text_save_as_name , sep='\t',index = False) 
        # Let user know
        sys.stderr.write("\n\nA text table of the data "
        "has been saved as a file\nin tab-delimited form.\n"
        "RESULTING TAB-SEPARATED TEXT FILE is stored as ==> '{}'\n".format(
            df_save_as_name ))

    
    # Return dataframe (optional)
    #---------------------------------------------------------------------------
    if return_df:
        sys.stderr.write("\n\nReturning a dataframe with the information "
                "as well.")
        return df[0]

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
    if df_save_as_name == 'no_pickling':
        kwargs['pickle_df'] = False
    kwargs['return_df'] = False #probably don't want dataframe returned if 
    # calling script from command line
    if sequence_file:
        kwargs['sequence_file'] = sequence_file
    else:
        kwargs['sequence_file'] = "NONE_WAS_PROVIDED"
    bendit_standalone_results_to_df(results,**kwargs)
    # using https://www.saltycrane.com/blog/2008/01/how-to-use-args-and-kwargs-in-python/#calling-a-function
    # to build keyword arguments to pass to the function above
    # (see https://stackoverflow.com/a/28986876/8508004 and
    # https://stackoverflow.com/a/1496355/8508004 
    # (maybe https://stackoverflow.com/a/7437238/8508004 might help too) for 
    # related help). Makes it easy to add more later.





#if __name__ == "__main__" and '__file__' in globals():
if __name__ == "__main__":
    """ This is executed when run from the command line """
    # Code with just `if __name__ == "__main__":` alone will be run if pasted
    # into a notebook. The addition of ` and '__file__' in globals()` is based
    # on https://stackoverflow.com/a/22923872/8508004
    # See also https://stackoverflow.com/a/22424821/8508004 for an option to 
    # provide arguments when prototyping a full script in the notebook.
    ###-----------------for parsing command line arguments-------------------###
    import argparse
    parser = argparse.ArgumentParser(prog='bendit_standalone_results_to_df.py',
        description="bendit_standalone_results_to_df.py \
        Takes raw output from command line-based standalone version of bendIt \
        & brings it into Python as a dataframe and saves a file of that \
        Pandas dataframe for use elsewhere. \
        **** Script by Wayne Decatur   \
        (fomightez @ github) ***")

    parser.add_argument("results_file", help="Name of file of bendIt results \
        file to parse.\
        ", metavar="RESULTS_FILE")
    # choices below based on https://docs.python.org/3/library/argparse.html and 
    # https://stackoverflow.com/a/20335589/8508004
    parser.add_argument("also_reported", choices=['bendability', 'complexity', 
        'gc_content'],help='also_reported, {%(choices)s}. This option is to \
        correspond to what the `-g` setting was set as in the call to bendIt \
        or what I call in my associated code as `report_with_curvature` \
        because it is what you choose to analyze and report along with \
        curvature when running bendIt.')
    parser.add_argument('-dfo', '--df_output', action='store', type=str, 
    default= df_save_as_name, help="OPTIONAL: Set file name for saving pickled \
    dataframe. If none provided, '{}' will be used. To force no dataframe to \
    be saved, enter `-dfo no_pickling` without quotes as output file \
    (ATYPICAL).".format(df_save_as_name))
    parser.add_argument("sequence_file", nargs='?', help="**OPTIONAL**Name of file \
        that has sequence in FASTA format. If not included, there simply won't\
        be a `Sequence` column in the resulting dataframe.", 
        metavar="SEQUENCE_FILE")
    # See
    # https://stackoverflow.com/questions/4480075/argparse-optional-positional-arguments 
    # and 
    # https://docs.python.org/2/library/argparse.html#nargs for use of `nargs='?'` 
    # to make output file name optional. Note that the square brackets
    # shown in the usage out signify optional according to 
    # https://stackoverflow.com/questions/4480075/argparse-optional-positional-arguments#comment40460395_4480202
    # , but because placed under positional I added clarifying text to help 
    # description.
    # IF MODIFYING THIS SCRIPT FOR USE ELSEWHERE AND DON'T NEED/WANT THE OUTPUT 
    # FILE TO BE OPTIONAL, remove `nargs` (& default?) BUT KEEP WHERE NOT
    # USING `argparse.FileType` AND USING `with open` AS CONISDERED MORE PYTHONIC.



    #I would also like trigger help to display if no arguments provided because 
    # need at least one for url
    if len(sys.argv)==1:    #from http://stackoverflow.com/questions/4042452/display-help-message-with-python-argparse-when-script-is-called-without-any-argu
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    results = args.results_file
    also_reported = args.also_reported
    df_save_as_name = args.df_output
    sequence_file = args.sequence_file


    main()

#*******************************************************************************
###-***********************END MAIN PORTION OF SCRIPT***********************-###
#*******************************************************************************
