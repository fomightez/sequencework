#!/usr/bin/env python
# bendit_server_results_to_df.py
__author__ = "Wayne Decatur" #fomightez on GitHub
__license__ = "MIT"
__version__ = "0.1.0"


# bendit_server_results_to_df.py by Wayne Decatur
# ver 0.1
#
#*******************************************************************************
# Verified compatible with both Python 2.7 and Python 3.6; written initially in 
# Python 3. 
#
#
# PURPOSE: Takes output from the online bend.it server & 
# brings it into Python as a Pandas dataframe and saves a file of that dataframe 
# for use elsewhere. 
# Optionally, it can also return that dataframe for use inside a Jupyter 
# notebook.
#
# This script is meant to be a utility script for working with the server based 
# version of bend.it , available at 
# http://pongor.itk.ppke.hu/dna/bend_it.html#/bendit_form, along with Python. 
# If you are working with the standalone version of bendIt, see 
# https://github.com/fomightez/sequencework/tree/master/bendit_standalone-utilities
# and
# https://github.com/fomightez/bendit-binder
# 
# It working necessitates ??????.
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
# function `def bendit_server_results_to_df()` to reflect if one or both?
# - update comment just below 'Read the bendIt results & process to get ready for going into dataframe:' if determine why positions listed twice
#
#
#
#
# TO RUN:
# Examples,
# Enter on the command line of your terminal, the line
#-----------------------------------
# python bendit_server_results_to_df.py RESULTS_FILE
#-----------------------------------
# Issue `bendit_server_results_to_df.py -h` for details.
# 
#
#
#
# To use this after importing/pasting or loading into a cell in a Jupyter 
# notebook, specify at least the results file (or results as a string) in the 
# call to the main function similar to below:
# df = bendit_server_results_to_df("test.out")
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
df = bendit_server_results_to_df("test.out")
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
df_save_as_name = 'bendit_server_pickled_df.pkl'#name for saving pickled dataframe
text_save_as_name = 'bendit_server_results.tsv'

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
def extract_column_info(line):
    '''
    Takes a line looking like following cases and returns the column number and 
    label in the line of text:

    Examples of input lines:
    -----------------------

    Column 1 = Position
    Column 2 = Sequence
    Column 3 = Predicted curvature
    Column 4 = Bendability

    '''
    column_no = int(line.split("Column")[1].split("=")[0].strip())
    label = line.split("=")[1].strip()

    return column_no, label

###--------------------------END OF HELPER FUNCTIONS---------------------------###
###--------------------------END OF HELPER FUNCTIONS---------------------------###

#*******************************************************************************
###------------------------'main' function of script---------------------------##

def bendit_server_results_to_df(results, return_df = True, pickle_df=True):
    '''
    Main function of script. 
    bendit server version results to Pandas dataframe(s).

    By results, I mean what you'd get if you highlighted text from results page
    from server and pasted it into a text document.

    The big things it does is:
    It reads the output header block and adds the columns correctly.
    Parses out the data block.
    

    It will take a text file of results from online bendit server and 
    make a dataframe that will be more useful with Python or Juputer contexts.

    Optionally also returns a dataframe(s) of the results data. Meant for use in 
    a Jupyter notebook.

    By default, the dataframe(s) are pickled but this can be opted out of as 
    well.
    '''
    import pandas as pd
    header_passed = False
    column_labels = []
    data_collected = ""
    with open(results, 'r') as input:
        # prepare to give feeback later by counting lines
        lines_processed = 0
        for line in input:
            lines_processed += 1
            if not header_passed:
                if line.startswith("Column"):
                    column_no, label = extract_column_info(line)
                    if " " in label:
                        label = label.replace(" ", "_")
                    column_labels.append(label)
                    assert(column_labels[column_no - 1]==label), ("The columns "
                        "appear to be listed out of order relative the listing "
                        "in the input file.")
                if line.startswith("------"):
                    header_passed = True
                    print(f"Header ended on line {lines_processed}.") #temp. to check if ever deviates
                    data_collected = "\t".join(column_labels)+"\n"
                    continue
            else:
                if line.startswith("------"):
                    break
                else:
                    data_collected += line


    # Make collected results into dataframe:
    #---------------------------------------------------------------------------
    try:
        from StringIO import StringIO
    except ImportError:
        from io import StringIO
    df = pd.read_csv(StringIO(data_collected), sep='\t', header=0)

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
    if pickle_df == False:
        sys.stderr.write("\n\nA dataframe of the data "
        "was not stored for use\nelsewhere "
        "because `no_pickling` was specified.\n")
    else:
        df.to_pickle(df_save_as_name )
        # Let user know
        sys.stderr.write("\n\nA dataframe of the data "
        "has been saved as a file\nin a manner where other "
        "Python programs can access it (pickled form).\n"
        "RESULTING DATAFRAME is stored as ==> '{}'".format(df_save_as_name ))
        df.to_csv(text_save_as_name , sep='\t',index = False) 
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
    if df_save_as_name == 'no_pickling':
        kwargs['pickle_df'] = False
    kwargs['return_df'] = False #probably don't want dataframe returned if 
    bendit_server_results_to_df(results,**kwargs)
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
    parser = argparse.ArgumentParser(prog='bendit_server_results_to_df.py',
        description="bendit_server_results_to_df.py \
        Takes output from the online bend.it server & \
        brings it into Python as a dataframe and saves a file of that \
        Pandas dataframe for use elsewhere. \
        **** Script by Wayne Decatur   \
        (fomightez @ github) ***")

    parser.add_argument("results_file", help="Name of file of text copied from \
        results page of job run on bend.it server.\
        ", metavar="RESULTS_FILE")
    parser.add_argument('-dfo', '--df_output', action='store', type=str, 
    default= df_save_as_name, help="OPTIONAL: Set file name for saving pickled \
    dataframe. If none provided, '{}' will be used. To force no dataframe to \
    be saved, enter `-dfo no_pickling` without quotes as output file \
    (ATYPICAL).".format(df_save_as_name))
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
    df_save_as_name = args.df_output


    main()

#*******************************************************************************
###-***********************END MAIN PORTION OF SCRIPT***********************-###
#*******************************************************************************