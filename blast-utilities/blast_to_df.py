#!/usr/bin/env python
# blast_to_df.py
__author__ = "Wayne Decatur" #fomightez on GitHub
__license__ = "MIT"
__version__ = "0.1.0"


# blast_to_df.py by Wayne Decatur
# ver 0.1
#
#*******************************************************************************
# Verified compatible with both Python 2.7 and Python 3.6; written initially in 
# Python 3. 
#
#
# PURPOSE: Takes output from command line-based BLAST and brings it into 
# Python as a dataframe and saves a file of that dataframe for use elsewhere. 
# Optionally, it can also return that dataframe for use inside a Jupyter 
# notebook.
#
# This script is meant to be a utility script for working with command 
# line-based BLAST and Python, see a demonstration of use in
# https://git.io/vh8M9
# 
# It working necessitates that the BLAST command run with `-outfmt "6 qseqid 
# sseqid stitle pident qcovs length mismatch gapopen qstart qend sstart send 
# qframe sframe frames evalue bitscore qseq sseq"` as illustrated at
# https://medium.com/@auguste.dutcher/turn-blast-results-into-a-presence-absence-matrix-cc44429c814
# With https://blastedbio.blogspot.com/2014/11/column-headers-in-blast-tabular-and-csv.html
# listing definitiions for those text codes.
#
#
# Written to run from command line or imported into/pasted/loaded inside a 
# Jupyter notebook cell. When doing in Jupyter (or IPython, I believe) you can 
# skip the file save intermediate, see https://git.io/vh8Mi for these advanced 
# examples.
#
#
#
# (Aside: this was originally called `BLAST_to_df.py` because NCBI uses all caps 
# on the software name everywhere but while in a debugging nightmare and trying 
# many things, I found https://stackoverflow.com/a/23030897/8508004 and link to 
# documentation that says "Modules ...[and]... Packages should have short, 
# all-lowercase names". So changed name.
# But it turned out my debugging nightmare was casued
# by trying to debug the imported function and it seems Jupyter wasn't fully
# importing even though I reran importing command--> see 
# https://support.enthought.com/hc/en-us/articles/204469240-Jupyter-IPython-After-editing-a-module-changes-are-not-effective-without-kernel-restart
# for confirmation and solution. The issue was when I thought I was re-running
# the import command, I wasn't and the traceback was trying to show the same
# error but it uses the code file it has to show the traceback and that is why
# sometimes it seems to point to irrelevant place, such as in comment.)
# 
#
#
#
#
# Developed by adapting simple `BLAST_to_df()` function at 
# https://gist.github.com/fomightez/baf668acd4c51586deed2a2c89fcac67 to be 
# more full-featured using backbone of `patmatch_results_to_df.py`.
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
#
#
#
#
# TO RUN:
# Examples,
# Enter on the command line of your terminal, the line
#-----------------------------------
# python blast_to_df.py RESULTS_FILE
#-----------------------------------
# Issue `blast_to_df.py -h` for details.
# 
# More examples from running from the command line are at the links below: 
# https://git.io/vh8M9
#
#
#
# To use this after importing/pasting or loading into a cell in a Jupyter 
# notebook, specify at least the results file (or results as a string) in the 
# call to the main function similar to below:
# df = blast_to_df("test.out")
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
df = blast_to_df("test.out")
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

## Settings and options for output plot 
df_save_as_name = 'BLAST_pickled_df.pkl' # name for saving pickled dataframe

#
#*******************************************************************************
#**********************END USER ADJUSTABLE VARIABLES****************************













#*******************************************************************************
#*******************************************************************************
###DO NOT EDIT BELOW HERE - ENTER VALUES ABOVE###

import sys
import os



###---------------------------HELPER FUNCTIONS---------------------------------###

def _(row):
    '''
    takes the row and .... PLACEHOLDER FOR NOW
    '''
    return _

###--------------------------END OF HELPER FUNCTIONS---------------------------###
###--------------------------END OF HELPER FUNCTIONS---------------------------###

#*******************************************************************************
###------------------------'main' function of script---------------------------##

def blast_to_df(results, return_df = True, pickle_df=True):
    '''
    Main function of script. 
    BLAST results to Pandas dataframe.
    based on https://medium.com/@auguste.dutcher/turn-blast-results-into-a-presence-absence-matrix-cc44429c814

    It will take a file of results (or the results as a Python string) from 
    command line-based BLAST and make a dataframe that will be more useful 
    with Python/other genetic-oriented Python scripts.
    Optionally also returns a dataframe of the results data. 
    Meant for use in a Jupyter notebook.

    The option to provide the results as a string is to handle where sending the
    data directly from shell to Python without a typical file 
    intermediate, see the advanced notebook at https://git.io/vh8M9 for 
    examples. The obvious use case for that is when working in the Jupyter 
    notebook environment.

    For the function to work, it necessitates that the BLAST command be run with
    `-outfmt "6 qseqid sseqid stitle pident qcovs length mismatch gapopen qstart
    qend sstart send qframe sframe frames evalue bitscore qseq sseq"` as
    illustrated at 
    https://medium.com/@auguste.dutcher/turn-blast-results-into-a-presence-absence-matrix-cc44429c814
    '''
    import pandas as pd
    col_names = ['qseqid', 'sseqid', 'stitle', 'pident', 'qcovs', 'length',
    'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'qframe',
    'sframe', 'frames', 'evalue', 'bitscore', 'qseq', 'sseq']
    # Bring in the necessary data and mke collected results into dataframe:
    #---------------------------------------------------------------------------
    # df = pd.read_csv(results, sep='\t', header=None, names=col_names)
    # AT FIRST I HAD ABOVE, because I mistakenly though documentation said
    # `pd.read_csv()` work with a string so I thought I wouldn't have to resort
    # to acrobatics I used in
    # `patmatch_results_to_df.py` to handle correctly whether provided a file
    # or string; however, seems that I do need to handle all that.
    try:
        df = pd.read_csv(results, sep='\t', header=None, names=col_names)
    except (TypeError,OSError,IOError) as e:
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
        try:
            from StringIO import StringIO
        except ImportError:
            from io import StringIO
        df = pd.read_csv(StringIO(results), sep='\t', header=None, names=col_names)
        # I need StringIO so string handled as file document. 
        

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
        "because `no_pickling` was specified.")
    else:
        df.to_pickle(df_save_as_name )
        # Let user know
        sys.stderr.write("\n\nA dataframe of the data "
        "has been saved as a file\nin a manner where other "
        "Python programs can access it (pickled form).\n"
        "RESULTING DATAFRAME is stored as ==> '{}'".format(df_save_as_name ))

    
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
    # calling script from command line
    blast_to_df(results,**kwargs)
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
    parser = argparse.ArgumentParser(prog='blast_to_df.py',
        description="blast_to_df.py \
        Takes output from command line-based BLAST and brings it into \
        Python as a dataframe and saves a file of that dataframe for use \
        elsewhere. Optionally, it can also return that dataframe for use \
        inside a Jupyter notebook. Meant to be a utility script for working \
        with command line-based BLAST and Python.\
        The BLAST command needs to use `-outfmt '6...` as ilustrated at \
        https://medium.com/@auguste.dutcher/turn-blast-results-into-a-presence-absence-matrix-cc44429c814. \
        **** Script by Wayne Decatur   \
        (fomightez @ github) ***")

    parser.add_argument("results_file", help="Name of file of BLAST results \
        file to parse.\
        ", metavar="RESULTS_FILE")
    parser.add_argument('-dfo', '--df_output', action='store', type=str, 
    default= df_save_as_name, help="OPTIONAL: Set file name for saving pickled \
    dataframe. If none provided, '{}' will be used. To force no dataframe to \
    be saved, enter `-dfo no_pickling` without quotes as output file \
    (ATYPICAL).".format(df_save_as_name))



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
