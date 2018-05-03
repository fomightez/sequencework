#!/usr/bin/env python
# patmatch_results_to_df.py
__author__ = "Wayne Decatur" #fomightez on GitHub
__license__ = "MIT"
__version__ = "0.1.0"


# patmatch_results_to_df.py by Wayne Decatur
# ver 0.1
#
#*******************************************************************************
# Verified compatible with both Python 2.7 and Python 3.6; written initially in 
# Python 3. (See below.)
#
#
# PURPOSE: Takes output from command line-based PatMatch and brings it into 
# Python as a dataframe and saves a file of that dataframe for use elsewhere. 
# Optionally, it can also return that dataframe for use inside a Jupyter 
# notebook.
#
# This script is meant to be a utility script for working with command 
# line-based PatMatch and Python, see a demonstration of use in
# https://github.com/fomightez/patmatch-binder/blob/master/notebooks/PatMatch%20with%20Python.ipynb
# 
# Assumes for nucleic acid patterns, it was run with `-c` flag and tries to 
# assign strand information.
# 
# See https://github.com/fomightez/patmatch-binder about PatMatch.
#
# Written to run from command line or pasted/loaded inside a Jupyter notebook 
# cell. 
#
#
#
#
#
#
# Verified compatible with both Python 2.7 and Python 3.6; written initially in 
# Python 3. 
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
# - test on protein sequences
# - incorporate in demo notebook in patmatch-binder; also(?) in that demo binder 
#   show how to bring into Python the web-based PatMatch data from the xls file?
#   Also include running from command line (curl to get to rep after uploaded; this way will always get most current) and using main function.
#   Incorporate  running on entire genome in another notebook.
# - Advanced notebook for demo: `Advanced: Sending PatMatch output directly to Python`, under `Additional topics`?
#
#
#
#
# TO RUN:
# Examples,
# Enter on the command line of your terminal, the line
#-----------------------------------
# python patmatch_results_to_df.py -RESULTS_FILE
#-----------------------------------
# Issue `patmatch_results_to_df.py -h` for details.
# 
#
# When using in a notebook, if you don't specify dataframe objects, , you must
# instead supply strings of file names for the pickled dataframes in the call
# to the main function. 
# To use this after pasting or loading into a cell in a Jupyter notebook, in
# the next cell specify the two dataframes then call the main function similar 
# to below:
# my_pattern= "DDWDWTAWAAGTARTADDDD"
# df = patmatch_results_to_df("test.out", pattern=my_pattern, name="promoter")
# df
#
#
#
#
# 
#
'''
CURRENT ACTUAL CODE FOR RUNNING/TESTING IN A NOTEBOOK WHEN LOADED OR PASTED IN 
ANOTHER CELL:
my_pattern= "DDWDWTAWAAGTARTADDDD"
df = patmatch_results_to_df("test.out", pattern=my_pattern, name="promoter")
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
df_save_as_name = 'patmatch_pickled_df.pkl' # name for saving pickled dataframe

#
#*******************************************************************************
#**********************END USER ADJUSTABLE VARIABLES****************************













#*******************************************************************************
#*******************************************************************************
###DO NOT EDIT BELOW HERE - ENTER VALUES ABOVE###

import sys
import os
import pandas as pd


###---------------------------HELPER FUNCTIONS---------------------------------###

def fix_pattern_ids(row):
    '''
    takes the row and makes the unique pattern id numbering better match order
    after sorting dataframe. Uses the numerical index of the row that can be
    accessed as `row.name`, see https://stackoverflow.com/a/26658301/8508004.
    '''
    #+1 to make identifier number match common language and not be zero-indexed
    return row.hit_id.rsplit('-',1)[0] + '-' + str((row.name+1))

###--------------------------END OF HELPER FUNCTIONS---------------------------###
###--------------------------END OF HELPER FUNCTIONS---------------------------###

#*******************************************************************************
###------------------------'main' function of script---------------------------##

def patmatch_results_to_df(
    results_file, pattern="?", name=None, return_df = True, pickle_df=True):
    '''
    Main function of script. 
    It will take a file of results from command line-based PatMatch and make
    a dataframe that will be more useful with Python/othergenetic-oriented 
    scripts.
    Optionally also returns a dataframe of the results data. Meant for use in 
    a Jupyter notebook.
    '''
    # Bring in the necessary data:
    #---------------------------------------------------------------------------

    with open(results_file, 'r') as the_file:
        results = the_file.read()

    # feedback
    sys.stderr.write("Provided results read...")


    # Parse:
    #---------------------------------------------------------------------------
    results = results.split('>')
    # remove blanks
    results = [x for x in results if x]

    # prepare to give some unique indentifiers to each match
    identifiers=[]
    if name:
        id_prefix = name + "-"
    elif pattern == "?":
        id_prefix = "match-"
    elif len(pattern) > 29:
        id_prefix = pattern[:15] + "...-"
    else:
        id_prefix = pattern + "-"

    # prepare additional lists and variables for collecting data
    source_ids =[]
    last_hit_source = ''
    hit_nums_per_source = []
    from collections import defaultdict
    hit_num_tracking_dict = defaultdict(int) #based on approach used in 
    # `plot_sites_position_relative_genes.py` to work when entried interwearved
    matching_patterns = []
    starts = []
    ends = []
    strand_info = []
    for indx,each in enumerate(results):
        first_line, matching_pattern = [x.strip() for x in each.split()][:2]

        # Parse numbers between the brackets in first line.(In preparation, I 
        # will split on the ':' and use last one (index -1) just in case the 
        # extracted first line has `[]` in first part). Since want source 
        # identifier too, can get while setting up parse between brackets too.
        # Instead of just skipping to `second_half = first_line.split(":")[-1]`. 
        fline_split_colon = first_line.rsplit(":",1)
        source_id, second_half = fline_split_colon[0], fline_split_colon[1]
        nums_str = second_half[second_half.find("[")+1:second_half.find("]")]
        start,end = nums_str.split(",")[:2]

        # Increment appropriate hit number; handling via dictionary so okay if
        # interweaved which seems to happen when multiple fasta entries and
        # using the `-c` flag with PatMatch
        hit_num_tracking_dict[source_id]+= 1

        # Determine strand. 
        # I noticed when `-c` flag is used the first number in the interval 
        # returned to indicate the location will be larger than the second for 
        # those on the negative strand. 
        # Obviously this is a pointless exercise for protein sequence, but I 
        # think easier to leave and let it scan anyway and then fix pointless
        # column in one step later. This could be revisied if optimization is
        # necessary for scaling up.
        if start > end:
            strand = -1
            # fix start and end so start is actually lowest value to be 
            # consistent with system I have been using of late (like Ensembl)
            start, end = end,start

        else:
            strand = 1
        assert start < end ,"The 'start' value should be lower; strand \
        information is handled by `strand` property."
        identifiers.append(id_prefix+str(indx+1)) #`+1` so numbering more 
        # tpyical than python zero-indexing
        source_ids.append(source_id)
        hit_nums_per_source.append(hit_num_tracking_dict[source_id])
        matching_patterns.append(matching_pattern)
        starts.append(start)
        ends.append(end)
        strand_info.append(strand)

    # Make collected results into dataframe and improve on it
    #---------------------------------------------------------------------------
    df = pd.DataFrame(list(zip(
        source_ids,hit_nums_per_source,identifiers, starts, ends,strand_info, 
        matching_patterns)),
        columns=['FASTA_id','hit_number','hit_id', 'start','end','strand',
        'matching pattern'])
    # add query pattern as a column
    df['query pattern'] = pattern

    # better re-order the columns(?)

    # If multiple sequences involved, sort based on source FASTA and then hit #
    if len(source_ids) > 1:
        df.sort_values(['FASTA_id','hit_number'], ascending=[1,1], inplace=True)
        df = df.reset_index(drop=True) # fix index of dataframe
        # fix the unique order order to match order now. They are totally 
        # arbitrary and so perfectly fine to adjust still and otherwise they 
        # will just looks ugly after sorting
        df['hit_id'] = df.apply(fix_pattern_ids, axis=1)
        # Obviously, just a varation on that last line could have been used to 
        # make the identifiers in the first place or just number them all later. 
        # But more efficient to do way I did because that will always be the 
        # case and have all the infroamtion and thus reserves performing this 
        # step on the whole datframe again for only cases where several 
        # sequences scanned.

    # Remove irrelevant column if results from scan of protein sequence(s)
    if protein_results:
        df.drop('strand', axis=1)


    # Reporting and Saving
    #---------------------------------------------------------------------------
    #print(updated_sites_df)#originally for debugging during development,added..
    # Document the full set of data collected in the terminal or 
    # Jupyter notebook display in some manner. 
    # Using `df.to_string()` because more universal than `print(df)` 
    # or Jupyter's `display(df)`.
    sys.stderr.write( "\nFor documenting purposes, the following lists the "
        "parsed data:\n")
    #with pd.option_context('display.max_rows', None, 'display.max_columns', None):
    #    display(df)
    sys.stderr.write(df.to_string())

    # Handle pickling the modified sites dataframe
    if pickle_df == False:
        sys.stderr.write("\n\nA dataframe of the parsed data shown above "
        "was not stored for use\nelsewhere "
        "because `no_pickling` was specified in place of the output file name.")
    else:
        df.to_pickle(df_save_as_name )
        # Let user know
        sys.stderr.write( "\n\nA dataframe of the parsed data shown above "
        "has been\nsaved as a file in a manner where other "
        "Python programs\ncan access it (pickled form).\n"
        "RESULTING DATAFRAME is stored as ==> '{}'".format(df_save_as_name ))

    
    # Return dataframe (optional)
    #---------------------------------------------------------------------------
    if return_df:
        sys.stderr.write( "\n\nReturning a dataframe with the information "
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
    if args.pattern:
        kwargs['pattern'] = args.pattern
    else:
        kwargs['pattern'] = '?'
    if args.pattern_name:
        kwargs['name'] = args.pattern_name
    if df_save_as_name == 'no_pickling':
        kwargs['pickle_df'] = False
    kwargs['return_df'] = False #probably don't want dataframe returned if 
    # calling script from command line
    patmatch_results_to_df(results_file,**kwargs)
    # using https://www.saltycrane.com/blog/2008/01/how-to-use-args-and-kwargs-in-python/#calling-a-function
    # to build keyword arguments to pass to the function above
    # (see https://stackoverflow.com/a/28986876/8508004 and
    # https://stackoverflow.com/a/1496355/8508004 
    # (maybe https://stackoverflow.com/a/7437238/8508004 might help too) for 
    # related help)





if __name__ == "__main__" and '__file__' in globals():
    """ This is executed when run from the command line """
    # Code with just `if __name__ == "__main__":` alone will be run if pasted
    # into a notebook. The addition of ` and '__file__' in globals()` is based
    # on https://stackoverflow.com/a/22923872/8508004
    # See also https://stackoverflow.com/a/22424821/8508004 for an option to 
    # provide arguments when prototyping a full script in the notebook.
    ###-----------------for parsing command line arguments-----------------------###
    import argparse
    parser = argparse.ArgumentParser(prog='patmatch_results_to_df.py',
        description="patmatch_results_to_df.py \
        Takes output from command line-based PatMatch and brings it into \
        Python as a dataframe and saves a file of that dataframe for use \
        elsewhere. Optionally, it can also return that dataframe for use \
        inside a Jupyter notebook. Meant to be a utility script for working \
        with command line-based PatMatch and Python.\
        Assumes for nucleic acid patterns, it was run with `-c` flag and tries \
        to assign strand. \
        **** Script by Wayne Decatur   \
        (fomightez @ github) ***")

    parser.add_argument("results_file", help="Name of file of PatMatch results \
        file to parse.\
        ", metavar="RESULTS_FILE")
    parser.add_argument('-patt', '--pattern', action='store', type=str, 
        help="**OPTIONAL** Pattern used to perform the pattern matching search \
        that generated the results. The resulting dataframe will more \
        informative if one is provided; however, it is not essential. To \
        provide the pattern simply enter the text after the flag. For example, \
        if the search had been for an EcoRI site, include `--pattern GAATTC`, \
        without quotes or ticks, in the call to the script.")
    parser.add_argument('-name', '--pattern_name', action='store', type=str, 
        help="**OPTIONAL** Identifier for the provided pattern. If provided, \
        it will be used to make a unique indentifer for each match; however, \
        it is not essential as a unique id will still be assigned in its \
        absence. To provide the identifier simply enter the text after the \
        flag. For example, if the search had been for patternn `GAATTAC`, you \
        can include `-name EcoRI`, \
        without quotes or ticks, in the call to the script.")
    parser.add_argument("-p", "--protein_results",help=
    "add this flag to indicate the data are from a pattern match of protein \
    sequences. Otherwise it assumed the results are from pattern matching on \
    nucleic acid sequences.", action="store_true")
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
    results_file= args.results_file
    df_save_as_name = args.df_output


    main()

#*******************************************************************************
###-***********************END MAIN PORTION OF SCRIPT***********************-###
#*******************************************************************************
