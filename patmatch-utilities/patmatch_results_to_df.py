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
# https://git.io/vpo13
# 
# Assumes for nucleic acid patterns, it was run with `-c` flag and tries to 
# assign strand information.
# 
# See https://github.com/fomightez/patmatch-binder about PatMatch.
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
# Dependencies beyond the mostly standard libraries/modules:
#
#
#
# VERSION HISTORY:
# v.0.1. basic working version

#
# To do:
# - ?
#
#
#
#
# TO RUN:
# Examples,
# Enter on the command line of your terminal, the line
#-----------------------------------
# python patmatch_results_to_df.py RESULTS_FILE
#-----------------------------------
# Best to add the `-p` flag if the results are from protein sequences.
# Issue `patmatch_results_to_df.py -h` for details.
# 
# More examples from running from the command line are at the links below: 
# https://github.com/fomightez/sequencework/blob/master/patmatch-utilities/
# https://git.io/vpo13
# https://git.io/vpo1s
# https://git.io/vpo1Z
# https://git.io/vpr7i
#
#
#
# To use this after importing/pasting or loading into a cell in a Jupyter 
# notebook, specify at least the results file (or results as a string) in the 
# call to the main function similar to below:
# my_pattern= "DDWDWTAWAAGTARTADDDD"
# df = patmatch_results_to_df("test.out", pattern=my_pattern, name="promoter")
# df
#
# - or-, for protein, (& if you like to import it) the lines would be similar to:
# !curl -O https://raw.githubusercontent.com/fomightez/sequencework/master/patmatch-utilities/patmatch_results_to_df.py
# from patmatch_results_to_df import patmatch_results_to_df
# aa_pattern= "TYEETGLQGHPS"
# prot_df = patmatch_results_to_df("pro.out", pattern=aa_pattern, protein_results= True)
# prot_df
#
#
# A more in-depth series of examples of using this script within a notebook 
# without need to save file intermediates is found at:
# https://git.io/vpr7i
#
# 
#
'''
CURRENT ACTUAL CODE FOR RUNNING/TESTING IN A NOTEBOOK WHEN IMPORTED/LOADED OR 
PASTED IN ANOTHER CELL:
my_pattern= "DDWDWTAWAAGTARTADDDD"
df = patmatch_results_to_df("test.out", pattern=my_pattern, name="promoter")
df

aa_pattern= "TYEETGLQGHPS"
prot_df = patmatch_results_to_df("pro.out", pattern=aa_pattern, protein_results= True)
prot_df
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
    results, pattern="?", name=None, return_df = True, pickle_df=True, 
    protein_results=False):
    '''
    Main function of script. 
    It will take a file of results (or the results as a Python string) from 
    command line-based PatMatch and make a dataframe that will be more useful 
    with Python/other genetic-oriented Python scripts.
    Optionally also returns a dataframe of the results data. Meant for use in a 
    Jupyter notebook.

    The option to provide the results as a string is to handle where sending the
    data directly from shell script to Python without a typical file 
    intermediate, see the advanced notebook at https://git.io/vpr7i for 
    examples. The obvious use case for that is when working in the Jupyter 
    # notebook environment.
    '''
    # Bring in the necessary data:
    #---------------------------------------------------------------------------

    try:
        with open(results, 'r') as the_file:
            results = the_file.read()
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
        df = df.drop('strand', axis=1)


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
        "because `no_pickling` was specified.")
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
    if args.protein_results:
        kwargs['protein_results'] = True
    if df_save_as_name == 'no_pickling':
        kwargs['pickle_df'] = False
    kwargs['return_df'] = False #probably don't want dataframe returned if 
    # calling script from command line
    patmatch_results_to_df(results,**kwargs)
    # using https://www.saltycrane.com/blog/2008/01/how-to-use-args-and-kwargs-in-python/#calling-a-function
    # to build keyword arguments to pass to the function above
    # (see https://stackoverflow.com/a/28986876/8508004 and
    # https://stackoverflow.com/a/1496355/8508004 
    # (maybe https://stackoverflow.com/a/7437238/8508004 might help too) for 
    # related help). Made it easy to add more as I thought of them.





if __name__ == "__main__" and '__file__' in globals():
    """ This is executed when run from the command line """
    # Code with just `if __name__ == "__main__":` alone will be run if pasted
    # into a notebook. The addition of ` and '__file__' in globals()` is based
    # on https://stackoverflow.com/a/22923872/8508004
    # See also https://stackoverflow.com/a/22424821/8508004 for an option to 
    # provide arguments when prototyping a full script in the notebook.
    ###-----------------for parsing command line arguments-------------------###
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
        flag. For example, if the search had been for patternn `GAATTC`, you \
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
    results = args.results_file
    df_save_as_name = args.df_output


    main()

#*******************************************************************************
###-***********************END MAIN PORTION OF SCRIPT***********************-###
#*******************************************************************************
