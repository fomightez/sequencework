#!/usr/bin/env python
# hhsuite3_results_to_df.py
__author__ = "Wayne Decatur" #fomightez on GitHub
__license__ = "MIT"
__version__ = "0.1.0"


# hhsuite3_results_to_df.py by Wayne Decatur
# ver 0.1
#
#*******************************************************************************
# Verified compatible with both Python 2.7 and Python 3.8; written initially in 
# Python 3. 
#
#
# PURPOSE: Takes output from HH-suite3 programs, such as HHblits and HHsearch, 
# in `.hhr` file format and brings it into Python as a dataframe and saves a 
# file of that dataframe for use elsewhere. 
# Optionally, it can also return that dataframe for use inside a Jupyter 
# notebook.
#
# This script is meant to be a utility script for working with command 
# line-based HH-suite3 and Python, see a demonstration of use in
# https://git.io/XXXX
# 
#
#
# Written to run from command line or imported into/pasted/loaded inside a 
# Jupyter notebook cell. When doing in Jupyter (or IPython, I believe) you can 
# skip the file save intermediate, see https://git.io/XXXX for these advanced 
# examples.
#
#
#
# 
#
#
#
#
# Developed by adapting backbone of `blast_to_df.py` to handle the `hhr` file 
# format gnerated by hhsearch and hhblits as part of HH-suite3. First, converted
# `.hhr` file to something resembling the table that is generated with `-outfmt`
# option 6 settings from blast that produce a table. Then I can edit the headers
# / column names in `col_names` to match whatever I customized beyond those in 
# the `blast_to_df.py` list.
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
# - make sure works with Python 2.7 <- sometimes I was lazy during development 
# used f-strings, & those need to be replaced with string formatting using 
# `.format()` because Python 2.7 never had f-strings, unless I add the use of 
# future_fstrings package, see https://stackoverflow.com/a/46182112/8508004
# - in the internal documentation and in the readme update the https://git.io/XXXX link to point at hhsuite3-binder` repo
#
#
#
#
# TO RUN:
# Examples,
# Enter on the command line of your terminal, the line
#-----------------------------------
# python hhsuite3_results_to_df.py RESULTS_FILE
#-----------------------------------
# Issue `hhsuite3_results_to_df.py -h` for details.
# 
# More examples from running from the command line are at the links below: 
# https://git.io/XXXX
#
#
#
# To use this after importing/pasting or loading into a cell in a Jupyter 
# notebook, specify at least the results file (or results as a string) in the 
# call to the main function similar to below:
# df = hhsuite3_results_to_df("test.out")
# df
#
#
# A more in-depth series of examples of using this script within a notebook 
# without need to save file intermediates is found at:
# https://git.io/XXXX
#
# 
#
'''
CURRENT ACTUAL CODE FOR RUNNING/TESTING IN A NOTEBOOK WHEN IMPORTED/LOADED OR 
PASTED IN ANOTHER CELL:
df = hhsuite3_results_to_df("results.hhr")
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

## Settings and options for output dataframe-as-file 
df_save_as_name = 'hhs3_results_pickled_df.pkl' #name for pickled dataframe file

#
#*******************************************************************************
#**********************END USER ADJUSTABLE VARIABLES****************************













#*******************************************************************************
#*******************************************************************************
###DO NOT EDIT BELOW HERE - ENTER VALUES ABOVE###

import sys
import os
import pandas as pd
# I need StringIO so string handled as file document. Also need to deal 
# with whether Python 3 or 2 because StringIO source differs for Python 2.
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO


###---------------------------HELPER FUNCTIONS---------------------------------###

def _(row):
    '''
    takes the row and .... PLACEHOLDER FOR NOW
    '''
    return _

def split_out_len(text_str):
    '''
    Takes a text string and splits out the number value from in between the 
    parantheses at the end.

    Returns in integer value

    Example input string:
    41-281 (291)

    Example return:
    291
    '''
    return int(text_str.split("(")[1].split(")")[0])


def tablicize_hhr(hhr_file_content):
    '''
    Takes HH-suite3 results in `.hhr` format and returns a string of the results
    arranged in a tab-delimited text table that can easily be read into Pandas.

    I am going to add some derived data as well. For now, the derived data to 
    be added will be the change in total size of the protein between the query 
    and each hit. 
    '''
    # Divide into three portions:
    #---------------------------------------------------------------------------
    # Structurally, the `.hhr` format text looks to break down to the following
    # three parts:
    # 1. top header
    # 2. simplified hit table
    # 3. Section with details for each hit.
    # The header is going to the source of two things:
    # a. The query identifier/description
    # b. size of the query protein as `Match_columns` value seems to correspond 
    # to that 
    # The simplified hit table will mostly be used to have a list of all the 
    # hits to mine out of the details section. The structure of the simplified 
    # hit table will also be used to collect the length of the sequence in the 
    # hit as the value and in parantheses at the end seems to correspon to that.
    # The simplified hit table contains many values redundant with those in the 
    # details section. Some of these redundant values could be used to verify
    # information; however, I haven't decided how much validation I want to 
    # perform using that information.
    # Most of the data will be mined from the third section as the names and
    # descriptions are not cut off and additional valuable data not present in 
    # simplified table are also there such as `identity` and `similarity`.
    simplified_hit_table_delimiter = " No Hit "
    details_section_delimiter = "No 1"

    header_text = hhr_file_content.split(simplified_hit_table_delimiter,1)[0]
    simpl_hit_table = simplified_hit_table_delimiter + hhr_file_content.split(
        simplified_hit_table_delimiter,1)[1].split(
        details_section_delimiter,1)[0]
    # the `simpl_hit_table` is an example of a table of fixed-width formatted 
    # lines and it can be read into a Pandas DataFrame using `pandas.read_fwf`
    # data and colspecs example of a table starting with Keanu Reeves can be 
    # found 
    # https://github.com/birforce/vnpy_crypto/blob/b9bb23bb3302bf5ba47752e93f8b23e04a9a2b27/venv/lib/python3.6/site-packages/pandas/tests/io/parser/test_read_fwf.py#L279
    colspecs = ((0, 3), (4, 34), (35, 40), (41, 48), (49, 57), (58,64), (65,69), 
        (70,74), (75,84), (85,100))
    # `pd.read_fwf()` step to be done next should have been the following but 
    # `skip_blank_lines = True` didn't seem to work and couldn't cast `Nan` to 
    # integer in reading in of the table then and so did as separate steps 
    # after:
    '''
    simp_df = pd.read_fwf(StringIO(simpl_hit_table), colspecs=colspecs, 
        skip_blank_lines = True, dtype={'No': 'int64','Cols': 'int64'})
    '''
    simp_df = pd.read_fwf(StringIO(simpl_hit_table), colspecs=colspecs)
    # Drop rows of all Nan or convert to integer in next step won't work
    simp_df = simp_df.dropna(how='all')  # should no longer needed after adding 
    # `skip_blank_lines = True` but didn't work in `pd.read_fwf()` for some 
    # reason
    # Fix `No` column back to integer by casting it
    simp_df.No = simp_df.No.astype(dtype='int64') # similarly, if 
    # `skip_blank_lines = True` had worked in `pd.read_fwf()` could use 
    # `dtype={'No': 'int64'}` in `pd.read_fwf()`, but didn't for an unknown 
    # reason and so easier to do separate to insude it happens for now.
    # Also fix `Cols` column back to integer by casting it
    simp_df.Cols = simp_df.Cols.astype(dtype='int64') 

    # Then extract hit length value that is in parantheses at the end of the 
    # last column out into it's own column.
    simp_df['hit_total_length'] = simp_df['Template HMM'].apply(split_out_len)

    hit_details_section_text = (
        details_section_delimiter + hhr_file_content.split(
        details_section_delimiter,1)[1])
    # make a list of metrics in each record in the details section
    metrics_list = hit_details_section_text.split("\n")[2].split()
    metrics = [x.split("=")[0] for x in metrics_list]
    hit_dict = {} # intialize a dictionary to begin to collect information for
    # each hit. The keys will be the hit number from the `No` column of the
    # simplified hits detail text table. Values will be a dictionary with the
    # following keys and values:
    # `qid`: query id
    # `qtitle` : query title
    # 'query_length' : query length
    # `hid` : hit id
    # `htitle` : hit title
    # 'hit_section' : text of extracted section for that hit record
    # 'hit_length' : total length in amino acids of that hit
    # 
    generalized_labels_from_alignment_section = [] # intialize a list to collect
    # all the generlized labels collected as a result of the alignment section.
    # I assume all these will be the same for each hit but this allows for that
    # not to be the case by collecting them all and then filtering down to 
    # unique ones when go to make header labels. This allows variation I see in
    # what different hhr files have for these labels to be handled more 
    # generally and not be hard coded, or have to account for any variant.
    # Use the simplified hit dataframe `No` (number) column to collect the 
    # specific sections for each hit from the details section at bottom of hhr
    num_col_as_list = simp_df['No'].tolist()
    for indx,row in enumerate(simp_df.itertuples()):
        hit_dict[row.No] = {}
        start_delimiter = 'No ' + str(row.No)
        if indx == len(simp_df) - 1:
            hit_section = (
                start_delimiter + hit_details_section_text.split(
                start_delimiter,1)[1])
        else:
            hit_section = (
                start_delimiter + hit_details_section_text.split(
                start_delimiter,1)[1]).split(
                'No ' + str(num_col_as_list[indx + 1]))[0]
        qid = (
            header_text.split("Query")[1].split("\n",1)[0].strip().split()[0])
        hit_dict[row.No]['qid'] = qid
        hit_dict[row.No]['qtitle'] = (
            header_text.split("Query")[1].split("\n",1)[0].strip())
        hit_dict[row.No]['query_length'] = (int(
            header_text.split("Match_columns ")[1].split("\n",1)[0].strip()))
        hit_dict[row.No]['hit_section'] = hit_section
        hit_dict[row.No]['hit_length'] = row.hit_total_length
        hid = hit_section.split(">")[1].split()[0]
        hit_dict[row.No]['hid'] = hid
        hit_dict[row.No]['htitle'] = (
            hit_section.split(">")[1].split('Probab')[0].strip())
        # interate over the extracted text and record the value for each metric
        for metric in metrics: 
            hit_dict[row.No][metric] = (
                hit_section.split("{}=".format(metric))[1].split()[0])

        

        # Extract the alignment section from the end of the text for the details
        # of the hit
        pw_aln_section = hit_section.split(metrics[-1])[1].split("\n",1)[1]
        

        # split the pairwise HMM alignments into chunks so can build each line 
        # as single strings, even if long
        pw_chunks = pw_aln_section.split("\n\n")
        # discard chunks that are empty, in other words just `'\n'`
        pw_chunks = [chunk for chunk in pw_chunks if chunk != '\n']
        # pw_chunks[0] = pw_chunks[0][1:] # remove the new line symbol right at 
        # start of first line of the chunk so that detecting the line in the 
        # middle of each chunk that doesn't have a label and corresponds to 
        # `pw_conservation` label below is easily managed. ACTUALLY NEED TO DO
        # FOR ALL FIRST LINES IF MATCH, next line does that:
        pw_chunks = [chunk[1:] for chunk in pw_chunks if chunk.startswith('\n')]
        # First, go through the first chunk of the pairwise alignment and make a 
        # list of the labels.
        #----------------------------------------------------------------------#
        # These will be used to collect the text from multiple chunks and then 
        # to later store the collected part of the alignment with the 
        # appropriate key.
        # the collected text once collected
        # Note that if the secondary structure wasn't added to the query or the 
        # values in the database, itseems the alignmnet section for each record 
        # will lack `Q ss_pred` & `T ss_pred` lines. So determine if that is the 
        # case and parse the alignment section appropriately. Some examples of 
        # this are the hhr files from the hits against the custom database.
        # Also the `Confidence` line can also be optional and IMPORTANTLY it can
        # be avilable for top ranked hits and not in the ones further below!!!
        # For these three, determine if they are present and adjust the key
        # values and collection of data appropriately.
        # The `alignment_section_label_relator` will assist in this process by
        # determining what labels the pairwise HMM alignments section has that 
        # correspond to what general categories of output labels. For example 
        # `Q Consensus` will be `qconsensus` in the output dataframe later. That
        # is reminiscent of some of the BLAST labels used exept I am using `h` 
        # for `hit` here instead of `s` for `subject`.
        alignment_section_label_relator = {
            'Q Consensus':'qconsensus',
            'T Consensus':'hconsensus',
            'Confidence':'aln_confidence',
            f'Q {qid}':'qseq',
            f'T {hid}':'hseq',
            'Q ss_pred':'qss_pred',
            'T ss_pred':'hss_pred',
            }
        # In each pairwise HMM alignment chunk there will be a line in the
        # middle that doesn't have a label, this will be the column score 
        # between the query and template amino acid distributions, and I'll call 
        # it the `pw_conservation` for `pairwise conservation`
        alignment_section_label_relator['pw_conservation'] = 'pw_conservation'
        alignment_section_labels = []
        for line in pw_chunks[0].split("\n"):
            # this `if line[:16].strip()` will skip line with no label in middle
            if not line[:16].strip():
                alignment_section_labels.append('pw_conservation')
            else:
                alignment_section_labels.append(line[:16].strip())


        # Now iterate over each chunk of the pairwise alignment and build each
        # long string from the individual lines
        #----------------------------------------------------------------------#
        # First initialize strings for each label in the chunks. Use a 
        # dictonary to keep them.
        aln_strings = {}
        for labl in alignment_section_labels:
            aln_strings[labl] = ""
        # Going to iterate through each chunk assigning content in columns 
        # 22 - 102 to the appropriate label, which will be kept in 
        # correspondence by using enumeration to get the correct label from the
        # alignment_section_labels. For this to work, each chunk has to have the
        # same number of lines as alignment_section_labels. Check that first.
        # all(len(i) == len(myList[0]) for i in myList)
        assert all(
            len(chunk.split("\n")) == len(
            alignment_section_labels) for chunk in pw_chunks),(
            "Chunks of alignments need to have same number of lines as labels "
            "detected.")
        for chunk in pw_chunks:
            for indx,line in enumerate(chunk.split("\n")):
                aln_strings[alignment_section_labels[indx]] += (
                    line[22:103].split()[0]) # The split is in order to get the 
                    # just the sequence for the end lines where it is shorter 
                    # width than ~100

        # Now with the entire lines collected assign these to entries in 
        # hit_dict. Use the `alignment_section_label_relator` to add in the 
        # proper but generalized keys. In other words, make it so these keys 
        # are the same for each hit and will work well in a table later. The 
        # labels collected from the chunk will have a unique id for the query
        # and template but want something more general. And might as will make
        # the general labels more BlAST like, see the larger comment above 
        # defining `alignment_section_label_relator`.
        for line_label,collected_str in aln_strings.items():
            generalized_label = alignment_section_label_relator[line_label]
            generalized_labels_from_alignment_section.append(generalized_label)
            hit_dict[row.No][generalized_label] = collected_str

    

    # Calculate and add derived data:
    #---------------------------------------------------------------------------
    # I am going to add some derived data as well. For now, the derived data to 
    # be added will be the change in total size of the protein between the query 
    # and each hit.
    for row_no, hit_details in hit_dict.items():
        hit_dict[row_no]['size_diff'] = abs(
            hit_dict[row_no]['query_length'] - hit_dict[row_no]['hit_length'])


    # Rearrange the hit data for each into a useful text table:
    #---------------------------------------------------------------------------
    # I am going to add some derived data as well. For now, the derived data to 
    # be added will be the change in total size of the protein between the query 
    # and each hit.
    # Add the labels on a header line that will then later get used by Pandas to
    # make collumn names. This makes things easier because I can use 
    # `generalized_labels_from_alignment_section` to easily handle those labels 
    # since I've seen them differ from different `hhr` files, for exmaple 
    # depending on if there are `Q ss_pred`/`T ss_pred` lines or `Confidence` 
    # lines in the pairwise alignments section.
    # Remember that
    columns_to_make_after_hit_num = (
        ['qid','qtitle','query_length','hid','htitle','hit_length'] + 
        metrics + ['size_diff'] + 
        list(set(generalized_labels_from_alignment_section)))
    header_line = "hit_num\t" + "\t".join(columns_to_make_after_hit_num)
    table_str = header_line + "\n"
    for row_no, hit_details in hit_dict.items():
        table_str += "{}\t".format(row_no)
        for l in columns_to_make_after_hit_num:
            # some well-ranked hits will have Confidence lines in the pairwise
            # alignments while poorly ranked hits may not, so just assign 
            # `np.nan` if the key isn't present not.
            try:
                table_str += "{}\t".format(hit_dict[row_no][l])
            except KeyError:
                if l == 'aln_confidence':
                    import numpy as np
                    table_str += "{}\t".format(np.nan)
                else:
                    raise
        table_str += "\n"
    return table_str

###--------------------------END OF HELPER FUNCTIONS---------------------------###
###--------------------------END OF HELPER FUNCTIONS---------------------------###

#*******************************************************************************
###------------------------'main' function of script---------------------------##

def hhsuite3_results_to_df(results, return_df = True, pickle_df=True):
    '''
    Main function of script. 
    HH-suite3 results to Pandas dataframe.
    It will take a file of results (or the results as a Python string) from 
    Hh-suite 3 and make a dataframe that will be more useful 
    with Python/other bioinformatics-oriented Python scripts.
    Optionally also returns a dataframe of the results data. 
    Meant for use in a Jupyter notebook.

    The option to provide the results as a string is to handle where sending the
    data directly from shell to Python without a typical file 
    intermediate, see the advanced notebook at https://git.io/XXXX for 
    examples. The obvious use case for that is when working in the Jupyter 
    notebook environment.

    Adapted from the main function in `blast_to_df.py`
    '''
    # Bring in the necessary data and prepare it in a manner for next step:
    #---------------------------------------------------------------------------
    # To make this easier, I'll try to follow moving the data into something 
    # along the lines of what BLAST does when you make a table using option 6, 
    # i.e., run with 
    # `-outfmt "6 qseqid sseqid stitle pident qcovs length mismatch gapopen 
    # qstart qend sstart send qframe sframe frames evalue bitscore qseq sseq"` 
    # as illustrated at 
    # https://medium.com/@auguste.dutcher/turn-blast-results-into-a-presence-absence-matrix-cc44429c814
    # I'll make it different in that I'll add column names as a header in 
    # the table made by tablicize_hhr()
    with open(results, 'r') as file_handle:
        results_tablicized_str = tablicize_hhr(file_handle.read())

    '''
    FROM BLAST:
    col_names = ['qseqid', 'sseqid', 'stitle', 'pident', 'qcovs', 'length',
    'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'qframe',
    'sframe', 'frames', 'evalue', 'bitscore', 'qseq', 'sseq']
    '''
    # Decided to pass the column names from the data converted into a data table
    # to Pandas using columns. BLAST didn't add labels/header line so that was
    # not possible there. This is easier, especially since some labels from
    # hhr pairwise alignment lines is optional. For example:
    # Keep in mind that `Q ss_pred` & `T ss_pred` are optional likes and so 
    # check if they are present and adjust the column names appropriately.
    # Also the `confidence` line can also be optional.  <=== This is now handled 
    # where I make the table here by adding headers. It will be easier than 
    # figuring out labels/headers in `tablicize_hhr()` and then figuring out 
    # here.
    #  

    # Make collected results table into dataframe:
    #---------------------------------------------------------------------------
    
    df = pd.read_csv(StringIO(
        results_tablicized_str), sep='\t', header=0, index_col=False)
        

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
    hhsuite3_results_to_df(results,**kwargs)
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
    parser = argparse.ArgumentParser(prog='hhsuite3_results_to_df.py',
        description="hhsuite3_results_to_df.py \
        Takes output from HH-suite3 programs, such as HHblits and HHsearch, \
        in `.hhr` file format and brings it into Python as a dataframe and \
        saves a file of that dataframe for use elsewhere. Optionally, it can \
        also return that dataframe for use inside a Jupyter notebook. \
        Meant to be a utility script for working \
        with HH-suite3 and Python.\
        **** Script by Wayne Decatur   \
        (fomightez @ github) ***")
    parser.add_argument("results_file", help="Name of file of .hhr results \
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
