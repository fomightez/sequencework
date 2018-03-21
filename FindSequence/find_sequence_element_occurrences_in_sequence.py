#!/usr/bin/env python
# find_sequence_element_occurrences_in_sequence.py
__author__ = "Wayne Decatur" #fomightez on GitHub
__license__ = "MIT"
__version__ = "0.1.0"


# find_sequence_element_occurrences_in_sequence.py by Wayne Decatur
# ver 0.1
#
#*******************************************************************************
# Verified compatible with both Python 2.7 and Python 3.6; written initially in 
# Python 3.
#
# PURPOSE: Takes a short sequence represented as a string and a FASTA-formatted
# sequence (either a file or URL to a file) and makes an accounting of the 
# occurrences of that short sequence element on both strands of the main 
# sequence. Exact matches only are counted as occurrences.
# The script is case-insensitive, meaning the case of either the element or 
# FASTA sequence provided do not matter and matches will be noted anyway.
#
# As written, the script expects and only processes one FASTA-formatted 
# sequence. If your, FASTA file has more than one sequence entry within it and 
# the one you want to have scanned is not the first, copy and paste it to a new 
# file and use that as the sequence to scan.
#
# Written to run from command line or pasted/loaded inside a Jupyter notebook 
# cell or imported. 
#
#
#
# This script based on work and musings developed in 
# `Resources for writing code to do GC cluster accounting.md`
#
# The script at https://www.biostars.org/p/209383/ (steve's answer) served as 
# the backbone for adding convenience and user-friendly features.
#
#
#
# Dependencies beyond the mostly standard libraries/modules:
# - Biopython
#
#
# VERSION HISTORY:
# v.0.1. basic working version
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
# python find_sequence_element_occurrences_in_sequence.py GAATTC sequences.fasta

#-OR-
# python find_sequence_element_occurrences_in_sequence.py GAATTC sequences.fasta URL
#-----------------------------------
# Issue `python find_sequence_element_occurrences_in_sequence.py -h` for details.
# 
#
# To use this after pasting or loading into a cell in a Jupyter notebook, in
# the next cell define critical variables and then call the main function 
# similar to below:
# source = "https://downloads.yeastgenome.org/sequence/S288C_reference/chromosomes/fasta/chrmt.fsa"
# element = "GAATTC"
# id_of_seq_element = "ele1" #Set to `None` without quotes or backticks to have defined automatically
# find_sequence_element_occurrences_in_sequence.py()
#
# Something similar would need to be done if importing the script into another.
# (`id_of_seq_scanned_hardcoded` can be assigned in a cell before calling the 
# function as well.)
#
# 
#
'''
CURRENT ACTUAL CODE FOR RUNNING/TESTING IN A NOTEBOOK WHEN LOADED OR PASTED IN 
ANOTHER CELL:
source = "https://downloads.yeastgenome.org/sequence/S288C_reference/chromosomes/fasta/chrmt.fsa"
element = "GAATTC"
id_of_seq_element = "EcoRI" #Set to `None` without quotes or backticks to have defined automatically
df = find_sequence_element_occurrences_in_sequence(return_dataframe = True)


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
output_file_name_prefix = "seq_account"
id_of_seq_scanned_hardcoded = None # replace `None` with what you want to use,
# with flanking quotes if something appropriate is not being extracted from the
# provided filepath/filename or URL to save as an indicator of the file scanned
# in the output file name. 



limit_of_name = 17 # number of bases of the sequence element to limit to using 
# if the sequence element sequence is used to make the name for the output file





#
#*******************************************************************************
#**********************END USER ADJUSTABLE VARIABLES****************************


















#*******************************************************************************
#*******************************************************************************
###DO NOT EDIT BELOW HERE - ENTER VALUES ABOVE###

import sys
import os
import re
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq # for reverse complement




###---------------------------HELPER FUNCTIONS---------------------------------###


def get_seq_element_representation(id_of_seq_element):
    '''
    Takes `id_of_seq_element` and returns a string represenation for using
    in the output file name and the dataframe.
    '''
    if id_of_seq_element:
        elem_id = id_of_seq_element
    elif len(element) > limit_of_name:
        elem_id = element[:limit_of_name]+"..."
    else:
        elem_id = element
    return elem_id



def generate_output_file_name(id_of_seq_element,id_of_seq_scanned):
    '''
    Takes a file name as an argument and returns string for the name of the
    output file. The generated name is based on a prefix that can be adjusted
    under the section ' USER ADJUSTABLE VALUES ', plus the provided text in 
    the function call.
    If there is no `id_of_seq_element` specified, the sequence element will be
    used, limited to the first number of bases specified in `limit_of_name`.

    Specific examples
    =================
    Calling function with
        ("elem1","chrmt")
    returns
        "seq_account_elem1_chrmt.tsv"

    Calling function with
        (None,"chrmt")
    returns
        "seq_account_GAATTC_chrmt.tsv"
    if `GAATTC` happened to be the provided sequence element.
    '''
    elem_id = get_seq_element_representation(id_of_seq_element)
    if elem_id.endswith("...") and (id_of_seq_element == None):
        # Because in order to use `get_seq_element_representation()` in all the 
        # places something similar is needed, I have it adding `...` if the
        # sequence of the element exceeds the limit size. Since I don't want
        # that in the file name I am going to remove if it seems I had added it,
        # and I will only have possibly added it if `id_of_seq_element` is None.
        # This extra check is mitigate chances I'll remove `...` in unlikely
        # possibility user wanted to include it in id for sequence element. 
        elem_id = elem_id[:-3] + "-" #add hyphen to indicate there is more 
    return "{prefix}_{elem_id}_{seqid}.tsv".format(
        prefix=output_file_name_prefix,elem_id=elem_id,
        seqid=id_of_seq_scanned)


def extract_id_of_seq_scanned(source):
    '''
    Take something like:
    https://downloads.yeastgenome.org/sequence/S288C_reference/chromosomes/fasta/chrmt.fsa

    -or- 

    /local/directory1/directory2/chrmt.fa

    -or-

    chrmt.fa

    And return:
    chrmt
    '''
    if "/" in source:
        last_bit = source.split("/")[-1]
    else:
        last_bit = source
    if '.' in last_bit:
        main_part_of_name, file_extension = os.path.splitext(
        last_bit) #from http://stackoverflow.com/questions/541390/extracting-extension-from-filename-in-python
        return main_part_of_name
    else:
        return last_bit


def get_seq_from_URL(url):
    '''
    takes a URL and gets the sequence
    '''
    try:
        from StringIO import StringIO
    except ImportError:
        from io import StringIO

    chromosomes_and_length = {}
    # Getting html originally for just Python 3, adapted from 
    # https://stackoverflow.com/a/17510727/8508004 and then updated from to 
    # handle Python 2 and 3 according to same link.
    try:
        # For Python 3.0 and later
        from urllib.request import urlopen
    except ImportError:
        # Fall back to Python 2's urllib2
        from urllib2 import urlopen
    html = urlopen(url)
    fasta_iterator = SeqIO.parse(StringIO(
        html.read().decode(encoding='UTF-8')), "fasta")
    # Use of `next()` on next line to get first FASTA -formatted sequence is 
    # based on http://biopython.org/DIST/docs/api/Bio.SeqIO-module.html
    # I think difference from `SeqIO.read()` in this approach is that it won't
    # give an error of more than one entry is in the html.
    # I found I needed `StringIO()` or got issues with trying to handle long file name.
    record = next(fasta_iterator)
    return record.seq


def search_strand(pattern, sequence_to_scan, strand=1):
    '''
    take a sequence pattern (element) and find occurrences of that on the 
    provided, larger 5'-->3' sequence.
    Assumes strand is first unless provided.

    Tracks the start and end points of each occurrence, returning a list of
    that information where each element is a tuple of the start and end points
    along with the strand.

    based on https://www.biostars.org/p/209383/ (specifically steve's answer)
    '''
    occurrences = []
    for match in re.finditer(pattern.upper(), str(sequence_to_scan.upper())):
        if strand == 1:
            start_pos = match.start() + 1
            end_pos = match.end() + 1
        else:
            start_pos = (len(sequence_to_scan) - match.start() ) + 1
            end_pos = (len(sequence_to_scan) - match.end() ) + 1
        # print (start_pos, '\t', end_pos, '\t',strand) # for debugging
        occurrences.append((start_pos, end_pos,strand))
    return occurrences 


def find_sequence_element_occurrences_in_sequence(return_dataframe = False):
    '''
    Main function of script. Scan a sequence and report on occurrences of a
    sub-sequence element in that sequence.

    Returns None
    Unless `return_dataframe = True`, and then it returns a dataframe of 
    accounting as well. That option being meant when using this script in a cell
    in a Jupyter notebook or importing it into another script.
    '''
    # get the fasta_seq to scan
    if source.lower().startswith("http"):
        fasta_seq = get_seq_from_URL(source)
    else:
        # Read sequence, treating source as a filepath.
        # Use of `with` on next line based on http://biopython.org/wiki/SeqIO , 
        # under "Sequence Input". Otherwise, backbone based on 
        # https://www.biostars.org/p/209383/, and fact `rU` mode depecated.
        with open(source, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                # print(record.seq) # for debugging
                fasta_seq = record.seq

    # With the approach in this next block, I can expose `id_of_seq_scanned` to 
    # setting for advanced use without it being required and without need to be 
    # passed into the function.
    if id_of_seq_scanned_hardcoded:
        id_of_seq_scanned = id_of_seq_scanned_hardcoded
    else:
        id_of_seq_scanned = extract_id_of_seq_scanned(source)
        


    #assert that element cannot be longer than fasta_seq
    assert len(element) < len(fasta_seq), (
    "the FASTA sequence has to be longer than the provided sequence element.\n"
    "The provided FASTA sequence, {0}, is {1} bases;\nthe provided sequence "
    "element '{2}' is {3} bases.".format(
        id_of_seq_scanned,len(fasta.seq),element,len(element)))

    occurrences = search_strand(element, fasta_seq) # first strand
    occurrences += search_strand(
        element, fasta_seq.reverse_complement(), strand=-1) #2nd strand
    if occurrences:
        # make it into a dataframe since provides convenient options for 
        # handling
        df = pd.DataFrame(occurrences, columns=['start_pos','end_pos','strand'])
        # add some useful information to the dataframe
        df["seq. element"] = get_seq_element_representation(id_of_seq_element)
        df = df[['seq. element','start_pos','end_pos','strand']] # 'seq.element'
        # column will show on far right otherwise

        # write to tab-delimited file
        output_file_name = generate_output_file_name(
            id_of_seq_element,id_of_seq_scanned )
        df.to_csv(output_file_name, sep='\t',index = False)
        sys.stderr.write( "\nThe information on the {0} occurrences of '{1}' "
                "has\nbeen saved as a file named"
                " '{2}'.".format(len(df),element,output_file_name))

        if return_dataframe:
            sys.stderr.write( "\n\nReturning a dataframe with the information "
                "as well.")
            return df
    else:
        sys.stderr.write( "\nNo occurrences of '{0}' "
                "found in the provided sequence.".format(element))
        if return_dataframe:
            sys.stderr.write( "\n\nNo data to return in a dataframe and so "
                "returning `None`.")
            return None


    










###--------------------------END OF HELPER FUNCTIONS---------------------------###
###--------------------------END OF HELPER FUNCTIONS---------------------------###















#*******************************************************************************
###------------------------'main' secion of script---------------------------###

def main():
    """ Main entry point of the script """
    # placing actual main action in a 'helper' script so can call that easily 
    # with a distinguishing name in Jupyter notebooks, where `main()` may get
    # assigned multiple times depending how many scripts imported/pasted in.
    find_sequence_element_occurrences_in_sequence()
        







if __name__ == "__main__" and '__file__' in globals():
    """ This is executed when run from the command line """
    # Code with just `if __name__ == "__main__":` alone will be run if pasted
    # into a notebook. The addition of ` and '__file__' in globals()` is based
    # on https://stackoverflow.com/a/22923872/8508004
    # See also https://stackoverflow.com/a/22424821/8508004 for an option to 
    # provide arguments when prototyping a full script in the notebook.
    ###-----------------for parsing command line arguments-----------------------###
    import argparse
    parser = argparse.ArgumentParser(prog='find_sequence_element_occurrences_in_sequence.py',
        description="find_sequence_element_occurrences_in_sequence.py takes \
        a short sequence represented as a string and a source FASTA-formatted \
        sequence (either a file or URL to a file) and makes an accounting of \
        the occurrences of that short sequence element (pattern) on both \
        strands of the main sequence. Exact matches only are counted as \
        occurrences. Matching is case-insensitive.    \
        **** Script by Wayne Decatur   \
        (fomightez @ github) ***")

    parser.add_argument("element", help="Sequence of element to search for. \
        For example, to search for an EcoRI site, provided `GAATTC`, without \
        any quotes or backticks, in the call to the script.\
        ", metavar="ELEMENT")
    parser.add_argument("source", help="Filepath or URL of FASTA-formatted \
        sequence to scan for occurrences of the sequence element. \
        ", metavar="SOURCE")
    parser.add_argument('-id', '--id_of_seq_element', action='store', type=str, 
        help="**OPTIONAL**Identifier \
        to reference sequence element. If none is provided, the sequence of \
        the user provided element will be used, limited to first {} bases if \
        exceeds that length.".format(limit_of_name))





    #I would also like trigger help to display if no arguments provided because 
    # need at least one for url
    if len(sys.argv)==1:    #from http://stackoverflow.com/questions/4042452/display-help-message-with-python-argparse-when-script-is-called-without-any-argu
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    element = args.element
    source= args.source
    id_of_seq_element= args.id_of_seq_element



    main()

#*******************************************************************************
###-***********************END MAIN PORTION OF SCRIPT***********************-###
#*******************************************************************************
