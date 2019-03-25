#!/usr/bin/env python
# extract_subsequence_from_FASTA.py
__author__ = "Wayne Decatur" #fomightez on GitHub
__license__ = "MIT"
__version__ = "0.1.0"


# extract_subsequence_from_FASTA.py by 
# Wayne Decatur
# ver 0.1
#
#*******************************************************************************
# Verified compatible with both Python 2.7 and Python 3.6; written initially in 
# Python 3. 
#
#
# PURPOSE:  Takes a sequence file (FASTA-format), an identifier, and start and 
# end range and extracts a sub sequence covering that region from the matching
# sequence in the provided FASTA file. Produces the sequence in FASTA format. 
# (The provided FASTA-formatted sequence file is assumed by default to be a 
# multi-FASTA, i.e., multiple sequences in the provided file, although it 
# definitely doesn't have to be. In case it is only a single sequence, the
# record id becomes moot, see below.) 
#
# Note that if there is only one record in the specified sequence file, the 
# record id is moot and you can instead provide any string for that parameter 
# as it will be ignored. This makes the script more flexible in cases where 
# sequence files aren't complex as the user doesn't need to provide an actual 
# record id.
# 
# I'll mention it here in case there are concerns from folks that know that 
# python/bopython are zero-indexed and wondering what numbering system might be
# used here. Every effort has been made to insure that numberings provided
# and returned assume typical, 'common language' mumbering formn where the first 
# residue is numbered one. In other words, the zero-indexing is kept 
# behind-the-scenes.
#
#
#
#
# If there is more than one sequence entry in the provided file than it will use
# the provided sequence identifier to parse the specified sequence, assuming 
# there is only one match. In other words it stops looking for the specified 
# sequence identifier after getting the first match. (ACtually because using
# pyfaidx, it might causes an error before that?)
#
#
# 
#
# Written to run from command line or imported into/pasted/loaded inside a 
# Jupyter notebook cell 
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
# pyfaidx
#
#
# VERSION HISTORY:
# v.0.1. basic working version.

#
# To do:
# - use chunking of 70 to save the FASTA sequence as multiline and not just
# one long line. (see gist for worked out code.) RERUN THE DEMO AND SAVE NEW VERSION.
#
#
#
#
# TO RUN:
# Examples,
# Enter on the command line of your terminal, the line
#-----------------------------------
# python extract_subsequence_from_FASTA.py SEQ_FILE 101-200
#-----------------------------------
#
# Issue `extract_subsequence_from_FASTA.py -h` for details.
# 
#
#
# To use this after importing/pasting or loading into a cell in a Jupyter 
# notebook, specify the sequence file (or results as a string) in the 
# call to the main function similar to below:
# extract_subsequence_from_FASTA("test.clustal","ScYBM",region_str="101-200")
# -OR-
# to also keep the original full description line in final output:
# extract_subsequence_from_FASTA("test.clustal","ScYBM",region_str="101-200", keep_description=True)
# -OR-
# to return the result as a string:
# result = extract_subsequence_from_FASTA("test.clustal","ScYBM",region_str="101-200", return_record_as_string=True)
#
# 
#
'''
CURRENT ACTUAL CODE FOR RUNNING/TESTING IN A NOTEBOOK WHEN IMPORTED/LOADED OR 
PASTED IN ANOTHER CELL:
result = extract_subsequence_from_FASTA("test.clustal", "ScYBM",region_str="101-200", return_record_as_string=True)
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


suffix_for_saving = "_extracted" #to be used for naming the output automatically
# when running script from command line to act on an input file

coordinates_delimiter_default = "-" #change to ":" to use a colon to specify the 
# positions range to span. Mainly meant for advanced/power users because for
# the command line you can just use the `--use_colon` (or `-uc`) flag. And if
# using Jupyter cell you can specify `use_colon = True` when calling the main 
# function.



#
#*******************************************************************************
#**********************END USER ADJUSTABLE VARIABLES****************************













#*******************************************************************************
#*******************************************************************************
###DO NOT EDIT BELOW HERE - ENTER VALUES ABOVE###

import sys
import os
#from Bio import SeqIO
#from Bio.Seq import Seq 
#from Bio.SeqRecord import SeqRecord 
#from Bio.Alphabet import generic_dna
from pyfaidx import Fasta



###---------------------------HELPER FUNCTIONS---------------------------------###




def generate_output_file_name(file_name,text_to_match, suffix_for_saving):
    '''
    Takes a file name as an argument and returns string for the name of the
    output file. The generated name is based on the text to match.


    Specific example
    =================
    Calling function with
        ("sequence.fa", "YDL140C", "_")
    returns
        "seq_YDL140C.fa"
    '''
    return "seq" + suffix_for_saving + text_to_match + ".fa"




###--------------------------END OF HELPER FUNCTIONS---------------------------###
###--------------------------END OF HELPER FUNCTIONS---------------------------###



#*******************************************************************************
###------------------------'main' function of script---------------------------##

def extract_subsequence_from_FASTA(
    sequence, id_, region_str, use_colon = False,
    keep_description = False,
    return_record_as_string = False,
    suffix_for_saving = suffix_for_saving, name_basis="alignment.clustal"):
    '''
    Main function of script. 
    It will take an sequence file in FASTA format , an identifier, and start 
    and end range and extracts a sub sequence covering that region from the 
    matching sequence in the provided FASTA file. Produces the sequence in FASTA 
    format. (The provided FASTA-formatted sequence file is assumed by default 
    to be a multi-FASTA, i.e., multiple sequences in the provided file, although 
    it definitely doesn't have to be. In case it is only a single sequence, the
    record id becomes moot, see below.) 

    The FASTA-formatted sequence file is assumed by default to be a 
    multi-FASTA, i.e., multiple sequences in the provided file, although it 
    definitely doesn't have to be. In case it is only a single sequence, the 
    record id becomes moot and users can provide anything for this parameter. 

    Has `return_record_as_string` option so it can be imported and used in
    IPython or Jupyter notebook and the result used directly without need to
    access generated file. `keep_description= True` can be used in calling 
    function to keep entire description line.
    '''
    # Bring in the sequence data:
    #---------------------------------------------------------------------------

    # get fasta records using pyfaidx
    records = Fasta(sequence)
    #for record in SeqIO.parse(sequence, "fasta"):
    #   records.append(record)
    single_record = False
    # Note on next line have to cast to list becaude pyfaidx records don't have 
    # `len`; otherwise get `TypeError: object of type 'Fasta' has no len()`.
    if len(list(records)) == 1:   
        id_ = records[0].name
        single_record = True
        # feedback
        sys.stderr.write("Single sequence with id of '{}' provided in the "
            "sequence file.\nIt will be used as the source of the sequence "
            "covering the provided positions.\n\n".format(id_))

    # feedback
    #sys.stderr.write("Sequence file read...")



    # Parse the region_str to get the start and end positions of the reference 
    # sequence to specify what corresponding segment to extract from each of 
    # the aligned sequences.
    #---------------------------------------------------------------------------
    if use_colon:
        coordinates_delimiter= ":"
    else:
        coordinates_delimiter = coordinates_delimiter_default
    region_str_parts = region_str.split(coordinates_delimiter)
    start, end = int(region_str_parts[0]), int(region_str_parts[1])
    # sanity checks
    assert start < end, (
    "The user-supplied 'start' ({}) must be less than "
    "'end' ({}).".format(start,end))




    # Collect the corresponding sequence.
    #---------------------------------------------------------------------------
    # If do not need to keep entire description line, can take advantage of
    # pyfaidx ability to use keys to get sequence without looping eand it
    # also  eturns FASTA format automatically with coordinates there!! 
    # No need to loop or edit description line to add coordinates or make 
    # FASTA. However cannot do that if want to keep the entire description line,
    # referred to as `.long_name` by pyfaidx.
    # Note: pyfaidix slices like normal Python list but it puts 
    # typical-based coordinates in description line when using keys.
    if keep_description: 
        seq_fa = "NOT_ANY_FOUND"
        for record in records:
            if id_ == record.name:
                if start-1 < 0:
                    start = 0
                seq_fa = (">" + record.long_name +":"+str(start)+
                    coordinates_delimiter_default+str(end)+ "\n"+
                    str(record[start-1:end]))+ "\n"
                # The way pyfaididx, works, start of 1 and end 50 will get, first 
                # through 50, using above line. Which means it sliced like normal
                # Python list.
                break
        if seq_fa == "NOT_ANY_FOUND":
            sys.stderr.write("**ERROR:No match to identifier found "
                "for ANY sequence record.  ***ERROR*** \nEXITING.\n")
            sys.exit(1)
    else:
        if id_ in records:
            seq_fa = records[id_][start-1:end]+ "\n"
        else:
            sys.stderr.write("**ERROR:No match to identifier found "
                "for ANY sequence record.  ***ERROR*** \nEXITING.\n")
            sys.exit(1)







    # Save the extracted FASTA record
    #---------------------------------------------------------------------------
    output_file_name = generate_output_file_name(
        sequence,id_,suffix_for_saving)
    with open(output_file_name, 'w') as output:
        if keep_description: 
            output.write(seq_fa)
        else:
            output.write(repr(seq_fa)) #this keeps FASTA format; whereas 
            # `str(seq_fa)` just outputs sequence portion as string
    # Feedback
    sys.stderr.write("\n\n*****************DONE**************************\n")
    sys.stderr.write("Extracted sequence saved in FASTA format as '{}'.".format(
            output_file_name))
    sys.stderr.write("\n*****************DONE**************************\n")


    # Return the sequence record as a string if requested
    if return_record_as_string:
        if keep_description: 
            return seq_fa
        else:
            return(repr(seq_fa))



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
    kwargs['region_str'] = region_str
    kwargs['use_colon'] = use_colon
    kwargs['keep_description'] = keep_description
    kwargs['suffix_for_saving'] = suffix_for_saving
    extract_subsequence_from_FASTA(sequence,id_,**kwargs)
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
        'extract_subsequence_from_FASTA.py',
        description="extract_subsequence_from_FASTA.py \
        takes a sequence file (FASTA-format), an identifier, and start and end \
        range and extracts a sub sequence covering that region from the \
        matching sequence in the provided FASTA file. Produces the sequence in \
        FASTA format. (The FASTA-formatted \
        sequence file is assumed by default to be \
        a multi-FASTA, i.e., multiple sequences in the provided file, \
        although it definitely doesn't have to be. In case it is only a \
        single sequence, the record id becomes moot, see below.) \
        **** Script by Wayne Decatur   \
        (fomightez @ github) ***")

    parser.add_argument("sequence_file", help="Name of sequence file to \
        use as input. Must be FASTA format. Can be a \
        multi-FASTA file, i.e., multiple sequences in FASTA format in one \
        file.", metavar="SEQUENCE_FILE")
    parser.add_argument("record_id", help="Specific identifier of sequence \
        entry in sequence file to mine. If the provided sequence file only \
        contains one sequence, that sequence will be mined and what is provided \
        for this parameter will be ignored. In other words, if the sequence \
        file is not a multi-FASTA file, you don't need to determine the \
        identifier and can instead just enter `blahblah` or any other \
        nonsensical string in this spot.", metavar="RECORD_ID")
    parser.add_argument("start_end", help="Positions to span in extracting. \
        Provide \
        start position and end position coordinates separated by the region \
        delimiter which is '{}' by default. You can use `--use_colon` flag to \
        change to a colon, for using something like, `201:405` instead of \
        `201{}405`. (Coordinates are meant to refer to 'common' numbering \
        scheme where first residue is numbered one, etc.)\
        ".format(coordinates_delimiter_default,coordinates_delimiter_default), 
        metavar="START_and_END")

    parser.add_argument("-uc", "--use_colon",help=
    "Add this flag to be able to specify that you want to use a colon in for \
    specifying the positions to extact from the corresponding sequences.",
    action="store_true")

    parser.add_argument("-kd", "--keep_description",help=
    "Add this flag to keep entire description line in output.",
    action="store_true")

    parser.add_argument('-os', '--output_suffix', action='store', type=str, 
    default= suffix_for_saving, help="OPTIONAL: Set a suffix for including in \
    file name of output. \
    If none provided, '{}' will be used.".format(suffix_for_saving))



    #I would also like trigger help to display if no arguments provided because 
    # need at least one for url
    if len(sys.argv)==1:    #from http://stackoverflow.com/questions/4042452/display-help-message-with-python-argparse-when-script-is-called-without-any-argu
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    sequence = args.sequence_file
    id_ = args.record_id
    region_str = args.start_end
    use_colon = args.use_colon
    keep_description = args.keep_description
    suffix_for_saving = args.output_suffix


    main()

#*******************************************************************************
###-***********************END MAIN PORTION OF SCRIPT***********************-###
#*******************************************************************************
