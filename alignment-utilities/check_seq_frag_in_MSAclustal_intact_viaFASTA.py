#!/usr/bin/env python
# check_seq_frag_in_MSAclustal_intact_viaFASTA.py
__author__ = "Wayne Decatur" #fomightez on GitHub
__license__ = "MIT"
__version__ = "0.1.0"


# check_seq_frag_in_MSAclustal_intact_viaFASTA.py by 
# Wayne Decatur
# ver 0.2
#
#*******************************************************************************
# Verified compatible with both Python 2.7 and Python 3.6; written initially in 
# Python 3. 
#
#
# PURPOSE: Takes a multiple sequence alignment in Clustal format and checks if
# for a specified identifier that the sequence there is contiguous by comparing
# with a provided FASTA sequence file. Doesn't need it to span the entire 
# length, just at least be an internal match.
# Meant to be used to validate no erroneoous deletions were introduced in the 
# course of manually editing a multiple sequence alignment that doesn't cover
# the entire span of a sequence. For example, if the MSA is meant to align
# a domain made up by a contiguous portion of the chain.
# All gaps are disregarded; the point is to check the sequence of residues.
# In contrast to the related script 
# `check_seq_in_MSAclustal_consistent_with_FASTA.py` which
# checked the entire sequence span in the MSA matched the entire sequence in the
# FASTA file, this script just check the represented residues in the sequence
# all match contiguously with a single region in the provided FASTA sequence 
# file. In other words, what is at the very start or very end of the sequence
# can be left off in MSA. 
# Optionally, the first residue number of the fragment in multiple sequence 
# alignment can also be returned. This value will be in common terms where the
# first position in a sequence would be number one, etc..
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
# If you are Wayne, see 
# `alignment of yeast XXX1 vs XXX2 vs XXX12 vs archael XXXX from Swiss Model.clw` 
# for specific impetus behind this script.
#
#
#
# Dependencies beyond the mostly standard libraries/modules:
# biopython
#
#
# VERSION HISTORY:
# v.0.1. basic working version
# v.0.2. adds returning first position of sequence within the MSA

#
# To do:
# - 
#
#
#
#
# TO RUN:
# Examples,
# Enter on the command line of your terminal, the line
#-----------------------------------
# python check_seq_frag_in_MSAclustal_intact_viaFASTA.py ALIGNMENT_TEXT_FILE seq_id seq.fa
#-----------------------------------
#
# Issue `check_seq_frag_in_MSAclustal_intact_viaFASTA.py -h` for details.
# 
#
#
# To use this after importing/pasting or loading into a cell in a Jupyter 
# notebook, specify at least the alignment file, id of sequence to analyze, and 
# ungapped FASTA file of the sequence in the call to the main function similar 
# to below:
# consistent = check_seq_frag_in_MSAclustal_intact_viaFASTA("alignment.clw", "VPH1", "S288C_YOR270C_VPH1_protein.fsa", return_TF = True)
# -OR-
# To also get the position numbering of the first residue of the sequence
# represented in the multiple sequence alignment, call function with 
# `report_start_residue_pos` set to true as well, like so:
# consistent, start_pos = check_seq_frag_in_MSAclustal_intact_viaFASTA("alignment.clw", "VPH1", "S288C_YOR270C_VPH1_protein.fsa", return_TF = True, report_start_residue_pos= True)
#
# 
#
'''
CURRENT ACTUAL CODE FOR RUNNING/TESTING IN A NOTEBOOK WHEN IMPORTED/LOADED OR 
PASTED IN ANOTHER CELL:
consistent = check_seq_frag_in_MSAclustal_intact_viaFASTA("alignment.clw", "VPH1", "S288C_YOR270C_VPH1_protein.fsa", return_TF = True)
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

gap_indicator = "-" #change to "." if using dots. Used by 
# `.ungap()`.



#
#*******************************************************************************
#**********************END USER ADJUSTABLE VARIABLES****************************













#*******************************************************************************
#*******************************************************************************
###DO NOT EDIT BELOW HERE - ENTER VALUES ABOVE###

import sys
import os
from collections import Counter
from Bio import SeqIO
from Bio import AlignIO





###---------------------------HELPER FUNCTIONS-------------------------------###








###--------------------------END OF HELPER FUNCTIONS---------------------------###
###--------------------------END OF HELPER FUNCTIONS---------------------------###



#*******************************************************************************
###------------------------'main' function of script--------------------------##

def check_seq_frag_in_MSAclustal_intact_viaFASTA(
    alignment, id_, fasta_fn, return_TF = False, 
    report_start_residue_pos = False):
    '''
    Main function of script. 
    It will take an alignment text file (CLUSTAL form) and check if for a \
    specified identifier that the sequence exists contiguously in an indicated 
    FASTA sequence file. It doesn't need to span end to end, just at least
    be an internal match.
    Returns whether the match is 'True' or 'False' to standard out if 
    `return_TF` is false, which is default setting. If `return_TF` is true, 
    whether it matches is passed back as as a boolean; this is intended to 
    handle when called from a function. 
    Optionally, also reports the position number of the start of the sequence
    fragment in the multiple sequence alignment. If the sequence fragment is not
    internally contiguous, it will report 'NA'. The setting of `return_TF` 
    determines in what manner the start position is reported, i.e., if it 
    printed to stdout or passed back when main function called. If 
    `report_start_residue_pos` is 
    false, which is the default setting, and `report_start_residue_pos` is true, 
    it will send the value to standard out.  If `report_start_residue_pos` is 
    true and `report_start_residue_pos` is true, the value will be passed back 
    as an integer. This is intended to handle when called from a function.
    '''
    # feedback about important option
    sys.stderr.write("\n**NOTE: gap indicator in this script is currently set "
        "to '{}'. If\nthat does not match what provided alignment uses to "
        "indicate gaps,\nplease change the setting within the script code "
        "under\n'USER ADJUSTABLE VALUES' around line 100 "
        "(give or take a few).**\n".format(gap_indicator))


    # Bring in the necessary alignment data:
    #---------------------------------------------------------------------------

    try:
        file_name = alignment
        alignment = AlignIO.read(alignment, "clustal") 
    except (TypeError,OSError,IOError) as e:
        file_name = name_basis
        try:
            from StringIO import StringIO
        except ImportError:
            from io import StringIO
        alignment = AlignIO.read(StringIO(alignment), "clustal") 
        # Note "FileNotFoundError is a subclass of OSError"
        # (https://stackoverflow.com/a/28633573/8508004)

    # feedback
    sys.stderr.write("Alignment file read...")



    # Get sequence user specified from alignment and ungap
    #---------------------------------------------------------------------------
    msa_ids = []
    for record in alignment:
        msa_ids.append(record.id)
        if record.id == id_:
            seq = record.seq
            sys.stderr.write("{} sequence collected from "
                "alignment...".format(id_))
    try: 
        if seq != None:
            pass
    except UnboundLocalError:
        sys.stderr.write("\n***    ERROR    ***   ERROR   ***\nCAUSE OF ERROR:")
        sys.stderr.write("\nNo matches to provided identifier were found in "
            "the multiple sequence alignment.\nUser specified id: {}\n"
            "Identifiers in the provided alignment: "
            "{}.\n".format(id_, repr(msa_ids)))
        sys.stderr.write("***  END OF ERROR DESCRIPTION ***\n")
        sys.exit(1)
    # Ungap the identified sequence
    seq = seq.ungap(gap_indicator)


    

    # Read in FASTA sequence and remove gaps if it contains any
    #---------------------------------------------------------------------------
    # in case mistakenly a multi-fasta file is provided, parst but try to grab
    # first one
    fasta_records = []  #intialize a list of records for the present FASTA file
    for fasta_record in SeqIO.parse(fasta_fn, "fasta"):
        fasta_records.append(fasta_record)
    fasta_seq = fasta_records[0].seq.ungap(gap_indicator)
    #Often '*' is included at the end of protein sequences to indicate stop and
    # we don't want that.
    if fasta_seq[-1] == '*':
        fasta_seq = fasta_seq[:-1]
    # feedback
    sys.stderr.write("FASTA file read...")


    # See if the sequence is in the FASTA sequence and return TF if called 
    # with `return_TF = True`. Otherwise, consider called from command line 
    # & return True or False to stdout
    # Also return the starting position depending on `report_start_residue_pos` 
    # setting.
    #---------------------------------------------------------------------------
    sys.stderr.write("Checking...\nIs sequence contiguous without regard to the"
        " ends?")
    if report_start_residue_pos:
        sys.stderr.write("Will return start position value, too.")
    sys.stderr.write("...  ...\n")
    if seq in fasta_seq:
        seq_in_fasta = True
        start_residue_pos = str(fasta_seq).index(str(seq)) + 1 # `+1` to convert 
        # to more common convention of first residue being number ONE. In other 
        # words, keep python/bipython zero indexing behind-the-scenes.
    else:
        seq_in_fasta = False
        start_residue_pos = 'NA'
    if return_TF and report_start_residue_pos:
        return seq_in_fasta, start_residue_pos
    elif return_TF and not report_start_residue_pos:
        return seq_in_fasta
    elif report_start_residue_pos:
        import time
        time.sleep(0.05) #putting delay here to insure stdout last thing printed
        print (str(seq_in_fasta), str(start_residue_pos))
    else:
        import time
        time.sleep(0.05) #putting delay here to insure stdout last thing printed
        print (str(seq_in_fasta))



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
    kwargs['report_start_residue_pos'] = report_start_residue_pos
    #kwargs['descr_source'] = descr_source
    #kwargs['suffix_for_saving'] = suffix_for_saving
    check_seq_frag_in_MSAclustal_intact_viaFASTA(alignment, id_, fasta_fn, **kwargs)
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
    ###-----------------for parsing command line arguments-----------------------###
    import argparse
    parser = argparse.ArgumentParser(prog=
        'check_seq_frag_in_MSAclustal_intact_viaFASTA.py',
        description="check_seq_frag_in_MSAclustal_intact_viaFASTA.py \
        takes a multiple sequence alignment in Clustal format and checks if \
        for a specified identifier that the sequence there is contiguous by \
        comparing with a provided FASTA sequence file. Doesn't need it to span \
        the entire length, just at least be an internal match. Meant to be \
        used to validate no \
        erroneous deletions were introduced in the course of manually editing \
        a multiple sequence \
        alignment that doesn't cover the entire span of a sequence. All gaps \
        are ignored.\
        **** Script by Wayne Decatur   \
        (fomightez @ github) ***")

    parser.add_argument("align_file", help="Name of file of alignment text \
        file (CLUSTAL format).\
        ", metavar="ALIGNMENT_FILE")

    parser.add_argument("id", help="Identifier that matches the sequence \
        in the MSA to check.\
        ", metavar="ID")

    parser.add_argument("fasta_file", help="Name of sequence file to check \
        against (FASTA format).\
        ", metavar="FASTA_FILE")

    parser.add_argument('-rsp', '--report_start_pos', help="Add this flag when \
        calling the script in \
        order to be able to also report the start position of the contiguous \
        fragment represented in the MSA. When the fragment isn't contiguous, \
        using this setting will result in 'NA'.\
        ", action="store_true")


    '''parser.add_argument('-rsp', '--report_start_pos', action='store', type=str, 
    default= None, help="OPTIONAL: Provide FASTA file with same \
    ids to use as source of alignment to use as source of descriptions. The \
    typical file is output as `*_extracted_ungapped.fa` files from the \
    `extract_regions_from_clustal_alignment.py` script\
    .")
    '''





    #I would also like trigger help to display if no arguments provided because 
    # need at least one for url
    if len(sys.argv)==1:    #from http://stackoverflow.com/questions/4042452/display-help-message-with-python-argparse-when-script-is-called-without-any-argu
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    alignment = args.align_file
    id_ = args.id
    fasta_fn = args.fasta_file
    report_start_residue_pos = args.report_start_pos


    main()

#*******************************************************************************
###-***********************END MAIN PORTION OF SCRIPT***********************-###
#*******************************************************************************
