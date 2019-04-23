#!/usr/bin/env python
# measure_intergenic_regions_in_mito_annotations.py
__author__ = "Wayne Decatur" #fomightez on GitHub
__license__ = "MIT"
__version__ = "0.1.0"


# measure_intergenic_regions_in_mito_annotations.py by Wayne Decatur
# ver 0.1
#
#*******************************************************************************
# Verified compatible with both Python 2.7 and Python 3.6; written initially in 
# Python 3. 
#
#
# PURPOSE: Takes takes a file of annotations (in gff form) of a mitochondrial 
# genome produced by MFannot and followed by conversion to gff via a script, the
# corresponding mitochondrial sequence in FASTA format, and 
# determines the intergenic spacing between each annotated gene (including rRNA,
# which oddly is labeled `rRNA` and not `gene` as well as `rps3` (ribosomal 
# small subunit protein 3 ) oddly labeled as `rRNA` and not a gene.). Assumes 
# the mitochondrial genome is circular.
# Note: Going from the MFAnnot master files (`.new` extesion) preferably use 
# `mfannot2gff3.pl` from
# https://github.com/yjx1217/LRSDAY/blob/master/scripts/mfannot2gff3.pl . I have 
# also used in the past `mfannot2gff.pl`from 
# https://github.com/kbseah/mitonotate/blob/master/mfannot2gff.pl` but noted it 
# suffered from an apparent one-off error for the RNAs at the start, see 
# https://github.com/kbseah/mitonotate/issues/5. And so I tried the LRSday one 
# (https://github.com/yjx1217/LRSDAY/blob/master/scripts/mfannot2gff3.pl) most
# recently. NOTE that there is a difference in output from these two. 
# mfannot2gff.pl includes notes on the translation table used in the mfannot run 
# on each feature.
# 
#
#
#
#
# Developed as a follow-up to 
# `Relating matches of cerevisiae mito promoters in yueomyces_sinensis to the 
# annotated genes and patterns in cerevisiae.md`.
#
#
# Dependencies beyond the mostly standard libraries/modules:
# sh
# ** NEEDS BLAST+ ALREADY INSTALLED IN ENVIRONMENT.
#
#
#
# VERSION HISTORY:
# v.0.1. basic working version

#
# To do:
# - verify works with Python 2 (can use https://github.com/fomightez/mcscan-blast-binder or mcscan-binder)
#
#
#
#
# TO RUN:
# Examples,
# Enter on the command line of your terminal, the line
#-----------------------------------
# python measure_intergenic_regions_in_mito_annotations.py mito_annotations.gff3 mito_seq.fa
#-----------------------------------
# Issue `measure_intergenic_regions_in_mito_annotations.py -h` for 
# details.
# 
# More examples from running from the command line are at the links below: 
# https://git.io/????????
#
#
#
# To use this after importing/pasting or loading into a cell in a Jupyter 
# notebook, specify at least the file of annotations:
# from measure_intergenic_regions_in_mito_annotations import measure_intergenic_regions_in_mito_annotations
# measure_intergenic_regions_in_mito_annotations("mito_annotations.gff3", "mito_seq.fa")
#
# 
#
'''
CURRENT ACTUAL CODE FOR RUNNING/TESTING IN A NOTEBOOK WHEN IMPORTED/LOADED OR 
PASTED IN ANOTHER CELL:
from measure_intergenic_regions_in_mito_annotations import measure_intergenic_regions_in_mito_annotations
df = measure_intergenic_regions_in_mito_annotations("mito_annotations.gff3", "mito_seq.fa", return_df = True)
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

text_to_add_to_name_of_intergenic_lengths = "_intergenic_gap_sizes"

#
#*******************************************************************************
#**********************END USER ADJUSTABLE VARIABLES****************************















#*******************************************************************************
#*******************************************************************************
###DO NOT EDIT BELOW HERE - ENTER VALUES ABOVE###

import sys
import os
import subprocess
import pandas as pd




###---------------------------HELPER FUNCTIONS---------------------------------###

def generate_output_file_name(file_name, text_to_add_to_name):
    '''
    Takes a file name as an argument and returns string for the name of the
    output file. The generated name is tagged to indicate difference.


    Specific example
    =================
    Calling function with
        ("mito_annotation.gff3", "_rnlFIXED")
    returns
        "mito_annotation_rnlFIXED.gff3"
    '''
    main_part_of_name, file_extension = os.path.splitext(
        file_name) #from http://stackoverflow.com/questions/541390/extracting-extension-from-filename-in-python
    if '.' in file_name:  #I don't know if this is needed with the os.path.splitext method but I had it before so left it
        return main_part_of_name + text_to_add_to_name  + ".tsv"
    else:
        return file_name + text_to_add_to_name + ".tsv"



###--------------------------END OF HELPER FUNCTIONS---------------------------###
###--------------------------END OF HELPER FUNCTIONS---------------------------###

#*******************************************************************************
###------------------------'main' function of script---------------------------##

def measure_intergenic_regions_in_mito_annotations(gff_file_name, 
        seq_file_name, return_df=False):
    '''
    Main function of script. 
    Takes takes a file of annotations (in gff form) of a mitochondrial 
    genome produced by MFAnnot and followed by conversion to gff via a script,
    the corresponding fungal mitochondrial genome (FASTA format), 
    and determines spacing between genes.

    Has options to return a dataframe that can be used further. This is
    meant for when calling the main function via IPython/Jupyter.

    Note: Going from the MFAnnot master files (`.new` extesion) preferably use 
    `mfannot2gff3.pl` from 
    https://github.com/yjx1217/LRSDAY/blob/master/scripts/mfannot2gff3.pl . I 
    have also used in the past `mfannot2gff.pl`from 
    https://github.com/kbseah/mitonotate/blob/master/mfannot2gff.pl` but noted 
    it suffered from an apparent one-off error for the RNAs at the start, see 
    https://github.com/kbseah/mitonotate/issues/5. And so I tried the LRSday one 
    (https://github.com/yjx1217/LRSDAY/blob/master/scripts/mfannot2gff3.pl) most
    recently.
    '''
    # Check sequence file present and then mine length for use later:
    #---------------------------------------------------------------------------
    assert os.path.isfile(seq_file_name), ("Specified sequence file '{}' "
            "not found.".format(seq_file_name))
    from pyfaidx import Fasta
    #assume genome file used to make annotations supplied and so first sequence
    # (and only) is genome
    records = Fasta(seq_file_name)
    # Note on next line have to cast to list becaude pyfaidx records don't have 
    # `len`; otherwise get `TypeError: object of type 'Fasta' has no len()`.
    single_record = True
    if len(list(records)) != 1:   
        single_record = False
        # feedback
        sys.stderr.write("Expected single genome file. '{}' contains {} "
            "sequences. Proceeding but only using first.\n\n".format(
            seq_file_name, len(list(records))))
    seq = records[0]
    seq_length = len(str(seq))
    sys.stderr.write("Provided genome '{}' is {} bps in length.".format(
        seq_file_name, seq_length))



    # Check annotation file present and then mine starts and ends of genes:
    #---------------------------------------------------------------------------
    assert os.path.isfile(gff_file_name), ("Specified annotation file '{}' "
            "not found.".format(gff_file_name))
    start_list = []
    end_list = []
    genes = []
    with open(gff_file_name, 'r') as gff_handler:
            for line in gff_handler:
                if not line.startswith("#") and len(line.split("\t")) > 1:
                #if (len(line.split("\t")) > 1 and line.split("\t")[2] == "gene") or (len(line.split("\t")) > 1 and line.split("\t")[2] == "rRNA"):
                #if len(line.split("\t")) > 1 and (line.split("\t")[2] == "gene" or line.split("\t")[2] == "rRNA" ):
                    if (line.split("\t")[2] == "gene" or 
                        line.split("\t")[2] == "rRNA"):
                        start_list.append(int(line.split("\t")[3]))
                        end_list.append(int(line.split("\t")[4]))
                        genes.append(line.split("ID=")[1].split(";")[0])
    # make dataframe with colleced information
    df = pd.DataFrame(list(zip(genes, start_list, end_list)),columns=[
        'gene','start', 'end'])

    #make a list of the gaps ACCOUNTING FOR CIRCULAR NATURE
    #---------------------------------------------------------------------------
    intergenic_gaps = []
    for row in df.itertuples():
        # Determine gap from last to first ACCOUNTING FOR CIRCULAR NATURE
        # on first row simply collect distance from start to first feature 
        # assuming first postion is numbered ONE. <-- I probably should check 
        # thatis true the way MFannot (and converting script) works.
        # On last row then take the distance from the last feature to the end of
        # the chromosome sequence and then add that to the 
        # distance_to_first_feature to get the intergenic distance that spans
        # the circularization point.
        if row.Index == 0:
            distance_to_first_feature = row.start
        elif row.Index == (len(df)-1):
            # calculate distance from last feature to this one first.
            gap = row.start - df.iloc[row.Index-1].end
            intergenic_gaps.append(gap)
            # now calculate the distance from the last feature around to the
            # first one because of circular nature of mitochondria
            last_gap = (seq_length - row.end) + distance_to_first_feature
            intergenic_gaps.append(last_gap)
        else:
            gap = row.start - df.iloc[row.Index-1].end
            intergenic_gaps.append(gap)
    from statistics import mean
    from statistics import median
    gap_mean_size = mean(intergenic_gaps)
    gap_median_size = median(intergenic_gaps)



    # Reporting and Saving continued
    #---------------------------------------------------------------------------

    #save a file of the gaps
    
    sys.stderr.write("\nThe intergenic gaps observed:\n{}\n\n".format(
        ",".join(str(x) for x in intergenic_gaps)))
    sys.stderr.write("The mean size of the gaps is {} bps.\n".format(
        gap_mean_size))
    sys.stderr.write("The median size of the gaps is {} bps.\n\n".format(
        gap_median_size))
    gap_df = pd.DataFrame(intergenic_gaps, columns=['gaps (bps)']) 
    output_fn = generate_output_file_name(gff_file_name,
        text_to_add_to_name_of_intergenic_lengths)
    gap_df.to_csv(output_fn, sep='\t',index = False) #add `,header=False` ?
    sys.stderr.write("A file listing the intergenic gaps has been "
        "saved as '{}'.\n".format(output_fn))
    if return_df:
        return gap_df
        sys.stderr.write("Additionally, a dataframe of the intergenic gaps "
            "was returned.\nThe mean gap size should be given by `df.mean(0)`"
            ".\n The median should be given by `df.median(0)`.\n")


        



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
    #kwargs['return_df'] = False #don't want dataframe returned if 
    # calling script from command line; don't need to define here because set as 
    # default in function now
    measure_intergenic_regions_in_mito_annotations(gff_file_name, 
        seq_file_name,**kwargs)
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
    parser = argparse.ArgumentParser(prog='measure_intergenic_regions_in_mito_annotations.py',
        description="measure_intergenic_regions_in_mito_annotations.py \
        takes output from MFannot that has been converted to a gff file by a\
        subsequent script, a corresponding fungal mitochondrial sequence file, \
        and determines intergenic gap sizes. Assumes CIRCULAR genome.\
        **** Script by Wayne Decatur   \
        (fomightez @ github) ***")

    parser.add_argument("gff_file", help="Name of annotaion results \
        file (gff format) from a mitochondrial genome to parse and possibly \
        fix.\
        ", metavar="GFF_FILE")
    parser.add_argument("sequence_file", help="Name of file containing the \
        fungal mitochondrial sequence (FASTA format) corresponding to the \
        annotation file.\
        ", metavar="SEQ_FILE")




    #I would also like trigger help to display if no arguments provided because 
    # need at least one for url
    if len(sys.argv)==1:    #from http://stackoverflow.com/questions/4042452/display-help-message-with-python-argparse-when-script-is-called-without-any-argu
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    gff_file_name = args.gff_file
    seq_file_name = args.sequence_file


    main()

#*******************************************************************************
###-***********************END MAIN PORTION OF SCRIPT***********************-###
#*******************************************************************************
