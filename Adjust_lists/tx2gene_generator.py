#!/usr/bin/env python
# tx2gene_generator.py by Wayne Decatur
__author__ = "Wayne Decatur" #fomightez on GitHub
__license__ = "MIT"
__version__ = "0.1.0"
#
#*******************************************************************************
# Compatible with both Python 2.7 and Python 3.6 (verified); written initially 
# in Python 2.7 to be compatible with Python 3.
#
# PURPOSE: Takes a `quant.sf` file from Salmon output for a situation where
# there is a 1:1 transcript: gene relationship & makes the `tx2gene.csv`
# file needed for using tximport to bring Salmon data into DESeq2. By default,
# the files read and saved will be those, respectively, but they can be 
# specified as arguments when executing the script.
# See the links below for information about tximport:
# https://www.bioconductor.org/help/workflows/rnaseqGene/
# and 
# http://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html
#
# In order to keep things simple with the argument parsing system, one has to 
# specify both input and output file if looking to not use the standard 
# default for output. 
#
# Even though it is stated as being made to work with Salmon's `quant.sf`-type
# output, it will actually work with any gene/transcript list where there is a
# single header/column names line first and the gene/transcript designation is 
# the first word or column (tab-separated) for pretty much anything,
# but comma-separated input. It could be easily edited to work to take a comma-
# separated list by editing the `split` command. The caveat here is this is made
# to work where there is a 1:1 relationship with transcripts to genes as there
# is generally assumed in practice with early branching eukaryotes, like for
# baker's/budding yeast.
#
# (Note to self: See
# `Downstream analysis of April 2017 data with tximport and DESeq2.md` and 
# `Checking for enrichment of NME1 in mito RNA-seq results quantified with Salmon.md`
# for info about this script).
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
#
# To do yet:
# - ?
#
#
#
# TO RUN:
# Example,
# Enter on the command line of your terminal, the line
#-----------------------------------
# python tx2gene_generator.py [INPUT_FILE][OUTPUT_FILE]
#-----------------------------------
# Specifying the input and output files is optional, run 
# `python tx2gene_generator.py -h` for details.
#
#
#*******************************************************************************
#


#*******************************************************************************
##################################
#  USER ADJUSTABLE VALUES        #

##################################
#
input_file_name_default = "quant.sf"
output_file_name_default = "tx2gene.csv"






#
#*******************************************************************************
#**********************END USER ADJUSTABLE VARIABLES****************************


















#*******************************************************************************
#*******************************************************************************
###DO NOT EDIT BELOW HERE - ENTER VALUES ABOVE###

import sys
import os



###---------------------------HELPER FUNCTIONS---------------------------------###


# NONE USED


###--------------------------END OF HELPER FUNCTIONS---------------------------###
###--------------------------END OF HELPER FUNCTIONS---------------------------###











#*******************************************************************************
###-----------------Actual Main portion of script---------------------------###


def main():
    """ Main entry point of the script """
    # prepare output file for saving so it will be open and ready
    with open(output_file_name, 'w') as output_file:

        # add the header line to output file
        output_file.write('transcript,gene')
        

        # read in the text file but skip header (see http://stackoverflow.com/questions/4796764/read-file-from-line-2-or-skip-header-row)
        with open(input_file_name, 'r') as input_handler:
            next(input_handler)


            # prepare to give feeback later
            lines_processed = 0

            for line in input_handler:
                lines_processed += 1
                info = line.split()
                # print info # ONLY FOR DEBUGGING
                transcript = info[0]


                # Send text to output
                output_file.write("\n"+transcript+","+ transcript)



        # Completed scan of input file and therefore close file and give feedback.
        input_handler.close() #not actually necessary due to use of `with open`
        sys.stderr.write( "\n"+ str(lines_processed) + " transcripts listed in '" + input_file_name + "' have been processed to make a table of transcript to gene conversions.")

        



    # Completed generating a file equivalnet to `tx2gene.csv` and so provide feedback.
    sys.stderr.write( "\nOutput to be used with tximport saved as '"+ output_file_name  + "'.\n")







if __name__ == "__main__" and '__file__' in globals():
    """ This is executed when run from the command line """
    # Code with just `if __name__ == "__main__":` alone will be run if pasted
    # into a notebook. The addition of ` and '__file__' in globals()` is based
    # on https://stackoverflow.com/a/22923872/8508004
    # See also https://stackoverflow.com/a/22424821/8508004 for an option to 
    # provide arguments when prototyping a full script in the notebook.
    ###-----------------for parsing command line arguments-----------------------###
    import argparse
    parser = argparse.ArgumentParser(prog='tx2gene_generator.py',
        description="tx2gene_generator.py  takes a `quant.sf` file from Salmon \
        output and makes the `tx2gene.csv` file needed for using tximport to bring \
        Salmon data into DESeq2.   By default, the files read and saved will be \
        those, respectively, but they can be specified as arguments when executing \
        the script.      \
        **** Script by Wayne Decatur   \
        (fomightez @ github) ***")

    parser.add_argument("input", nargs='?', help="**OPTIONAL**Name of the file \
        generated by Salmon \
        when run with your transcriptome of interest. Usually, this is \
        '"+input_file_name_default+"' &\
        if no input file name is provided then this will be used by \
        default.", default=input_file_name_default, metavar="INPUT_FILE")
    parser.add_argument("output", nargs='?', help="**OPTIONAL**Name of file to \
        save results. If BOTH input \
        and output file are not provided, '"+output_file_name_default+"', will \
        be used.", default=output_file_name_default, metavar="OUTPUT_FILE")
    # Note see https://stackoverflow.com/questions/18862836/how-to-open-file-using-argparse#comment35222484_18863004
    # for why not using `argparse.FileType` approach here.
    # See
    # https://stackoverflow.com/questions/4480075/argparse-optional-positional-arguments 
    # and 
    # https://docs.python.org/2/library/argparse.html#nargs for use of `nargs='?'` 
    # to make input and output file names optional. Note that the square brackets
    # shown in the usage out signify optional according to 
    # https://stackoverflow.com/questions/4480075/argparse-optional-positional-arguments#comment40460395_4480202
    # , but because placed under positional I added clarifying text to help 
    # description.
    # IF MODIFYING THIS SCRIPT FOR USE ELSEWHERE AND DON'T NEED/WANT THE INPUT AND
    # OUTPUT FILES TO BE OPTIONAL, remove `nargs` (& default?) BUT KEEP WHERE NOT
    # USING `argparse.FileType` AND USING `with open` AS CONISDERED MORE PYTHONIC.

    #Normally I like to trigger help to display if no arguments provided, but 
    # because both input and output file names optional, I won't do that here.


    args = parser.parse_args()
    # Note see https://stackoverflow.com/questions/18862836/how-to-open-file-using-argparse#comment35222484_18863004
    # for why not using `argparse.FileType` approach here.
    input_file_name = args.input
    output_file_name = args.output


    main()

#*******************************************************************************
###-***********************END MAIN PORTION OF SCRIPT***********************-###
#*******************************************************************************
