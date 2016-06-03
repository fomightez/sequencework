#!/usr/bin/env python
# makes_length_annotation_file.py by Wayne Decatur
# ver 0.1
#
#*******************************************************************************
# USES Python 2.7 but should be convertable via 2to3, see https://docs.python.org/3.0/library/2to3.html
#
# PURPOSE: Takes a file of `biomart_length.txt` and returns a file with gene_id
# and then gene length for use with NOISeq. See the User's Guide for NOISeq
# at http://www.bioconductor.org/packages/release/bioc/html/NOISeq.html .
# The provided file has the higher number always second, but I made the script
# so that it will work no matter the order the end or start is provided.
#
#
# Users will need to change the name of the file in the script to suit their
# needs.
# A file of the output will be saved in the same directory in which the script
# is run.
#
#
# (Note to self: see `making length annotation file from biomart info.md` for info about this script).
#
#
#
# Dependencies beyond the mostly standard libraries/modules:
#
#
#
# v.0.1. basic working version
#
#
# To do:
# -
#
# TO RUN:
# Example,
# Enter on the command line of your terminal, the line
#-----------------------------------
# makes_length_annotation_file.py
#-----------------------------------
#
#
#*******************************************************************************
#


#*******************************************************************************
##################################
#  USER ADJUSTABLE VALUES        #

##################################
#
biomart_input_file = "biomart_length.txt"






#
#*******************************************************************************
#**********************END USER ADJUSTABLE VARIABLES****************************


















#*******************************************************************************
#*******************************************************************************
###DO NOT EDIT BELOW HERE - ENTER VALUES ABOVE###

import sys
import os



###---------------------------HELPER FUNCTIONS---------------------------------###


def generate_output_file_name(file_name):
    '''
    Takes a file name as an argument and returns string for the name of the
    output file. The generated name is based on the original file
    name.

    Specific example
    ================
    Calling function with
        ("biomart_length.txt")
    returns
        "biomart_length_converted_for_NOISeq.txt"
    '''
    main_part_of_name, file_extension = os.path.splitext(
        file_name) #from http://stackoverflow.com/questions/541390/extracting-extension-from-filename-in-python
    if '.' in file_name:  #I don't know if this is needed with the os.path.splitext method but I had it before so left it
        return main_part_of_name + "_converted_for_NOISeq" + file_extension
    else:
        return file_name + "_converted_for_NOISeq"

def generate_output_file(provided_text):
    '''
    function text and saves it as a text file
    '''
    name_of_file_to_save = generate_output_file_name(biomart_input_file)
    data_file_stream = open(name_of_file_to_save , "w")
    data_file_stream.write(provided_text.rstrip('\r\n')) #rstrip to remove trailing newline
    # from http://stackoverflow.com/questions/275018/how-can-i-remove-chomp-a-newline-in-python
    data_file_stream.close()
    sys.stderr.write( "\nConversion sucessfully completed! Text saved as '{0}'.\n".format(name_of_file_to_save))


###--------------------------END OF HELPER FUNCTIONS---------------------------###
###--------------------------END OF HELPER FUNCTIONS---------------------------###










#*******************************************************************************
###-----------------Actual Main portion of script---------------------------###



# read in the text file but skip header (see http://stackoverflow.com/questions/4796764/read-file-from-line-2-or-skip-header-row)
with open(biomart_input_file , 'r') as input_handler:
    next(input_handler)

    # prepare to store information to be saved later
    text_to_save = ""

    # prepare to give feeback later
    lines_processed = 0

    for line in input_handler:
        lines_processed += 1
        info = line.split()
        # print info # ONLY FOR DEBUGGING
        gene_id = info[0]
        gene_start = int(info[1]) #cast to int to make calculation easy
        gene_end = int(info[2])
        gene_length = (abs(gene_end - gene_start)) + 1 # with addition of absolute value,
        # calculation future-proofed in case for Crick strand the actual end of
        # transcript (which is acutally a lower number) provided second. Input
        # script designed for originally didn't do that but it could be provided
        # that way in other places if TRUE start and TRUE end sites of transcript
        # considered.
        text_to_save += gene_id + "\t" + str(gene_length) + "\n"



# Completed scan of input file and therefore close file and give feedback.
input_handler.close()
sys.stderr.write( "\n"+ str(lines_processed) + " genes listed in '" + biomart_input_file + "' and processed to produce gene lengths.")



# Save results
generate_output_file(text_to_save)




#*******************************************************************************
###-***********************END MAIN PORTION OF SCRIPT***********************-###
#*******************************************************************************
