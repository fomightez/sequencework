#!/usr/bin/env python


# cufflinksIDs_list_to_systematic_names.py by Wayne Decatur



#*******************************************************************************
# USES Python 2.7 but should be convertable via 2to3, see https://docs.python.org/3.0/library/2to3.html
#
# PURPOSE: Using a `gtf` file that is output from Cufflinks, converts a list of
# `XLOC_` gene ids in a file to systematic yeast gene ids.
# The list should be gene ids each on a separate line of the file.
# It was written for yeast genes but may work with a gtf annotation file for any
# organism if the tags involved in the annotation file include `oId` and
# `gene_id` and , possibly, `nearest_ref`.
#
#
#
#
# see `Trying to convert Shafi's cufflinks list to something that works with other lists with systematic names.md`
# for info about this script
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
#
# TO RUN:
# Example,
# Enter on the command line of your terminal, the line
#-----------------------------------
# python cufflinksIDs_list_to_systematic_names.py list.txt
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
gtf_annotation_file_name = "merged.gtf"


cleared_as_being_converted_properly = ["HRA1", "PWR1", "Q0045", "tA(AGC)K2", "snR51" ] #there were some issues with
# some genes having 'CUFF.' as text being used for conversion and
# I used `nearest_ref` for those, but several didn't look like systematic yeast
# names even with using those. I checked and these ones confrimed converted fine
# and they won't have systematic names that match most since are non-coding RNAs
#  or on mitochondrial genome. YOU MAY NEED TO ADD MORE AS NEEDED IF GTF FILE
# UPDATED OR REPLACED WITH ONE FROM A DIFFERENT SOURCE.


#
#*******************************************************************************
#**********************END USER ADJUSTABLE VARIABLES****************************





















#*******************************************************************************
#*******************************************************************************
###DO NOT EDIT BELOW HERE - ENTER VALUES ABOVE###

import sys
import os
import argparse


###---------------------------HELPER FUNCTIONS---------------------------------###


def generate_output_file_name(file_name):
    '''
    Takes a file name as an argument and returns string for the name of the
    output file. The generated name is based on the original file
    name.

    Specific example
    ================
    Calling function with
        ("file1.txt")
    returns
        "file1_converted_to_sys_id.txt"
    '''
    main_part_of_name, file_extension = os.path.splitext(
        file_name) #from http://stackoverflow.com/questions/541390/extracting-extension-from-filename-in-python
    if '.' in file_name:  #I don't know if this is needed with the os.path.splitext method but I had it before so left it
        return main_part_of_name + "_converted_to_sys_id" + file_extension
    else:
        return file_name + "_converted_to_sys_id"

def generate_output_file(provided_text):
    '''
    function text and saves it as a text file
    '''
    name_of_file_to_save = generate_output_file_name(list_files_to_analyze_list[0])
    data_file_stream = open(name_of_file_to_save , "w")
    data_file_stream.write(provided_text.rstrip('\r\n')) #rstrip to remove trailing newline
    # from http://stackoverflow.com/questions/275018/how-can-i-remove-chomp-a-newline-in-python
    data_file_stream.close()
    sys.stderr.write( "\nOverlap identified! Shared items list saved as '{0}'.\n".format(name_of_file_to_save))

def list2text(a_list):
    '''
    a function that takes a lists and makes a string where each item in list
    is on a new line
    '''
    return "\n".join(a_list)

###--------------------------END OF HELPER FUNCTIONS---------------------------###
###--------------------------END OF HELPER FUNCTIONS---------------------------###














#*******************************************************************************
###-----------------for parsing command line arguments-----------------------###
parser = argparse.ArgumentParser(prog='cufflinksIDs_list_to_systematic_names.py',
    description=" cufflinksIDs_list_to_systematic_names.py uses a `gtf` file \
    that is output from Cufflinks to convert a list of `XLOC_` gene ids in a \
    file to systematic yeast gene ids. The list should be gene ids each on a \
    separate line of the file. \
    **** Script by Wayne Decatur   \
    (fomightez @ github) ***")

parser.add_argument("List", help="Names of file containing `XLOC` list to fix. REQUIRED.", type=argparse.FileType('r'), metavar="FILE")

#I would also like trigger help to display if no arguments provided because need at least one input file
if len(sys.argv)==1:    #from http://stackoverflow.com/questions/4042452/display-help-message-with-python-argparse-when-script-is-called-without-any-argu
    parser.print_help()
    sys.exit(1)
args = parser.parse_args()
XLOC_list_file = args.List












###-----------------Actual Main portion of script---------------------------###



# Go through gtf file making a dictionary if the conversion factors
conversion_resolving_dictionary = {}

# make some lists to note some problem gene_ids
cuff_id_problem_list = []
XLOC_doubles_problem_list = []

gtf_file_stream = open(gtf_annotation_file_name , "r")
for line in gtf_file_stream:
    line = line.strip() # don't want line endings so I can easily
        # manipulate later, hence the use of `.strip()`
    # need to get the `XLOC` gene id and the systematic one in the line
    if 'oId' in line:
        oId = line.split('oId "')[1].split('"')[0]
        #print oId
        systematic_id = oId
        if oId.startswith('CUFF.'):
            nearest_ref = oId = line.split('nearest_ref "')[1].split('"')[0]
            # next conditional isn't perfect but should be true if doesn't match
            # most systematic yeast names. Added endswith to try and avoid
            # not catching cases where begins with "Y" but doesn't match systematic
            # name. For example, "YRA1" should be caught and "YDR190C" should not.
            # NOTE IT WILL NOT CATCH THE SYSTEMATIC NAME OF THE MITOCHONDRIAL GENES
            # WHICH ALL BEGIN IN 'Q' AND END IN NUMBERS, SUCH AS, 'Q0140' and 'Q0045.
            if not (nearest_ref.startswith('Y') and nearest_ref.endswith(("W","C","W-A","W-B","C-A","C-B"))):
                # added an additional check for the ones I had seen that are okay
                if not (nearest_ref in cleared_as_being_converted_properly):
                    cuff_id_problem_list.append (nearest_ref)
                    print "Was {0} converted properly?".format(nearest_ref)

            systematic_id = nearest_ref

        #print systematic_id
        XLOC_gene_id = line.split('gene_id "')[1].split('"')[0]
        #print XLOC_gene_id
        # Now to add to dictionary with values being systematic_id and keys being XLOC_gene_id
        # But because I have seen instance where same XLOC_gene_id for different systematic_id
        # I am going to adjust to favor non-hyphenated ones and to note any of these duplicates

        # SO SHOULD HAVE BEEN ALL THAT WAS NEEDED WAS NEXT LINE
        #conversion_resolving_dictionary [XLOC_gene_id] = systematic_id

        # BUT more informative, safe route for now
        # if already stored key, add this XLOC_gene_id as a potential problem one
        if XLOC_gene_id in conversion_resolving_dictionary:
            # append both id of current and already stored to possible problem list
            XLOC_doubles_problem_list.append(systematic_id)
            XLOC_doubles_problem_list.append(conversion_resolving_dictionary [XLOC_gene_id])
            # to disfavor hyphenated ones, change out current stored systematic_id if
            # it has a hyphen, REPLACING with one just identified
            if '-' in conversion_resolving_dictionary [XLOC_gene_id]:
                conversion_resolving_dictionary [XLOC_gene_id] = systematic_id
        else:
            conversion_resolving_dictionary [XLOC_gene_id] = systematic_id


if cuff_id_problem_list:
    sys.stderr.write( "\nGiven, those you probably should use the `geneID_list_to_systematic_names.py`.\n")
# print conversion_resolving_dictionary






# CONVERT
# Now it is time to go through the list in the file provided when command called
# and convert from `XLOC` values to systematic ids based on the conversion_resolving_dictionary

new_file_text = ""
lines_processed = 0

# open input file and start reading
sys.stderr.write("\nReading input file and converting...")


# make some lists to note some problem gene_ids
overlap_with_cuff_id_problem_list = []
overlap_with_XLOC_doubles_problem_list = []


for line in XLOC_list_file:
    lines_processed += 1
    id_to_convert = line.strip() # don't want line endings so I can easily
    # manipulate and compare later, hence the use of `.strip()`
    new_file_text = new_file_text + conversion_resolving_dictionary[id_to_convert] + '\n'
    # check for trouble ids
    if id_to_convert in cuff_id_problem_list:
        overlap_with_cuff_id_problem_list.append(id_to_convert)
    if id_to_convert in XLOC_doubles_problem_list:
        overlap_with_XLOC_doubles_problem_list.append(id_to_convert)










# Completed scan of input file and therefore close file, alert user as to any
# issues, and write new file.
sys.stderr.write( "\n"+ str(lines_processed) + " lines read from '" + XLOC_list_file.name + "'.")

# alert user as to any potential issues
if overlap_with_cuff_id_problem_list:
    alert_string = "Your list contains genes that overlap with the following genes that were questionably converted because oId began with 'CUFF.' and id designated for conversion doesn't look like an SGD systematic gene id: "
    for each_item in overlap_with_cuff_id_problem_list:
        alert_string += each_item + ", "
    alert_string += "\b\b." #backspaces to delete last comma and space
    sys.stderr.write("\n"+alert_string+"\n")
    sys.stderr.write("You should probably check those were converted correctly.\n")
    sys.stderr.write( "\nGiven those, if they aren't RNAs like snR51, you probably should use the `geneID_list_to_systematic_names.py` script next if you are looking to compare with lists that use systematic names.\n")
if overlap_with_XLOC_doubles_problem_list:
    alert_string = "Your list contains genes that overlap with the following genes that involved two `XLOC` designations: "
    for each_item in overlap_with_XLOC_doubles_problem_list:
        alert_string += each_item + ", "
    alert_string += "\b\b." #backspaces to delete last comma and space
    sys.stderr.write("\n"+alert_string+"\n")
    sys.stderr.write("You should probably check those were converted correctly\n")


output_file_name = generate_output_file_name(XLOC_list_file.name)
output_file = open(output_file_name, "w")
output_file.write(new_file_text.rstrip('\r\n')) #rstrip to remove trailing newline
# from http://stackoverflow.com/questions/275018/how-can-i-remove-chomp-a-newline-in-python
output_file.close()
sys.stderr.write("\n\nOutput file named '" + output_file_name +"' created.\n")




#*******************************************************************************
###-***********************END MAIN PORTION OF SCRIPT***********************-###
#*******************************************************************************
