#!/usr/bin/env python


# systematic_names_to_standard_names_using_cufflinks_gtf.py by Wayne Decatur



#*******************************************************************************
# USES Python 2.7 but should be convertable via 2to3, see https://docs.python.org/3.0/library/2to3.html
#
# PURPOSE: Using a `gtf` file that is output from Cufflinks, converts a list of
# systematic yeast gene ids to standard (common names), where they exist.
# The list should be gene ids each on a separate line of the file.
# It was written for yeast genes but may work with a gtf annotation file for any
# organism if the tags involved in the annotation file include `oId` and
# `gene_id` and , possibly, `nearest_ref`.
# This was written as sort of a quick and dirty way to do this since it would
# be better to not require a gtf file. Ended up not being so 'quick' because
# there were some special cases I had built into the backbone I had built it on
# , see ` cufflinksIDs_list_to_systematic_names.py` for that backbone code,
# and it needed some special cases add itself, such as when making
# the conversion dictionary from the gtf file, don't place as a standard name
# one that looks like  systematic name as long as the current supposed systematic
# name looks actually not like a systematic name. This was to deal with a
# handful of edge cases, like ICR1, RUF5-1, RUF5-2, snR73, etc., where the GTF
# had systematic names as the `gene_name` value and more traditional `standard
# (common) names` in the `oId` position where the systematic name typically was.
#
#
#
#
# see end of
# `Trying to convert Shafi's cufflinks list to something that works with other lists with systematic names.md`
# for info about this script
#
#
# Dependencies beyond the mostly standard libraries/modules:
#
#
#
# VERSION HISTORY:
# v.0.1. basic working version.
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


# NOT BEING USED HERE. IT WAS FOR `cufflinksIDs_list_to_systematic_names.py`
cleared_as_being_converted_properly = ["HRA1", "PWR1", "Q0045", "tA(AGC)K2", "snR51" ] #there were some issues with
# some genes having 'CUFF.' as text being used for conversion and
# I used `nearest_ref` for those, but several didn't look like systematic yeast
# names even with using those. I checked and these ones confirmed converted fine
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
        "file1_converted_to_stf_id.txt"
    '''
    main_part_of_name, file_extension = os.path.splitext(
        file_name) #from http://stackoverflow.com/questions/541390/extracting-extension-from-filename-in-python
    if '.' in file_name:  #I don't know if this is needed with the os.path.splitext method but I had it before so left it
        return main_part_of_name + "_converted_to_std_id" + file_extension
    else:
        return file_name + "_converted_to_std_id"

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


def looks_like_systematic_name(id):
    '''
    taks a string of text and checks if looks like a yeast SGD systematic gene
    name.
    That is (for all but mitochondrial as best I know)
    begins with "Y"
    has "L" or "R" as the third character
    ends in "W","C","W-A","W-B","C-A","C-B"

    For example, "YRA1" should be `false` and "YDR190C" should be `true`.

    An issue is that it won't test true for mitochondrial genes like `VAR1`
    where  `Q0140` is the systematic name as shown at
    http://www.yeastgenome.org/locus/S000007275/overview
    Doesn't matter for what I want here but NOTE IT ISN'T COVERING ALL CASES YET.


    '''
    return (id.startswith("Y")) and (
        id[2] == "L" or id[2] == "R") and (
        id.endswith(("W","C","W-A","W-B","C-A","C-B")))
    # note a tuple of choices can be an argument of 'endswith'

###--------------------------END OF HELPER FUNCTIONS---------------------------###
###--------------------------END OF HELPER FUNCTIONS---------------------------###














#*******************************************************************************
###-----------------for parsing command line arguments-----------------------###
parser=argparse.ArgumentParser(prog="systematic_names_to_standard_names_using_cufflinks_gtf",
    description=" systematic_names_to_standard_names_using_cufflinks_gtf.py uses a `gtf` file \
    that is output from Cufflinks to convert a list of systematic gene ids in a \
    file to standard (common) gene names, where they exist. The list should be \
    gene ids each on a separate line of the file. \
    **** Script by Wayne Decatur   \
    (fomightez @ github) ***")

parser.add_argument("List", help="Names of file containing `systematic ids` list to convert. REQUIRED.", type=argparse.FileType('r'), metavar="FILE")

#I would also like trigger help to display if no arguments provided because need at least one input file
if len(sys.argv)==1:    #from http://stackoverflow.com/questions/4042452/display-help-message-with-python-argparse-when-script-is-called-without-any-argu
    parser.print_help()
    sys.exit(1)
args = parser.parse_args()
sys_ids_list_file = args.List












###-----------------Actual Main portion of script---------------------------###



# Go through gtf file making a dictionary if the conversion factors
conversion_resolving_dictionary = {}

# make some lists to note some problem gene_ids
cuff_id_problem_list = []
name_doubles_problem_list = []

gtf_file_stream = open(gtf_annotation_file_name , "r")
for line in gtf_file_stream:
    line = line.strip() # don't want line endings so I can easily
        # manipulate later, hence the use of `.strip()`
    # need to get the systematic one and the standard (common) one in the line
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
            if not (looks_like_systematic_name(nearest_ref)):
                # added an additional check for the ones I had seen that are okay
                if not (nearest_ref in cleared_as_being_converted_properly):
                    cuff_id_problem_list.append (nearest_ref)
                    print "Was {0} converted properly?".format(nearest_ref)

            # don't bother with nearest_ref if already have that as a key
            # in conversion dictionary because chances are it is one where the
            # oID began "CUFF." for one and already collected systematic name and standard name
            if (not nearest_ref in conversion_resolving_dictionary):
                systematic_id = nearest_ref

        #print systematic_id
        standard_gene_name = line.split('gene_name "')[1].split('"')[0]
        #print standard_gene_name

        # However, for some of entries in GTF file the oId names are actually
        # standard names and the gene_name name is systematic name. So for those
        # cases,it would be best to just take oId name as standard name too.
        # Also don't want a non-systematic-looking sys id if it begins with 'CUFF.'
        # Don't need to add anything specifically saying genes where systematic
        # name beings with 'Q' should skip this branch because the gene_name
        # position (standard_gene_name variable) for those will never look like
        # a typical yeast systematic gene that begeins with 'Y', has third letter
        # that is "R" or "L", etc.
        # Two steps:
        # STEP# 1: test for those cases
        if (not systematic_id.startswith('CUFF.')) and (
            not looks_like_systematic_name(systematic_id)) and (
            looks_like_systematic_name(standard_gene_name)):
            # STEP #2: In those cases, set the name harvested from the oId name
            # as the standard name in the conversion dictionary
            conversion_resolving_dictionary [systematic_id] = systematic_id
        else:
            # Now to add to dictionary with keys being systematic_id and values being standard_gene_name
            # But because when developing the script that went from XLOC id to systematic ID
            # I saw I have seen instance where same XLOC_gene_id for different systematic_id
            # I am going to be conservative and check for some issues I saw then
            # to adjust to favor non-CUFF ones and note any of these duplicates

            # SO SHOULD HAVE BEEN ALL THAT WAS NEEDED WAS NEXT LINE
            #conversion_resolving_dictionary [systematic_id] = standard_gene_name

            # BUT more informative, safe route for now
            # if already stored systematic id as a key, add this
            # standard_gene_name as a potential problem one.
            # Except Disregard those that begin with "CUFF" or if systematic_id = standard name
            # or , IMPORTANTLY,
            # the ones where the one currently collected equals the one already
            # stored.
            if (systematic_id in conversion_resolving_dictionary):
                # Will add as problem if also doesn't meet a few criteria as
                # mentioned above. There is no `else`. So if not added to
                # possible problem list, nothing is done with those cases for now.
                if (not systematic_id.startswith('CUFF')) and (
                    systematic_id != standard_gene_name) and (
                    standard_gene_name != conversion_resolving_dictionary [systematic_id]):
                    #print systematic_id
                    # append both id of current and already stored to possible problem list
                    name_doubles_problem_list.append(conversion_resolving_dictionary [systematic_id])
                    name_doubles_problem_list.append(systematic_id)
                # to disfavor ones starting with 'CUFF, change out current stored systematic_id
                if conversion_resolving_dictionary [systematic_id].startswith('CUFF'):
                    conversion_resolving_dictionary [systematic_id] = standard_gene_name
            else:
                conversion_resolving_dictionary [systematic_id] = standard_gene_name


if cuff_id_problem_list:
    sys.stderr.write( "\nGiven, those you probably should use the `geneID_list_to_systematic_names.py`.\n")
# print conversion_resolving_dictionary






# CONVERT
# Now it is time to go through the list in the file provided when command called
# and convert from systematic ids to standard (common) gene name based on the
# conversion_resolving_dictionary

new_file_text = ""
lines_processed = 0

# open input file and start reading
sys.stderr.write("\nReading input file and converting...")
input_file_stream = open(sys_ids_list_file, "r")


# make some lists to note some problem gene_ids
overlap_with_cuff_id_problem_list = []
overlap_with_name_doubles_problem_list = []


for line in input_file_stream:
    lines_processed += 1
    id_to_convert = line.strip() # don't want line endings so I can easily
    # manipulate and compare later, hence the use of `.strip()`

    # if it is a tRNA meaning it begins with lowercase 't', then it might not
    # have an entry in the gtf file with XLOCs and so just pass its value to new
    # file
    if id_to_convert.startswith('t'):
        new_file_text = new_file_text + id_to_convert + '\n'
    else:
        # want string "None" if key doesn't exist for the few cases
        # where other gene listing sources had genes not present in the gtf
        #, such as `tE(CUC)D`, `'YBL005W-A`, see http://stackoverflow.com/questions/9285086/access-dict-key-and-return-none-if-doesnt-exist
        # And so "None" is default if no key for command below
        converted = conversion_resolving_dictionary.get(id_to_convert, "None")
        new_file_text = new_file_text + converted + '\n'
    # check for trouble ids
    if id_to_convert in cuff_id_problem_list:
        overlap_with_cuff_id_problem_list.append(id_to_convert)
    if id_to_convert in name_doubles_problem_list:
        overlap_with_name_doubles_problem_list.append(id_to_convert)










# Completed scan of input file and therefore close file, alert user as to any
# issues, and write new file.
sys.stderr.write( "\n"+ str(lines_processed) + " lines read from '" + sys_ids_list_file.name + "'.")
input_file_stream.close()

# alert user as to any potential issues
if overlap_with_cuff_id_problem_list:
    alert_string = "Your list contains genes that overlap with the following genes that were questionably converted because oId began with 'CUFF.' and id designated for conversion doesn't look like an SGD systematic gene id: "
    for each_item in overlap_with_cuff_id_problem_list:
        alert_string += each_item + ", "
    alert_string += "\b\b." #backspaces to delete last comma and space
    sys.stderr.write("\n"+alert_string+"\n")
    sys.stderr.write("You should probably check those were converted correctly.\n")
    sys.stderr.write( "\nGiven those, if they aren't RNAs like snR51, you probably should use the `geneID_list_to_systematic_names.py` script next if you are looking to compare with lists that use systematic names.\n")
if overlap_with_name_doubles_problem_list:
    alert_string = "Your list contains genes that overlap with the following genes that involved two systematic id designations: "
    for each_item in overlap_with_name_doubles_problem_list:
        alert_string += each_item + ", "
    alert_string += "\b\b." #backspaces to delete last comma and space
    sys.stderr.write("\n"+alert_string+"\n")
    sys.stderr.write("You should probably check those were converted correctly\n")


output_file_name = generate_output_file_name(sys_ids_list_file.name)
output_file = open(output_file_name, "w")
output_file.write(new_file_text.rstrip('\r\n')) #rstrip to remove trailing newline
# from http://stackoverflow.com/questions/275018/how-can-i-remove-chomp-a-newline-in-python
output_file.close()
sys.stderr.write("\n\nOutput file named '" + output_file_name +"' created.\n")




#*******************************************************************************
###-***********************END MAIN PORTION OF SCRIPT***********************-###
#*******************************************************************************
