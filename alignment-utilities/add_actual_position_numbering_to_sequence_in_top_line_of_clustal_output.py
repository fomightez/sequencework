#!/usr/bin/env python
# add_actual_position_numbering_to_sequence_in_top_line_of_clustal_output.py
__author__ = "Wayne Decatur" #fomightez on GitHub
__license__ = "MIT"
__version__ = "0.1.0"


# add_actual_position_numbering_to_sequence_in_top_line_of_clustal_output.py by 
# Wayne Decatur
# ver 0.1
#
#*******************************************************************************
# Verified compatible with both Python 2.7 and Python 3.6; written initially in 
# Python 3. 
#
#
# PURPOSE: Takes a text document of an alignment in CLUSTAL format and adds 
# position indicators to the top line that is the actual numbering of the 
# residues in contiguous sequence. The added annotation makes it easier to 
# locate regions in the alignment if you know the location in the top sequence, 
# typically the reference sequence.
#
# Mview program at EMBL-EBI https://www.ebi.ac.uk/Tools/msa/mview/ will take
# a clustal formatted alignment and add numbers, but they will be the total 
# numbers of characters on the the line, including the dashes that indicate gaps
# and so they will simply count up 80 at a time. For locating residues relative
# the numbering in the contiguous sequence, it isn't very helpful for when there
# are large gaps, such as when aligning an enture chromosome of divergent 
# organisms.
# RELATED:
# Clustal Omega program at EMBL-EBI  https://www.ebi.ac.uk/Tools/msa/clustalo/
# generates alignment with character counting for each sequence in alignment if
# default `ClustalW with character counts` accepted. An additional note is that
# if you want to keep your reference sequence on the top, open `More options...` 
# and change to `ORDER` setting > `input` to keep the sequences in the alignment
# in the order if the input, instead of allowing the program to alter the order.
# PLUS, http://www.bioinformatics.org/sms2/group_dna.html will add numbers to
# FASTA format and FAST format can include gaps.
#
#
#
#
#
# If you are Wayne, see `Collecting yeast XXXX XXXXX (XXXXX) May 2018.md` for 
# impetus behind this script.
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
#
#
#
#
# TO RUN:
# Examples,
# Enter on the command line of your terminal, the line
#-----------------------------------
# python add_actual_position_numbering_to_sequence_in_top_line_of_clustal_output.py ALIGNMENT_TEXT_FILE
#-----------------------------------
#
# Issue `add_actual_position_numbering_to_sequence_in_top_line_of_clustal_output.py -h` for details.
# 
#
#
# To use this after importing/pasting or loading into a cell in a Jupyter 
# notebook, specify at least the results file (or results as a string) in the 
# call to the main function similar to below:
# add_actual_position_numbering_to_sequence_in_top_line_of_clustal_output("test.out")
# 
#-or-
# To specify file name to save as, specify `output_name` like on next line:
# add_actual_position_numbering_to_sequence_in_top_line_of_clustal_output("test.out",output_name="align_adjusted.txt")
#
#
# 
#
'''
CURRENT ACTUAL CODE FOR RUNNING/TESTING IN A NOTEBOOK WHEN IMPORTED/LOADED OR 
PASTED IN ANOTHER CELL:
add_actual_position_numbering_to_sequence_in_top_line_of_clustal_output("test.out")
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

output_name = 'alignment_adjusted.clustal' #to be used by default when running 
# in a Jupyter cell or imported, otherwise one based on alignment file is 
# generated automatically.

suffix_for_saving = "_ADJUSTED" #to be used for naming the output automatically
# when running script from command line to act on an input file

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
    =================
    Calling function with
        ("alignment.clustal")
    returns
        "alignment_ADJUSTED.clustal"
    '''
    main_part_of_name, file_extension = os.path.splitext(
        file_name) #from http://stackoverflow.com/questions/541390/extracting-extension-from-filename-in-python
    if '.' in file_name:  #I don't know if this is needed with the os.path.splitext method but I had it before so left it
        return main_part_of_name + suffix_for_saving  + file_extension
    else:
        return file_name + suffix_for_saving + ".clustal"

###--------------------------END OF HELPER FUNCTIONS---------------------------###
###--------------------------END OF HELPER FUNCTIONS---------------------------###



#*******************************************************************************
###------------------------'main' function of script---------------------------##

def add_actual_position_numbering_to_sequence_in_top_line_of_clustal_output(
    alignment, output_name = output_name):
    '''
    Main function of script. 
    It will take an alignment text file (or the alignment as a Python string) 
    and add position indicators to the top line that is the actual numbering of 
    the residues in contiguous sequence.

    The option to provide the results as a string is to handle where sending the
    data directly from shell script to Python without a typical file 
    intermediate, see the advanced notebook at https://git.io/vpr7i for 
    examples. The obvious use case for that is when working in the Jupyter 
    # notebook environment.
    '''
    # Bring in the necessary data:
    #---------------------------------------------------------------------------

    try:
        with open(alignment, 'r') as the_file:
            alignment = the_file.read()
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
    sys.stderr.write("Alignment read...")

    # Identify the first identifier in each alignment block
    first_words = []
    for line in alignment.split("\n"):
        if line:
            first_word = line.split()[0]
            if first_word in first_words:
                first_id = first_word
                break
            else:
                first_words.append(first_word)
    # feedback
    sys.stderr.write("top line identifier determined as '{}'...".format(first_id))

    # Go through each line and if the first line identifier recognized, add a
    # line with the growing count to the building output and then the line with
    # the identifier
    growing_ouput = []
    current_position = 0
    for line in alignment.split("\n"):
        if line.startswith(first_id):
            seq = line.split(first_id)[1].strip()
            char_num_without_gaps = len(seq.replace("-", ""))
            end_pos = current_position + char_num_without_gaps
            if char_num_without_gaps:
                current_position += 1 #because if line entirely gap,don't add 
                # on, but if anything on line it will be at least one above what
                # was on previous line
            if end_pos < 0:
                end_pos = 0
            spacing = len(seq) -  len(str(end_pos))
            if len(seq) > (len(str(current_position)) + len(str(end_pos))):
                num_text = '{:<{}}{}'.format(current_position,spacing,end_pos) #
                #based on https://stackoverflow.com/a/5676884/8508004 and comments
            else:
                if spacing < len(str(end_pos)):
                    spacing = len(seq) +  len(str(end_pos)) 
                num_text = '{:>{}}'.format(end_pos,spacing)
            start_spacing = line.index(seq)
            num_line = (" "* start_spacing) + num_text
            growing_ouput.append(num_line)
            growing_ouput.append(line.strip())
            current_position = end_pos
        else:
            growing_ouput.append(line.strip())
    output = "\n".join(growing_ouput)

        

    # Reporting and Saving
    #---------------------------------------------------------------------------
    with open(output_name, 'w') as output_file:
        output_file.write(output)
    # feedback
    sys.stderr.write("\nAnnotated alignment saved to '{}'.".format(output_name))



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
    kwargs['output_name'] = output_name
    add_actual_position_numbering_to_sequence_in_top_line_of_clustal_output(alignment,**kwargs)
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
        'add_actual_position_numbering_to_sequence_in_top_line_of_clustal_output.py',
        description="add_actual_position_numbering_to_sequence_in_top_line_of_clustal_output.py \
        takes output from a text document of an alignment in CLUSTAL format \
        and adds position indicators to the top line that is the actual \
        numbering of the residues in contiguous sequence. The added annotation \
        makes it easier to locate regions in the alignment if you know the \
        location in the top sequence, typically the reference sequence.\
        **** Script by Wayne Decatur   \
        (fomightez @ github) ***")

    parser.add_argument("align_file", help="Name of file of alignmnet text \
        file to add numbers according to the top line.\
        ", metavar="ALIGNMENT_FILE")

    parser.add_argument('-o', '--output', action='store', type=str, 
    default= output_name, help="OPTIONAL: Set file name for saving output. \
    If none provided, '{}' will be inserted in front of the file extension \
    of the alignment file name.".format(suffix_for_saving))



    #I would also like trigger help to display if no arguments provided because 
    # need at least one for url
    if len(sys.argv)==1:    #from http://stackoverflow.com/questions/4042452/display-help-message-with-python-argparse-when-script-is-called-without-any-argu
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    alignment = args.align_file
    if args.output == output_name:
        output_name = generate_output_file_name(alignment)
    else:
        output_name = args.output


    main()

#*******************************************************************************
###-***********************END MAIN PORTION OF SCRIPT***********************-###
#*******************************************************************************
