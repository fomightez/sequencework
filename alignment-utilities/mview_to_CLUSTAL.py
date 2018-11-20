#!/usr/bin/env python
# mview_to_CLUSTAL.py 
__author__ = "Wayne Decatur" #fomightez on GitHub
__license__ = "MIT"
__version__ = "0.1.0"


# mview_to_CLUSTAL.py  by 
# Wayne Decatur
# ver 0.1
#
#*******************************************************************************
# Verified compatible with both Python 2.7 and Python 3.6; written initially in 
# Python 3. 
#
#
# PURPOSE: Takes a multiple sequence alignment output by MView with options set 
# so mainly just the lines of sequence and not html or consenuss information and
# reformats it to CLUSTAL.
# What settings are needed for MView are shown in the demonstration notebook 
# for the script. (Also detailed in the docstring of main function of script.)
# About CLUSTAL format, see 
# http://meme-suite.org/doc/clustalw-format.html and
# http://scikit-bio.org/docs/0.4.2/generated/skbio.io.format.clustal.html
# and http://scikit-bio.org/docs/0.4.2/generated/skbio.io.format.clustal.html
# and https://www.ebi.ac.uk/Tools/msa/clustalw2/help/faq.html#18. A good 
# example can be viewed at http://wwwabi.snv.jussieu.fr/public/Clustal2Dna/clustal.html .
#
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
#
#
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
# to do:
#
#
#
#
# TO RUN:
# Examples,
# Enter on the command line of your terminal, the line
#-----------------------------------
# python mview_to_CLUSTAL.py mview_out.aln
#-----------------------------------
#
# Issue `mview_to_CLUSTAL.py  -h` for details.
# 
#
#
# To use this after importing/pasting or loading into a cell in a Jupyter 
# notebook, specify at least the results file (or results as a string) in the 
# call to the main function similar to below:
# mview_to_CLUSTAL("mview_out.txt")
# 
# 
#
'''
CURRENT ACTUAL CODE FOR RUNNING/TESTING IN A NOTEBOOK WHEN IMPORTED/LOADED OR 
PASTED IN ANOTHER CELL:
mview_to_CLUSTAL("mview_out.aln")
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


suffix_for_saving = "_clustalized.clw" #to be used for naming the output 
# automatically when running script from command line to act on an input file


#
#*******************************************************************************
#**********************END USER ADJUSTABLE VARIABLES****************************













#*******************************************************************************
#*******************************************************************************
###DO NOT EDIT BELOW HERE - ENTER VALUES ABOVE###

import sys
import os





###---------------------------HELPER FUNCTIONS---------------------------------###


def generate_output_file_name(file_name,suffix_for_saving):
    '''
    Takes a file name as an argument and returns string for the name of the
    output file. The generated name is based on the original file
    name.


    Specific example
    =================
    Calling function with
        ("mview_out.txt")
    returns
        "mview_out_clustalized.clw"
    '''
    main_part_of_name, file_extension = os.path.splitext(
        file_name) #from 
    #http://stackoverflow.com/questions/541390/extracting-extension-from-filename-in-python
    if '.' in file_name:  #I don't know if this is needed with the os.path.splitext method but I had it before so left it
        return main_part_of_name + suffix_for_saving
    else:
        return file_name + suffix_for_saving



###--------------------------END OF HELPER FUNCTIONS---------------------------###
###--------------------------END OF HELPER FUNCTIONS---------------------------###






#*******************************************************************************
###------------------------'main' function of script---------------------------##

def mview_to_CLUSTAL(mview_data, suffix_for_saving = suffix_for_saving):
    '''
    Main function of script. 
    It will take an alignment text file for MView (output with 'special' settings)
    and make it CLUSTAL.

    The 'special' settings for MView are essentially, in addition to choosing 
    'MVIEW' as `OUTPUT FORMAT`, 'ON' for `ALIGNMENT`, the desired width, most 
    settings are adjusted to 'OFF', in particular `HTML MARKUP`, `RULER`, and 
    `CONSENSUS`.

    Mainly it:
    - Discards the provided header line, if it is there.
    - Discards the provided footer line, if it is there.
    - Add header line that just starts `CLUSTAL`
    - Collects the IDs and the sequences
    - Truncates the ids to 30 if not already. (Based on 
    https://www.ebi.ac.uk/Tools/msa/clustalw2/help/faq.html#18)
    - To make things easy, add spaces to each indentifier to make them all 
    equal length with longest.
    - If length of longest identifier is less than 16 then left justify with 
    spaces to sixteen and then put sequence on each line at next 
    column  (because MUSCLE and example at 
    http://wwwabi.snv.jussieu.fr/public/Clustal2Dna/clustal.html seem to do 17th 
    column starts sequence). Otherwise, add six spaces to each equalized id and 
    then put sequence for that line.
    '''
    # Bring in the necessary MView output data:
    #---------------------------------------------------------------------------
    # Read entire file into memory
    with open(mview_data, 'r') as mviewfile:
        mviewlines=mviewfile.read().split('\n')


    # feedback
    sys.stderr.write("MView output read...")


    # Process header and footer, if there
    #---------------------------------------------------------------------------
    # The header or footer lines might not be there depending on what user 
    # actually copied and sotry to continue to work even if not there.
    if mviewlines[0].startswith("Reference sequence"):
        del mviewlines[0]
    if mviewlines[-1].startswith("MView "):
        del mviewlines[-1]

    # Add CLUSTAL header line and empty line after it
    #---------------------------------------------------------------------------
    header_to_add = ("CLUSTAL multiple sequence alignment by mview_to_CLUSTAL "
        "({})\n\n".format(__version__))
    output_str = header_to_add 

    # Process each line of actual alignment, parsing IDs and sequences
    #---------------------------------------------------------------------------
    # With the header and footer removed, the only lines with content should 
    # contain the id and sequence for each line at index 1 and index 4, 
    # respectively, if split line by default.
    # However, will need to determine the longest identifier so cannot just add
    # to output. Going to store each line data in a dictionary with line index
    # as key. Each value for lines with content will be a list with the id 
    # first and sequence for that line second. For empty lines, the value will
    # simply be contents of line with linebreak added.
    # While at it, truncate any identitiers long than 30 characters & collect 
    # information about identifiers to be used when making them consistent & 
    # ready for output.
    mviewlines_dict = {}
    identifiers = []
    for indx, line in enumerate(mviewlines):
        if len(line.strip()) > 2:
            line_parts = line.split()
            mviewlines_dict[indx] = [line_parts[1][:30], line_parts[4]+ "\n"]
            identifiers.append(line_parts[1][:30],)
        else:
            mviewlines_dict[indx] = line + "\n"
    sys.stderr.write("collected identifiers and sequences...")

    # Make IDs consistent and ready for output.
    #---------------------------------------------------------------------------
    # To make things easy, add spaces 
    # to each indentifier to make them all equal length with longest. Based on
    # https://stackoverflow.com/a/5676676/8508004
    identifiers = set(identifiers)
    len_longest_id = len(max(identifiers, key=len))
    for i in mviewlines_dict:
        if type(mviewlines_dict[i]) == list:
            mviewlines_dict[i][0] = mviewlines_dict[i][0].ljust(len_longest_id)

    # Make output string.
    #---------------------------------------------------------------------------
    # Work through line by line. For those lines where items stored as list
    # of id and sequence for that line, justify the id and add whitespace before
    # sequence how appropriate. If length of longest identifier is less than 
    # 16 then left justify with spaces to sixteen and then put sequence on each 
    # line at next column, i.e., the 17th  (because MUSCLE and example at 
    # http://wwwabi.snv.jussieu.fr/public/Clustal2Dna/clustal.html seem to do 
    # 17th column starts sequence). Otherwise, add six spaces to each equalized
    # id and then put sequence for that line.
    # Note that I use a range to iterate over the keys of mviewlines_dict 
    # because important to do line by line, starting with first. For Python 3, 
    # that should be order of dict already but not certain it will be for 2.7
    # and so doing that way will insure order no matter what version of Python.
    for n in range(len(mviewlines_dict)):
        if type(mviewlines_dict[n]) == list:
            if len_longest_id < 16: 
                the_line = (mviewlines_dict[n][0].ljust(16) + 
                    mviewlines_dict[n][1])
            else:
                the_line = (mviewlines_dict[n][0] + "      " + 
                mviewlines_dict[n][1])
            output_str += the_line 
        else:
            output_str += mviewlines_dict[n]
    sys.stderr.write("arranging for output...")




    
    # Save the modified version as a file.
    #---------------------------------------------------------------------------
    out_alignment_name = generate_output_file_name(mview_data,suffix_for_saving)
    # prepare output file for saving so it will be open and ready
    with open(out_alignment_name, 'w') as output_file:
        output_file.write(output_str)


    # Feedback
    sys.stderr.write("\n\nAlignment converted from MView to CLUSTAL "
        "saved as '{}'.".format(out_alignment_name))
    sys.stderr.write("\nFinished.\n")







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
    kwargs['suffix_for_saving'] = suffix_for_saving
    mview_to_CLUSTAL(mview_data,**kwargs)
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
        'mview_to_CLUSTAL.py ',
        description="mview_to_CLUSTAL.py  \
        takes a multiple sequence alignment output by MView with options set \
        so mainly just the lines of sequence and not html or consenuss \
        information and reformats it to CLUSTAL. What settings are needed for \
        MView are shown in the demonstration notebook for the script.\
        **** Script by Wayne Decatur   \
        (fomightez @ github) ***")

    parser.add_argument("mview", help="Name of file containing multiple sequence \
        alignment from MView output. REQUIRED. \
        ", metavar="MVIEW")

    parser.add_argument('-os', '--output_suffix', action='store', type=str, 
    default= suffix_for_saving, help="OPTIONAL: Set a suffix for including in file \
    name of output. \
    If none provided, '{}' will be used.".format(suffix_for_saving))



    #I would also like trigger help to display if no arguments provided because 
    # need at least one for url
    if len(sys.argv)==1:    #from http://stackoverflow.com/questions/4042452/display-help-message-with-python-argparse-when-script-is-called-without-any-argu
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    mview_data = args.mview
    suffix_for_saving = args.output_suffix


    main()

#*******************************************************************************
###-***********************END MAIN PORTION OF SCRIPT***********************-###
#*******************************************************************************
