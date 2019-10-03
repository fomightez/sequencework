#!/usr/bin/env python
# matches_a_patmatch_pattern.py
__author__ = "Wayne Decatur" #fomightez on GitHub
__license__ = "MIT"
__version__ = "0.1.0"


# matches_a_patmatch_pattern.py by Wayne Decatur
# ver 0.1
#
#*******************************************************************************
# Verified compatible with both Python 2.7 and Python 3.7; written initially in 
# Python 3. 
#
#
# PURPOSE: Takes a sequence pattern in PatMatch syntax, and checks if a provided 
# sequence (file name for sequence in FASTA format or string text can be 
# provided) contains a match. It reports True or False depending on that
# assessment. Optionally, it can be restricted to checking if the provided 
# sequence is a match to a sequence pattern. Note that the point, i.e., 
# checking that the provided sequence matches entirely to the pattern and is 
# therefore a specific instance of the more general pattern, was the 
# original impetus for writing this script. However, it seemed I might write a
# more general script and then include an option to do that.
# When the residue type is nucleic, both strands will be searched.
# Only concerned with first sequence if a multi-sequence file is provided. Loop
# over sequences and then call script with each if you need to do multiple.
#
# Currently, designed with comparing short sequence strings in mind. For 
# PatMatch syntax, see https://www.yeastgenome.org/nph-patmatch#examples .
# For those familiar with prepration of files for use with PatMatch, you don't 
# have to apply `unjustify.pl` yourself to your sequences. That preparation of 
# removing line endings will all be handled internally meaning 
# STANDARD FASTA FILE FORMAT IS ACCEPTABLE. 
#
# Note that with 'match_over_entirety', the pattern itself must be same length 
# in characters as the expected matching sequence, and therefore it is 
# recommended to avoid using it with fancier patterns using symbols as such 
# patterns will likely cause this not to work correctly. If using such fancier 
# search patterns, you'll want to screen the length of the sequence realtive the 
# possible outcomes prior and combine the results of the assessment of match 
# into your interpretation of what the report means. Or you'd need to edit the 
# script.
#
# Meant to be a 'yes' or 'no' answer. If you want to know about the location of 
# the match or matches you'll want to use `patmatch_results_to_df.py` from 
# https://github.com/fomightez/sequencework/tree/master/patmatch-utilities
# It actually uses that script to do the work here.
#
# Note this is conceptually, vaguely related to my script 
# `find_sequence_element_occurrences_in_sequence.py` found in 
# https://github.com/fomightez/sequencework/tree/master/FindSequence except that
# handles Regular Expression (REGEX) syntax, and here I want to be able to 
# handle more biological contexts without having to add many sets to the regex 
# code. By using PatMatch, I get the biological context without having to 
# re-implement it.
#
#
#
#
#
#
#
# Meant to be used intially in classifying many sequences obtained in some 
# searches where often nearby sequence elements were used in searching.
#
#
#
# Dependencies beyond the mostly standard libraries/modules:
# sh, pyfaidx
#
#
#
# VERSION HISTORY:
# v.0.1. basic working version

#
# To do:
# - verify actually works for Python 2 <--- I know I have BLAST with Python 2. Do I have PatMatch too?
#
#
#
#
# TO RUN:
# Examples,
# Enter on the command line of your terminal, the line
#-----------------------------------
# python matches_a_patmatch_pattern.py -n DDWDWTAWAAGTARTADDDD ATTGATATAAGTAATAGATA
#-----------------------------------
# Issue `matches_a_patmatch_pattern.py -h` for 
# details.
# 
#
#
#
# To use, in a Jupyter notebook or IPython console:
# from matches_a_patmatch_pattern import matches_a_patmatch_pattern
# matches_a_patmatch_pattern("DDWDWTAWAAGTARTADDDD","ATTGATATAAGTAATAGATA","nucleic",match_over_entirety=True);
#
# 
#
'''
CURRENT ACTUAL CODE FOR RUNNING/TESTING IN A NOTEBOOK WHEN IMPORTED:
from matches_a_patmatch_pattern import matches_a_patmatch_pattern
matches_a_patmatch_pattern("DDWDWTAWAAGTARTADDDD","ATTGATATAAGTAATAGATA","nucleic", match_over_entirety=True);
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




#
#*******************************************************************************
#**********************END USER ADJUSTABLE VARIABLES****************************






















#*******************************************************************************
#*******************************************************************************
###DO NOT EDIT BELOW HERE - ENTER VALUES ABOVE###

import sys
import os
import pandas
try:
    from pathlib import Path
except ImportError:
    from pathlib2 import Path
import tempfile 
import subprocess


###---------------------------HELPER FUNCTIONS--------------------------------###

from contextlib import contextmanager,redirect_stderr,redirect_stdout
from os import devnull

@contextmanager
def suppress_stdout_stderr():
    """
    A context manager that redirects stdout and stderr to devnull.
    From https://stackoverflow.com/a/52442331/8508004
    """
    with open(devnull, 'w') as fnull:
        with redirect_stderr(fnull) as err, redirect_stdout(fnull) as out:
            yield (err, out)

###--------------------------END OF HELPER FUNCTIONS--------------------------###
###--------------------------END OF HELPER FUNCTIONS--------------------------###

#*******************************************************************************
###------------------------'main' function of script--------------------------##

def matches_a_patmatch_pattern( pattern, sequence, residue_type,
    match_over_entirety = False):
    '''
    Takes a sequence pattern in PatMatch syntax, and checks if a provided 
    sequence (provided as a file_name of a FASTA file or string text of 
    sequence) contains a match. It reports True or False depending on 
    that assessment. Optionally, with 'match_over_entirety' it can be restricted 
    to checking if the provided sequence is a match to a sequence pattern. Note 
    that with 'match_over_entirety', the pattern itself must be same length in 
    characters as the expected matching sequence, and therefore avoid using it 
    with fancier patterns using symbols as such patterns will likely cause this 
    not to work right.

    Also need to specify residue type of sequence provided. Nucleic or protein.
    When the residue type is nucleic, both strands will be searched.

    Only concerned with first sequence if a multi-sequence file is provided. 
    Loop over sequences and then call this function with each if you need to 
    examine multiple sequences.
    '''

    # Quick check if the `match_over_entirety` option is used
    #--------------------------------------------------------------------
    # 
    # If the user provided sequences of unequal length, they may want to know
    # that isn't what it is for, and so in addition to reporting False, give 
    # feedback.
    # It is straightforward to check of the length of the two strings (or 
    # sequence in the file) and so might as well exclude doing all the 
    # additional steps if clear that size isn't same.
    if match_over_entirety:
        if os.path.isfile(sequence):
            from pyfaidx import Fasta
            sequence_to_check = Fasta(sequence)[0]
            sequence_to_check = str(sequence_to_check) #This has the added 
            # bonus that pyfaidx handles removing all line endings from sequence
        else:
            sequence_to_check = "".join(sequence.splitlines()) # I suspect if 
            # providing a sequence as a string, line breaks won't be included 
            # but this removes them just in case. Based on 
            # https://stackoverflow.com/a/52088151/8508004 .
        if len(pattern) != len(sequence_to_check):
            sys.stderr.write("****Error??*** You called the script with "
                "the `match_over_entirety` option;\nhowever, your provided "
                "sequences are not equal in size and thus don't\nmatch over "
                "their entirety. Feel free to ignore this concern.\n\n")
            if __name__ == "__main__":
                print("False")
                return
            else:
                return False

    # Verify PatMatch installed in environment.
    #--------------------------------------------------------------------
    # 
    sys.stderr.write("Checking PatMatch software present"
        "...\n")
    # Tried the following but couldn't get a non-error with reasonable return 
    # because no version reported or anything. Maybe just check the script is 
    # there?
    '''
    try:
        cmd="perl patmatch.pl"
        result = subprocess.check_output(cmd, shell=True) # based on 
        # https://stackoverflow.com/a/18739828/8508004 
        result = result.decode("utf-8") # so it isn't bytes
        if "Must give residue type" in result:
            sys.stderr.write("Detected patmatch...\n")
    except CalledProcessError:  #or should it be `subprocess.CalledProcessError` like in https://stackoverflow.com/a/29824645/8508004
        sys.stderr.write("\nPatMatch not detected. Please include PatMatch software or "
        "run in an\nenvironment launched from "
        "https://github.com/fomightez/patmatch-binder.\n**EXITING !!**.\n")
        sys.exit(1)
    '''
    from sh import find
    patmatch_path = find(".","-type","f","-name","patmatch.pl")
    if patmatch_path:
        patmatch_path = patmatch_path.split()[0].strip()
    if patmatch_path == "":
        # try going up in hierarchy to home directory (based on 
        # https://stackoverflow.com/a/4028943/8508004 )
        home = str(Path.home())
        patmatch_path = find(home,"-type","f","-name","patmatch.pl")
        if patmatch_path:
            patmatch_path = patmatch_path.split()[0].strip()

    if patmatch_path == "":
        sys.stderr.write("\nPatMatch not detected. Please add PatMatch "
        "software or "
        "run in an\nenvironment launched from "
        "https://github.com/fomightez/patmatch-binder.\n**EXITING !!**.\n")
        sys.exit(1)



    # Get the script containing the function needed for sending PatMatch data
    # to pandas dataframe and import the function
    #--------------------------------------------------------------------
    #
    # Because I want to build in a way to allow this to be done manually
    # check if already there first and then otherwise try to get.
    file_needed = "patmatch_results_to_df.py"
    if not os.path.isfile(file_needed):
        sys.stderr.write("Obtaining script containing function to pass "
            "PatMatch results to Python"
            "...\n")
        # based on http://amoffat.github.io/sh/
        from sh import curl
        curl("-OL",
            "https://raw.githubusercontent.com/fomightez/sequencework/"
            "master/patmatch-utilities/patmatch_results_to_df.py")
        # verify that worked & ask for it to be done manually if fails
        if not os.path.isfile(file_needed):
            github_link = ("https://github.com/fomightez/sequencework/blob/"
                "master/patmatch-utilities/patmatch_results_to_df.py")
            sys.stderr.write("\n'patmatch_to_df.py' not found. Please add it to "
            "your current working\ndirectory from {}"
            ".\n**EXITING !!**.\n".format(github_link))
            sys.exit(1)
    from patmatch_results_to_df import patmatch_results_to_df

    # Note I was going to automate detection of residue type but then realized
    # that if someone wanted to search a pattern that is purely glycine, 
    # cysteine, threonine, and alanine, they'd be out of luck.


    # Create a temp file of the sequence.
    #--------------------------------------------------------------------
    #
    # If a file name is provided, then the temp file is the sequence without the 
    # line endings. Meant to be equivalent to running the `unjustify.pl` script 
    # for how PatMatch is typically used. If a sequuence string is provided then
    # the tempfile is made with that so that PatMatch can access it since
    # PatMatch needs a file.
    # `tempfile.mkstemp` based on 
    # https://github.com/caltechads/deployfish/blob/master/deployfish/mysql.py .
    # Had to use that so perl could access it; Python's `tempfile.TemporaryFile` 
    # gets deleted upon closure.
    if os.path.isfile(sequence):
        # use pyfaidx to get the sequence. This has the added bonus that pyfaidx 
        # handles removing all line endings from sequence.
        from pyfaidx import Fasta
        seq_record = Fasta(sequence)[0]
        tmphandle, tmppath = tempfile.mkstemp(text=True,suffix=".fa")
        with open(tmphandle, 'w') as output:
            output.write(">seq\n{}".format(str(seq_record)))
    else:
        tmphandle, tmppath = tempfile.mkstemp(text=True,suffix=".fa")
        with open(tmphandle, 'w') as output:
            output.write(">seq\n{}".format("".join(sequence.splitlines()))) # I 
            # suspect if providing a sequence as a string, line breaks won't be 
            # included but `"".join(sequence.splitlines())` handles removing 
            # them, just in case. Based on 
            # https://stackoverflow.com/a/52088151/8508004 .
    prepared_sequence = tmppath


    # Format the command and run
    #--------------------------------------------------------------------
    # 
    if residue_type == "nucleic":
    	type_setting = "c"
    elif residue_type == "protein":
    	type_setting = "p"
    cmd = 'perl {} -{} {} {}'.format(
    	patmatch_path,type_setting,pattern,prepared_sequence)
    result = subprocess.check_output(cmd, shell=True)
    #result = subprocess.check_output('perl {} -{} {} {}'.format(
    #   patmatch_path,type_setting,pattern,StringIO(sequence)), shell=True)
    result = result.decode("utf-8")
    sys.stderr.write("Sending the PatMatch results to Python"
        "...\n")
    with suppress_stdout_stderr():
        pm_df = patmatch_results_to_df(result)
    #print(pm_df) #for debugging


    # Clean up temp file (prepared sequence file) now that it has been used
    #----------------------------------------------------------------------
    #
    os.remove(prepared_sequence)


    #OLD TESTING I WILL BE ABLE TO RUN patmatch LATER since I had problem earlier when not really running it
    # # BASICALLY ALL WORKED ONCE I FIXED `patmatch_path`! <-- that was issue!
    #cmd = 'perl {} -c "CAAAGGAA" patmatch_1.2/test/ATH1_cdna_test.prepared 1 ids'.format(patmatch_path)
	#result = subprocess.check_output(cmd, shell=True)
	#result = result.decode("utf-8") 
    # THIS WORKED BEFORE I FIXED `patmatch_path` because it didn't use it. <-- `patmatch_path` was issue!
    #cmd = "perl ./patmatch_1.2/patmatch.pl -c CAAAGGAA patmatch_1.2/test/ATH1_cdna_test.prepared 1 ids"
    #res = subprocess.check_output(cmd, shell=True)
    #res = res.decode("utf-8") 
    #print(res)
    # THIS WORKED but requires me to save a file
    #cmd = 'perl {} -c "CAAAGGAA" patmatch_1.2/test/ATH1_cdna_test.prepared 1 ids > another_test.txt'.format(patmatch_path)
    #subprocess.run(cmd, shell=True) 
    # NEXT COMMAND WORKS TOO!
    #subprocess.call(["sh", "-c", "perl ./patmatch_1.2/patmatch.pl -c CAAAGGAA patmatch_1.2/test/ATH1_cdna_test.prepared 1 ids >results.txt"]) #based on https://stackoverflow.com/a/29264627/8508004 works but `subprocess.check_output` wasn't working as written on next few lines
    # Okay, so saving a file of results works but how to do it without a file intermediate!?!?
    # NEXT TRY ---> https://stackoverflow.com/questions/18822036/python-get-output-of-child-process
    # BELOW THIS DIDN'T WORK - returned non-zero exit status 127<--- LATER I DETERMINED BECAUSE I WAS NOT POST-PROCESSING `patmatch_path` from a 'sh.RunningCommand' object to correct string properly!!
    #-------------------
    #cmd = 'perl {} -c "CAAAGGAA" patmatch_1.2/test/ATH1_cdna_test.prepared 1 ids'.format(patmatch_path)
    #result = subprocess.check_output(cmd, shell=True) # based on 
    # https://stackoverflow.com/a/18739828/8508004
    #result = result.decode("utf-8") # so it isn't bytes
    #print(result)
    #cmd = 'perl {} -c "CAAAGGAA" patmatch_1.2/test/ATH1_cdna_test.prepared 1 ids'.format(patmatch_path)
    #result = subprocess.check_output(cmd, shell=True) # based on 
    # https://stackoverflow.com/a/18739828/8508004
    #result = result.decode("utf-8") 



    # Analyze the results and report
    #--------------------------------------------------------------------
    #
    # Provide feedback  using `if __name__ == "__main__"` to customize 
    # if returned via Python or command line depending if script called from 
    # command line.
    if not match_over_entirety:
        if len(pm_df)>=1:
            call =  True
        elif pm_df.empty:
            call = False
    else:
        if len(pm_df)==1 and len(
            pm_df.iloc[0]["matching pattern"]) == len(pattern):
            call = True
        else:
            call = False
    sys.stderr.write("Reporting: {}.\n\n".format(call))
    if __name__ == "__main__":
        print(call)
        return
    else:
        return call




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
    kwargs['match_over_entirety'] = match_over_entirety
    matches_a_patmatch_pattern(pattern,sequence,
        residue_type ,**kwargs)
    # using https://www.saltycrane.com/blog/2008/01/how-to-use-args-and-kwargs-in-python/#calling-a-function
    # to build keyword arguments to pass to the function above
    # (see https://stackoverflow.com/a/28986876/8508004 and
    # https://stackoverflow.com/a/1496355/8508004 
    # (maybe https://stackoverflow.com/a/7437238/8508004 might help too) for 
    # related help). Makes it easy to add more later.





if __name__ == "__main__":
    ###-----------------for parsing command line arguments-------------------###
    import argparse
    parser = argparse.ArgumentParser(prog='matches_a_patmatch_pattern.py',
        description="matches_a_patmatch_pattern.py \
        Takes a sequence pattern in PatMatch syntax, and checks if a provided \
        sequence contains a match. It reports True or False depending \
        on that assessment. Optionally, it can be restricted to checking if \
        the provided sequence is a match to a sequence pattern.\
        **** Script by Wayne Decatur   \
        (fomightez @ github) ***")

    
    parser.add_argument("pattern", help="Sequence pattern in PatMatch syntax. \
        For example, to search for a S. cerevisiae mitochondrial promoter, \
        provide `DDWDWTAWAAGTARTADDDD`, without \
        any quotes or backticks, in the call to the script.\
        ", metavar="PATTERN")

    parser.add_argument("sequence", help="Filename for FASTA sequence file or \
        text of sequence string to examine for presence of the pattern. **Only \
        the first sequence of a multi-entry FASTA file is considered.**\
        ", metavar="SEQUENCE")

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-n','--nucleic',action='store_true',help='Specify the \
        residue type as nucleic. MUST SPECIFY A RESIDUE TYPE.')
    group.add_argument('-p','--protein',action='store_true',help='Specify the \
        residue type as protein. MUST SPECIFY A RESIDUE TYPE.')


    parser.add_argument("-moe", "--match_over_entirety",help=
        "add this flag to force the match test to have the sequence have to \
        match over an entirety to the pattern. In other words, force testing \
        if the sequence is a specific example of the pattern. **CAUTION: The \
        pattern itself must be same length in characters as the expected \
        matching sequence, and therefore it is recommended you not use this \
        with fancier patterns using symbols as it will likely cause this not \
        to work correctly.**",
        action="store_true")




    #I would also like trigger help to display if no arguments provided because 
    # need at least one for url
    if len(sys.argv)==1:    #from http://stackoverflow.com/questions/4042452/display-help-message-with-python-argparse-when-script-is-called-without-any-argu
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    pattern = args.pattern
    sequence = args.sequence
    match_over_entirety = args.match_over_entirety
    if args.nucleic:
        residue_type = "nucleic"
    elif args.protein:
        residue_type = "protein"





    main()

#*******************************************************************************
###-***********************END MAIN PORTION OF SCRIPT***********************-###
#*******************************************************************************
