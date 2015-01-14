#ConvertFASTAdnaSEQtoRNA.py  by Wayne Decatur
#ver 0.1
#
#
# To GET HELP/MANUAL, enter on command line:
# python ConvertFASTAdnaSEQtoRNA.py  --help
#
#*************************************************************************
#USES Python 2.7
# Purpose: Takes as input a file of nucleotide sequences in FASTA format
# from NCBI and changes all 'T's to 'U's in the sequence, effectively
# changing it to mRNA sequence. Importantly, it doesn't alter the description
# line of the FASTA entries, i.e., the line beginning '>' before the sequence.

# Note there is no check if these are protein or nucleic sequences. Thus it will
# change T to U in protein if you try to give it protein sequences.
#
#
#
#
# TO RUN:
# For example, enter on the command line, the line
#-----------------------------------
# python ConvertFASTAdnaSEQtoRNA.py
#-----------------------------------
#
#
#
#*************************************************************************




##################################
#  USER ADJUSTABLE VALUES        #
##################################
#
#
#
#*******************************************************************************
#*******************************************************************************












#*******************************************************************************
#*******************************************************************************
###DO NOT EDIT BELOW HERE - ENTER VALUES ABOVE###


import os
import sys
import logging
import argparse
from argparse import RawTextHelpFormatter
from Bio import Entrez
import urllib
import re
import time
#import gzip

#DEBUG CONTROL
#comment line below to turn off debug print statements
#logging.basicConfig(stream=sys.stderr, level=logging.DEBUG)


#path_to_folder_with_file = "" # LEFT HERE FOR POSSIBLE USE IN DEBUGGING
#FASTA_protein_sequence_records_file = "30test.faa" # LEFT HERE FOR USE IN DEBUGGING



###---------------------------HELPER FUNCTIONS---------------------------------###



###--------------------------END OF HELPER FUNCTIONS---------------------------###










###-----------------Actual Main function of script---------------------------###
###------------------GET FASTA FILE AND PREPARE TO PARSE---------------------###
#file to be provided as a argument when call program.
#argparser from http://docs.python.org/2/library/argparse.html#module-argparse and http://docs.python.org/2/howto/argparse.html#id1
parser = argparse.ArgumentParser(
    prog='ConvertFASTAdnaSEQtoRNA.py',description="ConvertFASTAdnaSEQtoRNA takes as input a file of nucleotide sequences in \nFASTA format from NCBI and changes all 'T's to 'U's in the sequence,\neffectively changing it to mRNA sequence.\nImportantly, it doesn't alter the description line of the FASTA entries.\n\n\nWritten by Wayne Decatur --> Fomightez @ Github or Twitter.  \n\n\nActual example what to enter on command line to run program:\npython GetmRNAforProtein.py input_file.txt\n ", formatter_class=RawTextHelpFormatter
    )
#learned how to control line breaks in description above from http://stackoverflow.com/questions/3853722/python-argparse-how-to-insert-newline-the-help-text
#DANG THOUGH THE 'RawTextHelpFormatter' setting seems to apply to all the text for argument choices. I don't know yet if that is what really what I wanted.
parser.add_argument("InputFile", help="name of data file containing a list of FASTA DNA sequences\nfor which you want mRNA sequences. REQUIRED.")
#I would also like trigger help to display if no arguments provided because need at least input file
if len(sys.argv)==1:    #from http://stackoverflow.com/questions/4042452/display-help-message-with-python-argparse-when-script-is-called-without-any-argu
    parser.print_help()
    sys.exit(1)
args = parser.parse_args()

if os.path.isfile(args.InputFile):
    #root_path = path_to_folder_with_file # LEFT HERE FOR USE IN DEBUGGING
    #fasta_file = open(root_path + FASTA_protein_sequence_records_file , "r")# LEFT HERE FOR USE IN DEBUGGING; JUST UNCOMMENT THIS AND ABOVE LINE AND COMMENTOUT NEXT LINE
    fasta_file_name = args.InputFile
    logging.debug(fasta_file_name)



    #open the input file for reading
    fasta_file = open(fasta_file_name , "r")

    #Set_up_for_ouput_file
    TheFileNameMainPart, fileExtension = os.path.splitext(fasta_file_name) #from http://stackoverflow.com/questions/541390/extracting-extension-from-filename-in-python
    if '.' in fasta_file_name:  #I don't know if this is needed with the os.path.splitext method but I had it before so left it
        output_fileName=  TheFileNameMainPart+"_mRNAconv"+fileExtension
    else:
        output_fileName= fasta_file_name +"_mRNAconv"

    #Open the output file for writing
    output_file = open(output_fileName, 'w')

    #Next..
    # step through the lines of the input file and replace T with U in sequence,
    # ignoring the decription line that begins with '>'. Make a new file
    # with the unaltered (in case of description line) line or altered sequence.

    NumberofFastaEntries = 0

    for line in fasta_file:
        line = line.strip ();#this way I have better knowledge of first character and control of ends ultimately becaue I know there isn't any unless I add
        if line.startswith('>'): # only need those that start with the greater-than symbol as they have the accessions and need to be modified
            NumberofFastaEntries += 1
            output_line = line
        else:
            #Could convert to uppercase first so no need to deal with catching
            # lower or upper case versions with line below.
            #output_line = (line.upper()).replace('T','U')

            #However, maybe some information encoded in case by user? Maybe
            # best just to leave alone. Easy in this case to handle both.
            output_line = line.replace('T','U')
            output_line = output_line.replace('t','u')


        output_file.write(output_line + "\n")




    #done with user file read in and making output file; close them
    fasta_file.close()
    output_file.close()

    #Give user some feedback
    sys.stderr.write("\nResults written to file '"+ output_fileName +"'.")




    #some additional cleaning up
    #next line is just to clean up so stdout is on next line at end
    sys.stderr.write("\n")





else:
    sys.stderr.write("SORRY. " + args.InputFile + " IS NOT RECOGNIZED AS A FILE.\n\n")
    parser.print_help()
    sys.exit(1)



#*******************************************************************************
#*******************************************************************************




