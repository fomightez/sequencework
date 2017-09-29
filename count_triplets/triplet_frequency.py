#! /usr/bin/env python

#triplet_frequency.py by Wayne Decatur
#ver 0.1
#
#
# To GET HELP/MANUAL, enter on command line:
# python triplet_frequency.py --help
#
#*************************************************************************
#USES Python 2.7
# Purpose: Takes a set of RNA sequences and analyzes the frequency (actually
# total counts) of triplet nucleotide occurences among the entire set.
#
#
# v.0.1. Started
#
# To do:
#      -deal with spaces in input - web form should have text form control to remove or leave
#             - THE BASIC DEFAULT APPROACH WILL BE TO REMOVE THEM
#      - with fasta / input format for web version (add file handling for standalone)
#
#
# TO RUN:
# For example, enter on the command line, the line
#-----------------------------------
# python triplet_frequency.py
#-----------------------------------
#
#
#*************************************************************************




##################################
#  USER ADJUSTABLE VALUES        #

##################################
#
sequences = (
    'GGUUGAACAGUUUCGGCUUACAUUUUGAUCCUAUAGUUUGCUAGCGCCCUUUUUACAUAUUGUUUACCUGGAAUU',
    'GGUCAAAUUUUUUUACUAAGAGCUUUAAUUCCAGGUUUUAUCGUCUUCCUUUCAGACGAUUUUUUCUAUAGGAUC',
    'GGAAGACGAUUUUCAAUAUGUAAUUUCUAGACAUAAUUUAAAUUUGACCUUUAUGUAAGCCGUUUGAAUGAAGCU',
    'GGGCGCUAGCUUUAAAUCGUCUGUUUAGCUUCAUUCUUUCUGUUCAACCUUUGCUCUUAGUAUUUUUAUGUCUAG',
    )  # CASE DOES NOT MATTER

standard_RNA_nucloetides = ("A","U","G","C") #should not need changing
#
#*************************************************************************
#*************************************************************************




















#*************************************************************************
#*************************************************************************
###DO NOT EDIT BELOW HERE - ENTER VALUES ABOVE###
import os
import sys
from stat import *
#from urllib import urlopen
import logging
import argparse
#from argparse import RawTextHelpFormatter

#DEBUG CONTROL
#comment line below to turn off debug print statements (or make all the debug
# statements comments)
logging.basicConfig(stream=sys.stderr, level=logging.DEBUG)




#argparser from http://docs.python.org/2/library/argparse.html#module-argparse and http://docs.python.org/2/howto/argparse.html#id1
parser = argparse.ArgumentParser(prog='triplet_frequency.py',description="triplet_frequency.py \
    takes a set of RNA sequences and analyzes the frequency of triplet nucleotide \
    occurences among the entire set.")
#parser.add_argument("InputData", help="Set of RNA sequences. REQUIRED.")
# #I would also like trigger help to display if no arguments provided because need at least input file
# if len(sys.argv)==1:    #from http://stackoverflow.com/questions/4042452/display-help-message-with-python-argparse-when-script-is-called-without-any-argu
#     parser.print_help()
#     sys.exit(1)
args = parser.parse_args()











###-----ACQUIRE NEEDED INFORMATION FROM USER IF NOT PROVIDED ON COMMAND LINE----###
processing_preamble = """The program blah blah blah <--DELETE IF NOT NEEDED.
"""

#SHOULDN'T BE MUCH HERE BECAUSE JUST NEED ONCE ARGUMENT.

###--------------------------END OF INFO ACQUISITION---------------------------###











###---------------------------HELPER FUNCTIONS---------------------------------###
def gen_all_sequences(outcomes, length):
    """
    Iterative function that enumerates the set of all sequences of
    outcomes of given length
    """

    ans = set([()])
    for dummy_idx in range(length):
        temp = set()
        for seq in ans:
            for item in outcomes:
                new_seq = list(seq)
                new_seq.append(item)
                temp.add(tuple(new_seq))
        ans = temp
    return ans

###--------------------------END OF HELPER FUNCTIONS---------------------------###




###-----------------Actual Main function of script---------------------------###
###------------------GET DATA AND PREP BASED ON IT---------------------------###
#Examine input and make sure nothing odd about it
converted_Ts = False #flag if input contained Ts changed to Us
contains_non_nts_in_seq = False #flag if input contains anything that isn't
# standard nucleotides A,G,C,U, or T. Case doesn't matter.

could_be_used_in_case_of_converted_sequences = [] #will need to collect all
# sequences, and if convert some nts from T to U will need to swap this
# for starting sequences

for sequence in sequences:
    if 't' in sequence.lower():
        converted_Ts = True
        #logging.debug(sequence)
        sequence = sequence.upper().replace('T', 'U')
        #logging.debug(sequence)

    #impt to collect all sequences since later one might be the one with
    #converted nts. Need to collect after letters changed though.
    could_be_used_in_case_of_converted_sequences.append(sequence)

    #don't want to fault lines beginning with '>' for having letters other than
    #nucleotides becase it probably means it is a fasta formatted sequence and
    #then the line with the caret at the start can have letters other than
    #nucleoctides.
    if not sequence.startswith('>'): #assume fasta format if starts with '>'
        for letter in sequence:
            if letter not in standard_RNA_nucloetides and letter.lower() != 't':
                contains_non_nts_in_seq = True

#Now swap in the converted set for sequences if any were converted.
if converted_Ts:
    #logging.debug(sequences)
    sequences = could_be_used_in_case_of_converted_sequences
    #logging.debug(sequences)

#Next, enumerate all outcomes of length 3 for standard RNA nucleotides
possible_triplets_tuples = gen_all_sequences(standard_RNA_nucloetides,3)

#FOR DEBUGGING
#logging.debug(len(possible_triplets_tuples ))
#logging.debug(possible_triplets_tuples )

# Since gen_all_sequences I copied from coursera principles of computing produces
# tuples. Convert each tuple to a string using join. from
# http://stackoverflow.com/questions/19641579/python-convert-tuple-to-string
possible_triplets = [''.join(the_triplet_tuples
    ) for the_triplet_tuples in possible_triplets_tuples]


#FOR DEBUGGING
#logging.debug(len(possible_triplets))
#logging.debug(possible_triplets)




# Next make a dictionary with all those triplets initialized with zero as each value
triplet_freq_dict = {key: 0 for key in possible_triplets} #used dict comprehension
# see http://stackoverflow.com/questions/2241891/how-to-initialize-a-dict-with-keys-from-a-list-and-empty-value-in-python
# Seeed even easier than dict.fromkeys or method used by cookie_clicker_help.py

#FOR DEBUGGING
#logging.debug(triplet_freq_dict)

# Scan through each sequence updating the dictionary of values according to the
# frequency of that triplet in sequences.
for sequence in sequences:
    for triplet in possible_triplets:
        #logging.debug(triplet)
        #logging.debug(sequence)
        sequence = sequence.lower() # Just in case it is in caps or mixed case since counting
        # triplets on next line in lowercase.
        triplet_freq_dict[triplet] += sequence.count(triplet.lower())
#logging.debug(triplet_freq_dict)

# Return the triplet frequencies
### adapted from Alex Martelli's answer on
#### http://stackoverflow.com/questions/991350/counting-repeated-characters-in-a-string-in-python
for each in sorted(triplet_freq_dict, key=triplet_freq_dict.get, reverse=True):
  print '%s %6d' % (each, triplet_freq_dict[each])

#GiveFeedback on possible issues
if converted_Ts:
    sys.stderr.write("\nT's in the sequence have been converted to U's for counting.\n")
if contains_non_nts_in_seq:
    sys.stderr.write("\n***NOTICE*****************NOTICE*****************NOTICE*** \n")
    sys.stderr.write("Please be aware the sequence \
examined contains letters\nother than standard nucleotides.\n\
IT MAY SIGNAL THAT THE TEXT PROVIDED AND ANALYZED WAS NOT\n\
EXACTLY WHAT YOU INTENDED. Of course, it may not.\n\
The characters not corresponding to standard nucleotides\n\
were not considered, but it is being mentioned here\nin case the input \
possessing additional characters\nwas possibly an oversight and not your intention.\n")
    sys.stderr.write("**END OF NOTICE*********END OF NOTICE*********END OF NOTICE** \n \n")

