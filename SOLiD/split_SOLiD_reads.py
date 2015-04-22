#split_SOLiD_reads.py  by Wayne Decatur
#ver 0.1
#
#
# To GET HELP/MANUAL, enter on command line:
# python split_SOLiD_reads.py  --help
#
#*************************************************************************
#USES Python 2.7
# Purpose:
# Separates out reads that match a barcode in Applied Biosystems SOLiD
# data, such as that in Rhee HS and Pugh BF, 2011*.

# For example, the metada associated with AB SOLiD System 3.0 reads for data
# SRR346368, run #1 at http://www.ncbi.nlm.nih.gov/sra/SRX098212[accn]),
# says the fastq file contains both Reb1 and Gal4 data, each having a diffferent
# barcode for  the  first 6  bp of the reverse mate of the paired-end reads.
# In the fastq file the pairs are interleaved and  so this program handles
# that type of data.
#
# *Full Citation:
# Rhee HS and Pugh BF. 2011. Comprehensive genome-wide protein-DNA
# interactions detected at single-nucleotide resolution.
# Cell. 2011 Dec 9;147(6):1408-19. doi: 10.1016/j.cell.2011.11.013.
# PMID: 22153082 http://www.ncbi.nlm.nih.gov/pubmed/22153082.
# Reb1- and Gal4-associated data downloaded
# http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?cmd=dload&run_list=SRR346368&format=fastq
# empoloying this approach -->
# http://user.list.galaxyproject.org/does-Galaxy-have-a-tool-for-converting-sra-files-to-fastq-files-td4655179.html
#
#
# LIMITATIONS:
# Unlike other programs for demultiplexing read data, this one is rather naive
# and does not allow for mismatches or bp in front of the barcode. See
# https://www.biostars.org/p/82513/ for examples that do.
#
# Dependencies:
# Nothing but the fairly standard modules such as os, sys, and argparse.
#
#
#
# v.0.1. Started
#
# To Do: The -I option of SRA toolkit insures a read id is indicated after
# each after spot id as `accession.spot.readid` on defline of the read infor and
# I'd like to add checking that value matches so that doesn't collect any where
# the sync got knocked off between the the tag and the proper genomic read mate.
#
#
# TO RUN:
# Adjust the 'USER ADJUSTABLE VALUES' to match the settings you need.
#
# then, enter on the command line, the line
#-----------------------------------
# python split_SOLiD_reads.py FASTQ_FILE_NAME BARCODE BARCODE_CONTAINING_READ_LENGTH
#-----------------------------------
#
# where FASTQ_FILE_NAME is the name of the file of reads and BARCODE is
# the sequence of the barcode. (For the Gal4 data in Rhee and Pugh, 2011 this is
# TAGCGT.) BARCODE_CONTAINING_READ_LENGTH is the length of
# the reverse read containing the barcode. (In the case of Rhee and Pugh, 2011
# this is 10.)
#
# So for Rhee and Pugh, 2011 Gal4 data.
#-----------------------------------
# python split_SOLiD_reads.py sraSRR34638.fastq TAGCGT 10
#-----------------------------------
#
# Example input and output using above command:
# INPUT:
#
# @SRR346368.1.1 0176_20090623_2_Specific_Factors_926_144_1410 length=35
# ACATAACTAATTTACCGCTTTTTAAACATACCCCC
# +SRR346368.1.1 0176_20090623_2_Specific_Factors_926_144_1410 length=35
# #&&####+$*%#$##%##%$%$*'%#$)&#%%+#%
# @SRR346368.1.2 0176_20090623_2_Specific_Factors_926_144_1410 length=0
#
# +SRR346368.1.2 0176_20090623_2_Specific_Factors_926_144_1410 length=0
#
# @SRR346368.122860.1 0176_20090623_2_Specific_Factors_947_1888_393 length=35
# AGGGATGGCAAAAAAAGGGCAATCTAAAGATAAAG
# +SRR346368.122860.1 0176_20090623_2_Specific_Factors_947_1888_393 length=35
# (6<7<)?+-&')%,:-+551#+$3='-(%84,-2&
# @SRR346368.122860.2 0176_20090623_2_Specific_Factors_947_1888_393 length=10
# TGTTTCCCCC
# +SRR346368.122860.2 0176_20090623_2_Specific_Factors_947_1888_393 length=10
# '$.<4)<1;<
# @SRR346368.4889.1 0176_20090623_2_Specific_Factors_946_46_1420 length=35
# CAAAAGAGAAAGACAAAGGCCGGGGGAAAAGAAAA
# +SRR346368.4889.1 0176_20090623_2_Specific_Factors_946_46_1420 length=35
# =:59??;<8@?=7=1<<'8=<4.>52227;,<74;
# @SRR346368.4889.2 0176_20090623_2_Specific_Factors_946_46_1420 length=10
# TAGCGTTCTC
# +SRR346368.4889.2 0176_20090623_2_Specific_Factors_946_46_1420 length=10
# 9;6+8898<:
#
#
# OUTPUT:
#
# @SRR346368.4889.1 0176_20090623_2_Specific_Factors_946_46_1420 length=35
# CAAAAGAGAAAGACAAAGGCCGGGGGAAAAGAAAA
# +SRR346368.4889.1 0176_20090623_2_Specific_Factors_946_46_1420 length=35
# =:59??;<8@?=7=1<<'8=<4.>52227;,<74;

#
#*************************************************************************


##################################
#  USER ADJUSTABLE VALUES        #
##################################
#

#
#*******************************************************************************
#*******************************************************************************










#*******************************************************************************
#*******************************************************************************
###DO NOT EDIT BELOW HERE - ENTER ANY VALUES ABOVE###

import os
import sys
import logging
import argparse
from argparse import RawTextHelpFormatter
#import urllib
#import re


#DEBUG CONTROL
#comment line below to turn off debug print statements
#logging.basicConfig(stream=sys.stderr, level=logging.DEBUG)








###---------------------------HELPER FUNCTIONS---------------------------------###


def devise_file_name(file_name, barcode):
    '''
    Takes a file name as an argument and returns string for the output file
    name based on the original file name, adding the barcode as part of the
    new name.

    Specific example
    ================
    Calling function with
        ("sra_dataSRR346368.fastq", "TAGCGT")
    returns
        sra_dataSRR346368_TAGCGT.fastq
    '''
    main_part_of_name, file_extension = os.path.splitext(
        file_name) #from http://stackoverflow.com/questions/541390/extracting-extension-from-filename-in-python
    if '.' in file_name:  #I don't know if this is needed with the os.path.splitext method but I had it before so left it
        return main_part_of_name +"_" + barcode + file_extension
    else:
        return file_name+"_" + barcode



def extract_reads(fastq_input_file_name, the_barcode, the_barcode_seq_len):
    '''
    This function takes a file containing AB SOLiD seqence reads as interleaved
    paired-ends and sorts out those with the matching barcode. Saving them to a
    new file.

    Arguments for the function are as follows:
        * the fastq file to scan
        * the barcode corresponding to the reads sought
        * the length of the sequence that contains the barcode

    For example, in Rhee and Pugh, the last value is 10 bp because in the case
    of SRR346368, run #1 at http://www.ncbi.nlm.nih.gov/sra/SRX098212[accn],
    involving Gal4 and Reb1 data, the reverse read of the
    mated read pairs is 10 bases while the main genome forward read is 35 bases.

    The function returns the following:
        * number_of_paired_reads
        * number_of_reads_collected
        * the name of the created output_file
    '''
    # in preparation use the fastq_input_file_name and barcode to generate
    # a name for the output file
    output_file_name = devise_file_name(fastq_input_file_name, the_barcode)


    #initialize values
    the_FASTQ_file = open(fastq_input_file_name , "r")
    the_new_file = open(output_file_name , "w")
    paired_reads_tally = 0
    reads_collected_tally = 0
    # initialize the storage list with values so always has 5 strings in it
    # even at start of file reading
    previous_five_lines_storage_list = ["nada","zip","zero","nilch","zed"]
    #step through file mining the information
    for line in the_FASTQ_file:
        line = line.strip ();#this way I have better control of ends ultimately
        # ecaue I know there isn't any unless I add something.

        # Since want to tally paired-reads for feedback to user, I'll count
        # when I see the symbol '@' at the beginning of a line in the file
        # since first line (sequence identifier and an optional description)
        # of each read in fastq format begins with that symbol. See
        # http://en.wikipedia.org/wiki/FASTQ_format
        # Count as 0.5 since counting pairs. Concerned with pairs here because
        # in Rhee and Pugh, 2011 Gal4 data, the short read of the mate pairs is
        # is mainly for the barcode and the long read was used for mapping.
        if line.startswith('@'):
            paired_reads_tally += 0.5

        # Need to identify those starting with barcode but that aren't
        # the longer reads. So want the read to match size of the barcode
        # -contaning reads as provided by user.
        if (line.startswith(the_barcode)) and (
            len(line) == the_barcode_seq_len) :
            # Now want to collect the read info from lines 1 through 4 of the
            # previous five lines & write it to file, plus
            # add 1 to the collected reads count collected count
            for index, item in enumerate(previous_five_lines_storage_list[:4]):
                the_new_file.write(item + "\n")
            reads_collected_tally += 1

        # Prepare for reading in next line by shifting stored lines information.
        # First, add the new line to the stored lines
        previous_five_lines_storage_list.append(line)
        # Next, delete the first element stored to reset to only storing last 5.
        del previous_five_lines_storage_list[0]
        assert (len(previous_five_lines_storage_list) == 5), (
            "previous_five_lines_storage_list should not contain more than \
            five items after preparing for reading next line of file.")

    #Completed scan of input file and therefore close files and return results.
    the_FASTQ_file.close()
    the_new_file.close()
    return (paired_reads_tally, reads_collected_tally, output_file_name)


###--------------------------END OF HELPER FUNCTIONS---------------------------###










###-----------------Actual Main function of script---------------------------###
###----------------------GET FILE AND PREPARE TO PARSE-----------------------###
#file to be provided as a argument when call program.
#argparser from http://docs.python.org/2/library/argparse.html#module-argparse and http://docs.python.org/2/howto/argparse.html#id1
parser = argparse.ArgumentParser(
    prog='split_SOLiD_reads.py',description="split_SOLiD_reads.py is designed to demultiplex based on barcode AB SOLiD fastq\nfiles with interleaved paired-ends reads, such as from Rhee and Pugh, 2011.\n\n\nWritten by Wayne Decatur --> Fomightez @ Github or Twitter.  \n\nActual example what to enter on command line to run program:\npython split_SOLiD_reads.py sraSRR34638.fastq TAGCGT 10\n ", formatter_class=RawTextHelpFormatter
    )
#learned how to control line breaks in description above from http://stackoverflow.com/questions/3853722/python-argparse-how-to-insert-newline-the-help-text
#DANG THOUGH THE 'RawTextHelpFormatter' setting seems to apply to all the text for argument choices. I don't know yet if that is what really what I wanted.
parser.add_argument("InputFile", help="name of file containing barcoded reads in fastq format. REQUIRED.")
parser.add_argument("barcode", help="sequence of barcode of interest. REQUIRED.")
parser.add_argument("barcode_seq_len", help="Length of sequence which contains barcode. REQUIRED.", type=int)
#I would also like trigger help to display if no arguments provided because need at least input file
if len(sys.argv)==1:    #from http://stackoverflow.com/questions/4042452/display-help-message-with-python-argparse-when-script-is-called-without-any-argu
    parser.print_help()
    sys.exit(1)
args = parser.parse_args()

if os.path.isfile(args.InputFile):
    #root_path = path_to_folder_with_file # LEFT HERE FOR USE IN DEBUGGING
    #fastq_file = open(root_path + FASTA_protein_sequence_records_file , "r")# LEFT HERE FOR USE IN DEBUGGING; JUST UNCOMMENT THIS AND ABOVE LINE AND COMMENTOUT NEXT LINE
    fastq_file = args.InputFile
    logging.debug(fastq_file)
    number_of_paired_reads = 0 #initiate with zero as value of number of paired reads analyzed
    number_of_reads_collected = 0 #initiate with zero as value of reads collected
    barcode = args.barcode
    barcode_seq_len = args.barcode_seq_len

    #Read FASTQ read entries, keeping track of total, and write those of
    # interest to a file. THIS FUNCTION CALL IS THE MAIN POINT OF THIS PROGRAM.
    sys.stderr.write("Reading in your FASTQ file...")
    number_of_paired_reads, number_of_reads_collected, output_file = (
        extract_reads(fastq_file, barcode, barcode_seq_len))

    #FOR DEBUGGING
    logging.debug(number_of_paired_reads)
    logging.debug(number_of_reads_collected)
    logging.debug(output_file)

    #give user some stats and feeback
    sys.stderr.write("\nConcluded. \n"+ str(number_of_paired_reads)+
        " read pairs analyzed. " + "A grand total of "+
        str(int(number_of_paired_reads * 2))+" reads.\n" +
        str(number_of_reads_collected)+" reads collected based on barcode.")
    sys.stderr.write("\nThe file "+output_file+" has been created in same directory as the input file.\n\n")


else:
    sys.stderr.write("SORRY. " + args.InputFile + " IS NOT RECOGNIZED AS A FILE.\n\n")
    parser.print_help()
    sys.exit(1)


#*******************************************************************************
#*******************************************************************************
