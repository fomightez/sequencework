#! /usr/bin/env python

#add_source_organism_info_to_FASTA.py  by Wayne Decatur
#ver 0.1
#
#
# To GET HELP/MANUAL, enter on command line:
# python add_source_organism_info_to_FASTA.py  --help
#
#*************************************************************************
#USES Python 2.7
# Purpose:
# Adds source organism information in brackets to end of descrption line of each
# FASTA-formatted record for protein or nucleic acid in a file. The source
# organism information added will be the genus and species plus any other
# information such as strain and variety associated with that sequence entry.
# The sequences in the provided file must be all of the same kind, either all
# nucleic or all protein. The script will try to detect which one
# based on the first sequence in the provided file. Subsequently, the script
# will fail if the sequences in the file are a mix of protein and nucleic.
# It will leave the description line untouched if it already has some text
# occuring between brackets as it assumes that is the source organism designation
# as is typical in many protein entry description lines of FASTA-formatted
# sequences at NCBI.
# Adding source organism in brackets to description line for mRNAs is important
# because for example even though the protein entry
# http://www.ncbi.nlm.nih.gov/protein/685422211?report=fasta
# has the source organism in between brackets the mRNA entry does not. See it at
# http://www.ncbi.nlm.nih.gov/nuccore/685422210?report=fasta
# The addition of this will make the produced FASTA-formatted entries more
# consistent with those of proteins and more easily parsed by scripts.
#
# As an example of making it more easily parsed, I point out that running this
# script on a standard FASTA formatted list of nucleotide
# sequences would prepare that list of sequences for use with my
# `compare_organisms_in_two_files_of_fasta_entries.py` script.
#
# Dependencies:
# Biopython
# In addition to the fairly standard modules such as os, sys, and argparse.
#
#
#
# v.0.1. Started
#
# To Do:
# - create demo notebook in https://github.com/fomightez/cl_sq_demo-binder 
#
#
# TO RUN:
# enter on the command line, the line
#-----------------------------------
# python add_source_organism_info_to_FASTA.py file.fa
#-----------------------------------
#
# where file.fa is the name of the file of fasta records.
#
# Example input and output using above command:
# INPUT:
'''
>gi|685422210|ref|XM_009231696.1| Gaeumannomyces graminis var. tritici R3-111a-1 hypothetical protein partial mRNA
ATGACGACGGACAGCACCAACCGCGCCTCGCCGCTGGCCAGGCTCGAGCAGGCCCGGCAGCTGCTCGCCA
...
GCTCGCTGATGTAG
>gi|47075178|ref|NM_211569.1| Ashbya gossypii ATCC 10895 AGL160Wp (AGOS_AGL160W) mRNA, complete cds
ATGTCGGACAAGGCTTTGCGCGCTGGTGAGGATGGCACGGAGATCCGTAATGCGCTTCGCAGCCTACAGC
...
CGACATTTTCGGCTGA
'''
# (where '...' represents sequence not shown for simplicty sake.)
#
#
# OUTPUT:
'''
>gi|685422210|ref|XM_009231696.1| Gaeumannomyces graminis var. tritici R3-111a-1 hypothetical protein partial mRNA [Gaeumannomyces graminis var. tritici R3-111a-1]
ATGACGACGGACAGCACCAACCGCGCCTCGCCGCTGGCCAGGCTCGAGCAGGCCCGGCAGCTGCTCGCCA
...
GCTCGCTGATGTAG
>gi|47075178|ref|NM_211569.1| Ashbya gossypii ATCC 10895 AGL160Wp (AGOS_AGL160W) mRNA, complete cds [Ashbya gossypii ATCC 10895]
ATGTCGGACAAGGCTTTGCGCGCTGGTGAGGATGGCACGGAGATCCGTAATGCGCTTCGCAGCCTACAGC
...
CGACATTTTCGGCTGA
'''
#
#

#
#*************************************************************************


##################################
#  USER ADJUSTABLE VALUES        #
##################################
#
#
User_Email = "wdecatur@yahoo.com" #PUT YOUR E-MAIL HERE or NCBI's SERVER WILL COMPLAIN
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
from Bio import Entrez
from Bio import Seq
from Bio import SeqIO
#import urllib
#import re


#DEBUG CONTROL
#comment line below to turn off debug print statements
#logging.basicConfig(stream=sys.stderr, level=logging.DEBUG)








###---------------------------HELPER FUNCTIONS---------------------------------###


def generate_output_file_name(file_name):
    '''
    Takes a file name as an argument and returns string for the name of the
    output file. The generated name is based on the original file
    name.

    Specific example
    ================
    Calling function with
        ("sequences.fa")
    returns
        "sequences_with_source.fa"
    '''
    main_part_of_name, file_extension = os.path.splitext(
        file_name) #from http://stackoverflow.com/questions/541390/extracting-extension-from-filename-in-python
    if '.' in file_name:  #I don't know if this is needed with the os.path.splitext method but I had it before so left it
        return main_part_of_name + "_with_source" + file_extension
    else:
        return file_name + "_with_source"


def guess_if_nucleic(seq, thresh = 0.90, nucleic_letters = ['G', 'A', 'T', 'C', 'U']):
    """Guess if the given sequence is nucleic acid.

    It's considered nucleic if more than 90% of the sequence is GATCs (or Us).
    The threshold is configurable via the thresh parameter. nucleic_letters can
    be used to configure which letters are considered DNA; for instance,
    adding N might be useful if you are expecting data with ambiguous bases.

    function based on code found at
    http://biopython.org/pipermail/biopython/2009-January/004868.html
    As noted there it will fail for pathological cases like short amino acids
    with lots of Gly, Ala, Cys, or Thrs
    """
    if isinstance(seq, Seq.Seq):
        seq = seq.data
    elif isinstance(seq, type("")) or isinstance(seq, type(u"")):
        seq = str(seq)
    else:
        raise ValueError("Do not know provided type: %s" % seq)
    seq = seq.upper()
    nucleic_alpha_count = 0
    for letter in nucleic_letters:
        nucleic_alpha_count += seq.count(letter)
    if (len(seq) == 0 or float(nucleic_alpha_count) / float(len(seq)) >= thresh):
        return True
    else:
        return False


def extract_GI_number(string_fragment):
    '''
    uses string split to extract GI number for typical FASTA formatted sequence

    based on http://biopython.org/pipermail/biopython/2010-January/006143.html

    CAVEAT:
    ***NOTE SINCE AS POSTED AT ftp://ftp.ncbi.nih.gov/genbank/gbrel.txt GI numbers
    ARE BEING PHASED OUT, I IMAGINE USE NOT VIABLE LONG TERM.
    KEEP EYES OPEN. ***
    '''
    return string_fragment.split("|")[1]


def make_list_of_unique_identifiers (fasta_file):
    '''
    takes a file of sequence records in fasta form and reads each one extracting
    the unique identifier to use later with ELink to get Genbank formatted entry
    for sequence. Returns them as a list.

    Also it returns best gues if first entry is nucleic or protein to be later
    used to access correct database.

    Want GI number as unique indentifier because ELink works with that. I'd have
    to add another step to get it to work with accession.version identifier. I
    originally couldn't get it with the indentifier that has the letters like CAC82173.
    See code block #3 of 'Using ELink to convert protein sequences to mRNA.ipynb'
    (publically available at https://www.biostars.org/p/64078/) where
    I did it starting with the accession.version identifier. However, if have full
    Genbank FASTA record, just easier to get from that.

    CAVEAT:
    ***NOTE SINCE AS POSTED AT ftp://ftp.ncbi.nih.gov/genbank/gbrel.txt GI numbers
    ARE BEING PHASED OUT, I IMAGINE ELink MAY EVENUTALLY WORK WITH
    Accession.Version identifiers??? (Or does it already and I missed how?)
    KEEP EYES OPEN. ***

    '''
    uid_list = []
    nucleic_or_protein = "NOT YET ASSESSED"
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        seq_info  = str(seq_record.id)
        the_gi_number = extract_GI_number(seq_info)
        uid_list.append(the_gi_number)
        # for first sequence determine if seems to be nucleic or protein
        if nucleic_or_protein == "NOT YET ASSESSED":
            the_sequence = str(seq_record.seq) #based on sections 2.4.2  and 3.4 of http://biopython.org/DIST/docs/tutorial/Tutorial.html
            if (guess_if_nucleic(the_sequence)):
                nucleic_or_protein = "nucleic"
            else:
                nucleic_or_protein = "protein"
    return uid_list, nucleic_or_protein



def EPost_uids(unique_ids, db_to_use):
    '''
    Sends a list of unique identifiers to the NCBI Entrez history server in
    preparation for fetching step, using correct database.

    Returns the two items needed to access them off the server (webenv, query_key)
    subsequently.
    '''
    #EPost
    epost_handle = Entrez.epost(db=db_to_use, id=",".join(unique_ids))
    epost_result = Entrez.read(epost_handle)
    epost_handle.close()

    webenv = epost_result["WebEnv"]
    query_key = epost_result["QueryKey"]
    return webenv, query_key



def EFetch_sequences_GB_FORMAT(
    uids_list_length, the_webenv, the_query_key, db_to_use):
    '''
    Fetches the sequences in GENBANK format for a list posted prior to the
    NCBI Entrez history server.

    I am trying to make it so it is more general and handles both nucleotide or
    protein and so I have added a way to specify `db_to_use`. For nucleic acid,
    it would be `nuccore`. For protein it would be... ????

    Returns those records as a list.
    '''
    #EFetch
    batch_size = 20
    the_records = ""
    the_records_list = []
    for start in range(0, uids_list_length, batch_size):
        end = min(uids_list_length, start + batch_size)
        sys.stderr.write("Fetching records %i thru %i..." % (start + 1, end))
        fetch_handle = Entrez.efetch(db=db_to_use,
                                     rettype="gb", retmode="text",
                                     retstart=start, retmax=batch_size,
                                     webenv=the_webenv,
                                     query_key=the_query_key)
        data = fetch_handle.read()
        fetch_handle.close()
        the_records = the_records + data
    # Genbank entries START WITH THE TEXT `LOCUS`
    the_records_list = filter(None, "!LOCUS".join(the_records.split('LOCUS')).split('!')) # from http://www.gossamer-threads.com/lists/python/python/71964 Splitting on a regex w/o consuming delimite
    # and http://stackoverflow.com/questions/3845423/remove-empty-strings-from-a-list-of-strings; the filter removes the blank one from the beginning of '  print ",NCCO".join(str.split('NCCO')).split(',')   '
    # WOW, I WISH I HAD FOUND THE ABOVE COMBINATION OF SOLUTION EARLIER BECAUSE I HAD OTHER PLACES I COULD HAVE USED THESE, ESPECIALLY 'FILTER'

    return the_records_list

def parse_source_organism_info(records_list):
    '''
    returns the list of source organims corresponding to the provided
    genbank-formatted records
    '''
    source_organism_info_list = []
    for record in records_list:
        source_organism_list = record.split('ORGANISM')
        source_organism_info = source_organism_list[1].splitlines()[0]
        extracted_source_organism  = source_organism_info.strip()
        source_organism_info_list.append(extracted_source_organism)
    return source_organism_info_list


def containsAll(str, set):
    '''
    checks if a string contains all characters in a set.

    from
    http://code.activestate.com/recipes/65441-checking-whether-a-string-contains-a-set-of-chars/
    '''
    for c in set:
        if c not in str: return 0;
    return True;

def add_source_data_to_file(provided_fasta_file, id_and_source_dict, fixed_output_file):
    '''
    This function takes a fasta file and uses the provided dictionary of ids as
    keys and the source organism as values to add the source organism info
    between brackets to the description line for each entry, if it seems not
    to be there already.
    It saves a file where the entries with the modified description line named
    with the provided name.
    This file returns the number of entries fixed.

    '''
    #initialize values
    entries_fixed_count = 0
    input_file_stream = open(provided_fasta_file , "r")
    output_file_stream = open(fixed_output_file , "w")
    for line in input_file_stream:
        line = line.strip () #for better control of ends of lines
        # see if it is a description line
        if line.startswith(">"):
            # See if description line has brackets. If it does,
            # assume they contain the source organism info and just save to
            # output file.
            if containsAll(line, '[]'):
                output_file_stream.write(line + "\n")
            else:
                # Since no brackets, need to add source organism info between
                # brackets. Extract the the GI number from the current
                # description line string.
                gi_number = extract_GI_number(line)
                # Use the GI number as the key to pull out the source organism
                # based on matching GI number.
                source_organism_info = id_and_source_dict[gi_number]
                # Insert source organism info flanked by brackets to end of
                # description line and save modified line to the output file
                # being produced.
                output_file_stream.write(
                    line + " ["+ source_organism_info + "]" + "\n")
                entries_fixed_count += 1
        #if not a description line just write text to output file
        else:
            output_file_stream.write(line + "\n")

    #Completed. So close input and output file streams.
    input_file_stream.close()
    output_file_stream.close()
    return entries_fixed_count

###--------------------------END OF HELPER FUNCTIONS---------------------------###










###-----------------Actual Main function of script---------------------------###
###----------------------GET FILE AND PREPARE TO PARSE-----------------------###
def main():
    #file to be provided as a argument when call program.
    #argparser from http://docs.python.org/2/library/argparse.html#module-argparse and http://docs.python.org/2/howto/argparse.html#id1
    parser = argparse.ArgumentParser(
        prog='add_source_organism_info_to_FASTA.py',description=
        "add_source_organism_info_to_FASTA.py adds to each description line the source organism information in between brackets if it appears not present based on lack of brackets. Written by Wayne Decatur --> Fomightez @ Github or Twitter."
        )
    parser.add_argument("InputFile", help="file of sequences in FASTA format. REQUIRED.")
    #I would also like trigger help to display if no arguments provided because need at least input file
    if len(sys.argv)==1:    #from http://stackoverflow.com/questions/4042452/display-help-message-with-python-argparse-when-script-is-called-without-any-argu
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()

    if os.path.isfile(args.InputFile):
        #root_path = path_to_folder_with_file # LEFT HERE FOR USE IN DEBUGGING
        fasta_file = args.InputFile
        logging.debug(fasta_file)
        # APPROACH OUTLINE:
        # Using unique identifiers extracted for FASTA-formatted sequences, get
        # corresponding Genbank entries, parse out source organism info from
        # each as a list, and use `dict(zip(x, y))` , see
        # http://stackoverflow.com/questions/15183084/how-to-create-a-dictionary-using-two-lists-in-python ,
        # to make a dictionary with the two lists that can then be employed
        # to pull the source organism info back up later.
        # Then go through each line of fasta file and
        # when a description line is encountered, check if it already has brackets.
        # If it does, add that line to the output file being produced and move
        # to next line. All lines that aren't description lines (i.e. don't
        # begin with '>') simply add to the output file being produced.
        # For those lines that don't have brackets already present, extract the
        # the GI number from the description line string. Use that as the key to
        # pull out the source organism based on matching GI number and append
        # source organism info flanked by brackets to end of descption line and
        # modified line to the output file being produced.
        #----------------------------------------------------------------------

        # extact unique identifier from FASTA-formatted sequences to use to get
        # corresponding
        sys.stderr.write("Reading in your FASTA file...")
        uids_list, nucleic_or_protein = make_list_of_unique_identifiers (fasta_file)

        # determine database to use when getting GENBANK-formatted record
        if nucleic_or_protein == "nucleic":
            db_to_use = "nuccore"
        else:
            db_to_use = "protein"

        # Get the GENBANK-formatted entries. Want Genbank-formatted entries
        # because easy to parse `SOURCE` from that format because it is a specified
        # field in the record.

        # TAKES TWO STEPS TO ACTUALLY GET THE ENTRIES
        # STEP 1: ONCE I HAVE UNIQUE INDENTIFIERS...
        # In preparation for fetching via EFetch follow recommended practice &
        # upload the list of uids needed to the NCBI Entrez history server via EPost
        sys.stderr.write("\nSending NCBI identifiers to fetch records...")
        #EPost step
        Entrez.email = User_Email     # Always tell NCBI who you are
        webenv, query_key = EPost_uids(uids_list, db_to_use)

        # STEP 2: Fetch the mRNA sequences in GENBANK format via EFetch,
        # following recommneded practice and using the NCBI Entrez history
        # server. The EPost step just prior set the stage for this.
        #EFetch step
        records_list = EFetch_sequences_GB_FORMAT(
            len(uids_list), webenv, query_key, db_to_use)

        # parse out source organism info from from the records list
        source_organism_info_list = parse_source_organism_info(records_list)

        # make a dictionary with the two lists that can then be employed
        # to easily recall the source organism info later.
        # key for each entry will be the id and the value will be the source organism
        # see
        # http://stackoverflow.com/questions/15183084/how-to-create-a-dictionary-using-two-lists-in-python
        id_and_source_dict = dict(zip(uids_list, source_organism_info_list))


        # prepare for saving the fixed data to file
        output_file = generate_output_file_name(fasta_file)

        # Now go through the provided file, fixing description lines that need
        # source organism added in brackets, and save a file with the fixed
        # FASTA-formatted sequences.
        entries_fixed_count = add_source_data_to_file(
            fasta_file, id_and_source_dict, output_file)

        #FOR DEBUGGING
        logging.debug(output_file)

        #give user some stats and feeback
        sys.stderr.write("\nConcluded.\n")
        sys.stderr.write(str(entries_fixed_count)+" entries had the source organism added.\n")
        sys.stderr.write("The new file '"+output_file+"' has been created in same directory as the input file.\n\n")


    else:
        sys.stderr.write("SORRY. " + args.InputFile + " IS NOT RECOGNIZED AS A FILE.\n\n")
        parser.print_help()
        sys.exit(1)


# Did program set-up the way it is below so that I could import this into
# nosetests. See `test_split_SOLiD_reads_usingNOSE_ACTUAL.py` for a reminder.
if __name__ == "__main__":
    main()

#*******************************************************************************
#*******************************************************************************

