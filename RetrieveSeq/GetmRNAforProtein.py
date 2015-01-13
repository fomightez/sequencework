#GetmRNAforProtein.py  by Wayne Decatur
#ver 0.6
#
#
# To GET HELP/MANUAL, enter on command line:
# python GetmRNAforProtein.py  --help
#
#*************************************************************************
#USES Python 2.7
# Purpose: Takes as input a file of protein sequences in FASTA format
# from NCBI (such as via BLAST) and gets the corresponding mRNA sequences
# (actually cDNAs as the Ts in the sequence indicate) in FASTA format. Producing
# a file with those sequences.
#
# CAVEAT: This will not work for protein sequences from PDB files like those available
# from PDBsum. The reason is that those files don't have Entrez uids that relate
# to Genbank sequence entries. Therefore cannot get mRNA. Would need to get
# same sequence from Genbank first.
#
#
# v.0.6. This works for all cases where there is a one-to-one relationship
# between protein sequence and an mRNA under the ' db="nuccore",
# LinkName="protein_nuccore_mrna"' part of an Entrez link
# (see http://www.ncbi.nlm.nih.gov/entrez/query/static/entrezlinks.html).
# However, it fails with those where the sequence is in a large section of genome
# where a lot of 'cds' are described for example in a genomic sequence. The
# specific 'cds' will be just one of many under the 'protein_nuccore' LinkName.
# It lets you know of its failure with a message and makes a list of those.
# ***See 'GetmmRNAorCDSforProtein.py' to get the mRNA or at least a CDS of protein
# sequences in FASTA format. That program tries the approach here to attempt to
# get the mRNA, and if that fails it tries successively lessening efficient
# methods to get at least the CDS that corresponds to the protein sequence. ***
#
#
#
# TO RUN:
# For example, enter on the command line, the line
#-----------------------------------
# python GetmRNAforProtein.py
#-----------------------------------
#
#
# Note that I worked out this script using my IPython notebook entitled '
# Using ELink to convert protein sequences to mRNA.ipynb'. See that script
# for further documentation and links related to substeps.
#*************************************************************************




##################################
#  USER ADJUSTABLE VALUES        #
##################################
#
User_Email = "wdecatur@yahoo.com" #PUT YOUR E-MAIL HERE or NCBI's SERVER WILL COMPLAIN
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
logging.basicConfig(stream=sys.stderr, level=logging.DEBUG)


#path_to_folder_with_file = "" # LEFT HERE FOR POSSIBLE USE IN DEBUGGING
#FASTA_protein_sequence_records_file = "30test.faa" # LEFT HERE FOR USE IN DEBUGGING



###---------------------------HELPER FUNCTIONS---------------------------------###



def ExtractGI_numbers(a_FASTA_file):
    '''
    This function takes a file containing FASTA entries from Genbank, for example
    collected from BLAST or using Batch Entrez, and mines the GI numbers,
    counting the total number of FASTA records in the process.

    Example, entries in file start with :
    >gi|16580628|emb|CAC82173.1| FSH receptor [Podarcis siculus]
    >gi|37778925|gb|AAO72730.1| follicle-stimulating hormone receptor [Bothrops jararaca]

    Returns total and the list of the GI_numbers:
    2, [1658062, 37778925]

    Want GI number because ELink works with that. I'd have to add another step to get
    to work with accession.version identifier. I originally couldn't get it with
    the indentifier that has the letters like CAC82173.
    See code block #3 of 'Using ELink to convert protein sequences to mRNA.ipynb'
    (publically available at https://www.biostars.org/p/64078/) where
    I did it starting with the accession.version identifier. However, if have full
    Genbank FAST record, just easier to get from that.
    '''
    #initialize values
    the_FASTA_file = open(a_FASTA_file , "r")
    FASTA_entries_tally = 0
    GI_numbers_list = []
    #step through file mining the information
    for line in the_FASTA_file:
        line = line.strip ();#this way I have better control of ends ultimately becaue I know there isn't any unless I add
        if len(line) > 0:
            if line.startswith('>'): # only need those that start with the greater-than symbol as they have the accessions and need to be modified
                FASTA_entries_tally += 1
                sys.stderr.write(".")
                words = line.split()    # In case I need words as I go through
                InfoAndUIDs = line.split("|") #makes it easier to grab designation of type and accession number than words
                FormatOfFasta = (InfoAndUIDs[0])[1:]
                GI_number = InfoAndUIDs[1]
                if GI_number not in GI_numbers_list:  #Do not bother adding if already there because it will save us time later by avoiding repeating lookups
                    GI_numbers_list.append(GI_number) #this will later be submitted Entrez as at http://www.bio-cloud.info/Biopython/en/ch8.html
    #done with user file read in
    the_FASTA_file.close()
    return (FASTA_entries_tally, GI_numbers_list)



def GetmRNA_uids(protein_GI_numbers):
    '''
    returns a list of mRNA sequence identifiers for each GI number in a list of
    GI numbers for protein sequence.
    '''
    the_retrieved_mRNA_uids = []
    difficult_GI_mumbers = []
    #ELink
    handle = Entrez.elink(dbfrom="protein", db="nuccore", LinkName="protein_nuccore_mrna", id=protein_GI_numbers)
    result = Entrez.read(handle)
    handle.close()
    #print result
    for each_record in result:
        #because I found not all protein sequences have entries in the nuccore
        # database I needed to add checking for that first
        if each_record["LinkSetDb"] == []:
            #it seems to return the same number for IdList if it finds nothing
            # and so I can use that it the error information
            record_with_issue = each_record["IdList"][0]
            sys.stderr.write("\n*********************ERROR**********************")
            sys.stderr.write("\n************************************************")
            sys.stderr.write("\nFASTA entry with GI number " + record_with_issue +
             " failed to return mRNA from the nuccore database.")
            sys.stderr.write("\n************************************************")
            sys.stderr.write("\n************END OF THIS ERROR NOTICE************")
            difficult_GI_mumbers.append(record_with_issue)
        else:
            mrna_id = each_record["LinkSetDb"][0]["Link"][0]["Id"]
            the_retrieved_mRNA_uids.append(mrna_id)
    logging.debug(the_retrieved_mRNA_uids)
    logging.debug(difficult_GI_mumbers)
    return the_retrieved_mRNA_uids, difficult_GI_mumbers

def EPost_mRNA_uids(mRNA_uids):
    '''
    Sends a list of uids to the NCBI Entrez history server in preparation for
    fetching step.

    Returns the two items needed to access them off the server (webenv, query_key)
    subsequently.
    '''
    #EPost
    epost_handle = Entrez.epost(db="nuccore", id=",".join(mRNA_uids))
    epost_result = Entrez.read(epost_handle)
    epost_handle.close()

    webenv = epost_result["WebEnv"]
    query_key = epost_result["QueryKey"]
    return webenv, query_key

def EFetch_mRNA_sequences(mRNA_uids_list_length, the_webenv, the_query_key):
    '''
    Fetches the mRNA sequences in FASTA format for a list posted prior to the
    NCBI Entrez history server.

    Returns those records.
    '''
    #EFetch
    batch_size = 20
    the_records = ""
    for start in range(0, mRNA_uids_list_length, batch_size):
        end = min(mRNA_uids_list_length, start + batch_size)
        sys.stderr.write("Fetching records %i thru %i..." % (start + 1, end))
        fetch_handle = Entrez.efetch(db="nuccore",
                                     rettype="fasta", retmode="text",
                                     retstart=start, retmax=batch_size,
                                     webenv=the_webenv,
                                     query_key=the_query_key)
        data = fetch_handle.read()
        fetch_handle.close()
        the_records = the_records + data
    return the_records



def Write2File_orPrint_Lines_Of_Output_List(a_list,OutputFile):
    '''
    Prints contents of a list to lines in file.
    Adapted from http://stackoverflow.com/questions/4675728/redirect-stdout-to-a-file-in-python, see MARCOG's answer and mgold's comments
    Developed as 'Print_Lines_Of_Output_List'
    for SPARTAN08_Fixerv0.7 so I could simply develop with output
    going to stdout and then later easily switch to it going to a file ---
    in the end to send instead from stdout to file just needed to add the
    next two lines nd after close stream and restore process.
    '''
    stdout=sys.stdout
    sys.stdout = open(OutputFile, 'w')
    for the_line in a_list:
        print the_line
    sys.stdout.close()
    sys.stdout = stdout



def list_for_user_not_directly_obtained_mRNAs(not_directly_linked_uids):
    '''
    A lot of protein squences won't have an mRNA under the ' db="nuccore",
    LinkName="protein_nuccore_mrna"' part of an Entrez link. List those in
    summary at end.
    '''
    problem_uid_feedback_string = ""
    sys.stderr.write("\nThe following are GI_numbers for which an mRNA was not obtained due to lack of an mRNA entry at nuccore:\n")
    for uid in not_directly_linked_uids:
        problem_uid_feedback_string = problem_uid_feedback_string + uid + ', '
    sys.stderr.write(problem_uid_feedback_string[0:-2] + ".\n")




###--------------------------END OF HELPER FUNCTIONS---------------------------###










###-----------------Actual Main function of script---------------------------###
###------------------GET FASTA FILE AND PREPARE TO PARSE---------------------###
#file to be provided as a argument when call program.
#argparser from http://docs.python.org/2/library/argparse.html#module-argparse and http://docs.python.org/2/howto/argparse.html#id1
parser = argparse.ArgumentParser(
    prog='GetmRNAforProtein.py',description="GetmRNAforProtein is designed to get mRNA sequence for each \nentry in a list of protein sequence in FASTA format,\nsuch as from BLAST.\n\n\nWritten by Wayne Decatur --> Fomightez @ Github or Twitter.  \n\nActual example what to enter on command line to run program:\npython GetmRNAforProtein.py input_file.txt\n ", formatter_class=RawTextHelpFormatter
    )
#learned how to control line breaks in description above from http://stackoverflow.com/questions/3853722/python-argparse-how-to-insert-newline-the-help-text
#DANG THOUGH THE 'RawTextHelpFormatter' setting seems to apply to all the text for argument choices. I don't know yet if that is what really what I wanted.
parser.add_argument("InputFile", help="name of data file containing a list of FASTA protein sequences\nfor which you want mRNA sequences. REQUIRED.")
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
    number_of_FASTA_entries = 0 #initiate with zero as value of Number of FASTA entries identified


    #Read FASTA entries, keeping track of total, and mine GI numbers
    sys.stderr.write("Reading in your FASTA entries..")
    Entrez.email = User_Email     # Always tell NCBI who you are
    number_of_FASTA_entries, list_of_GI_numbers  = ExtractGI_numbers(fasta_file_name)

    #FOR DEBUGGING
    logging.debug(number_of_FASTA_entries)
    logging.debug(list_of_GI_numbers)

    #give user some feedback stats
    sys.stderr.write(" Concluded. \n"+ str(number_of_FASTA_entries
        )+ " FASTA entries read in with " + str(len (list_of_GI_numbers
        ))+ " different Accession numbers.")



    #NEXT STEP...
    #get id for mRNAs for all the protein sequences via ELink
    sys.stderr.write("\nObtaining information on mRNA sequence records related to entries...")
    #ELink step
    retrieved_mRNA_uids, not_directly_linkable_uids = GetmRNA_uids(list_of_GI_numbers)

    #NEXT...
    #In preparation for fetching via EFetch follow recommended practice and upload
    # the list of uids needed to the NCBI Entrez history server via EPost
    sys.stderr.write("\nSending NCBI the list of mRNA unique ids for the records to fetch...")
    #EPost step
    webenv, query_key = EPost_mRNA_uids(retrieved_mRNA_uids)

    #NEXT...
    #Fetch the mRNA sequences in FASTA format via EFetch, following recommneded
    # practice and using the NCBI Entrez history server.
    sys.stderr.write("\nFetching the mRNA sequence records...")
    #EFetch step
    the_records = EFetch_mRNA_sequences(
        len(retrieved_mRNA_uids), webenv, query_key)






    #print the_records #FOR DEBUGGING ONLY

    #FINALLY
    # Save the mRNA sequences to an output file and give feeback to user.
    #preparation for saving and for providing feedback.
    TheFileNameMainPart, fileExtension = os.path.splitext(fasta_file_name) #from http://stackoverflow.com/questions/541390/extracting-extension-from-filename-in-python
    if '.' in fasta_file_name:  #I don't know if this is needed with the os.path.splitext method but I had it before so left it
        OutputFileName=  TheFileNameMainPart+"_mRNAseq"+fileExtension
        Issues_FileName=  TheFileNameMainPart+"_failed"+fileExtension
    else:
        OutputFileName= fasta_file_name +"_mRNAseq"
        Issues_FileName=  TheFileNameMainPart+"_failed"

    #right before saving real results, give user some feedback stats and info
    # about issues
    sys.stderr.write(
        "\n"+ str(len(retrieved_mRNA_uids)
        )+ " mRNA sequence records in FASTA format retrieved from NCBI's Entrez server.")

    #if there were any the program didn't obtain, list those
    if not_directly_linkable_uids != []:
        list_for_user_not_directly_obtained_mRNAs(not_directly_linkable_uids)
        Write2File_orPrint_Lines_Of_Output_List(not_directly_linkable_uids,Issues_FileName)
        sys.stderr.write("\nList of GI_numbers for which no mRNA was obtained written to file '"+ Issues_FileName +"'.\n")


    #Save the good results
    fastafileoutput = open(OutputFileName, "w")
    fastafileoutput.write(the_records)
    sys.stderr.write("\nResults written to file '"+ OutputFileName +"'.")
    fastafileoutput.close()

    #some additional cleaning up
    #next line is just to clean up so stdout is on next line at end
    sys.stderr.write("\n")





else:
    sys.stderr.write("SORRY. " + args.InputFile + " IS NOT RECOGNIZED AS A FILE.\n\n")
    parser.print_help()
    sys.exit(1)



#*******************************************************************************
#*******************************************************************************




