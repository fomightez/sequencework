#GetmRNAorCDSforProtein.py  by Wayne Decatur
#ver 0.8
#
#
# To GET HELP/MANUAL, enter on command line:
# python GetmRNAorCDSforProtein.py  --help
#
#*************************************************************************
#USES Python 2.7
# Purpose: Takes as input a file of protein sequences in FASTA format
# from NCBI (such as via BLAST) and gets the corresponding mRNA sequences
# (actually cDNAs as the Ts in the sequence indicate) or at least the CDS
# (coding sequence) in FASTA format. Producing a file with those sequences.
# It will give the output in the same order with the exception of redudant entries.
# For redundant entries, the mRNA will be listed as output in the same order as
# the first instanct in the input file.
#
# CAVEAT: This will not work for protein sequences from PDB files like those available
# from PDBsum. The reason is that those files don't have Entrez uids that relate
# to NCBI sequence entries. Therefore cannot get mRNA. Would need to get
# same sequence from NCBI first.
#
#
#
# v.0.8. Added more improvements to handling of mining contigs to make much more
# efficient. Did two things:
# 1) I changed so instead of making a string representing parts proceeding those
# involved with N's, it just substracts the value it would have put from
# what to extract from the assembled sequence in the end. So this saves on calls
# early on and has the added advantage of not having a huge string being built
# and kept in memory.  For example, for protein 355751990 the string would be
# in excess of 9 million. Much easier to just not get that part and then account
# for that in the part extracted later.
# 2) Now if a contig or scaffold was assembled, it stores that as the next part of
# of the cds super assembly coding region is checked. If that involveds the
# same previously assembled contig, it just recycles and uses that instead of
# assembling the sub-components yet again. This also should greatly enhance such
# instances as those protein 355751990, where for example via the old
# method assembling the 'coded by' region of 'join(CM001289.1:9619502..9619604,
# CM001289.1:9619791..9619921,CM001289.1:9620070..9620123,
# CM001289.1:9620240..9620341)' would have involved assembling CM001289.1 for
# each instance separated by the commas, instead of just doing it once.
# NOTE: That even with these improvements in efficiencies, it took me three
# runs on PythonAnywhere of the same script before I got results back when I
# included sequence for protein 355751990 as last in a list of 30 FASTA protein
# sequences. The first two attempts failed with server errors that seemed to be
# from NCBI and not the program. This was confirmed by fact third time worked
# without program being changed. (The time that worked, I also closed the console
# window by going to the the Dashboard on PythonAnywhere. I don't know if it was
# a case of third time being the charm or not having the console window open
# lowered demand. Just noting in case it helps later.)
# Also fixed handling of gaps to handle bothe types listed in GG663364.1  GI:225561288
# related to protein  EEH09737.1  GI:225561457
# http://www.ncbi.nlm.nih.gov/protein/225561457 . It had the 'gap(629)' type I
# knew of and originally designed for. But additionally, it has 'gap(unk100)'
# instances. It has four of them and I assume since all numbered as 100 or the
# overall counting of contig positions it assumes 100 bp for the unknown bp size.
# Also added to function to try coded_by info handling examples like
# http://www.ncbi.nlm.nih.gov/protein/226287726 where complement is first and
# join inside like
'''
            /coded_by="complement(join(KN275970.1:36527..37040,
                         KN275970.1:37237..37550))"
'''
# Under the function 'mine_the_cds_sequence_from_nuccore_seq' I added a try and
# handling of a key error because found gi number 523422945 gives keyerror for
# 'product' with line
'''
seq.id ="gi|"+ genbank_seq_to_mine + ":"+ str(start) + "-" + str(end
    ) + strand + " " + record.description + feature.qualifiers['product'][0]
'''
# so rewrote to handle based on http://stackoverflow.com/questions/2247412/python-best-way-to-check-for-existing-key
# v.0.7.5. Added improvements to handling of mining contigs but making all parts
# before region containing start position of cds for protein as 'N's. This should
# save time and decrease calls to NCBI Entrez server, drastically in some cases.
# v.0.7. Added handling cases where there is not a direct mRNA under
# the ' db="nuccore", LinkName="protein_nuccore_mrna"' part of an Entrez link.
# v.0.7. Added handling cases where there is not a direct mRNA under
# the ' db="nuccore", LinkName="protein_nuccore_mrna"' part of an Entrez link.
# This goes beyond 'GetmmRNAforProtein.py' to get the mRNA or at least a
# CDS corresponding to the protein sequences in FASTA format.
# This version of the program tries the approach here to attempt to
# get the mRNA, and if that fails it tries successively lessening efficient
# methods to get at least the CDS that corresponds to the protein sequence.
# v.0.6. is called 'GetmRNAforProtein.py' and only handles  cases where there is
# a one-to-one relationship between protein sequence and an mRNA under
# the ' db="nuccore", LinkName="protein_nuccore_mrna"' part of an Entrez link
# (see http://www.ncbi.nlm.nih.gov/entrez/query/static/entrezlinks.html).
#
#
#
#
# To do:
# - add diaplay of warning if user tries to run during the non-allowed time a job
# with  more than 25 sequences because need treat as 4 times that given possible
# number of calls to NCBI for difficult sequences to obtain.
#
# - Possibly add a test to see if cds begins with ATG (or GTG?) and ends with
# three nts corresponding to a stop codon. Print a warning for thoee that do not.
# mRNAs obtained by the first method from NCBI should fine and won't necessarily
# correspond to begining with ATG or ending with stop codon, since they are more
# than ORF encoding part often.
#
#
#
# TO RUN:
# For example, enter on the command line, the line
#-----------------------------------
# python GetmRNAorCDSforProtein.py
#-----------------------------------
#
#
# Note that I worked out the basic, first version of this script using my
# IPython notebook entitled 'Using ELink to convert protein sequences to
# mRNA.ipynb'. See that script for further documentation and links related
# to substeps.
#*************************************************************************




##################################
#  USER ADJUSTABLE VALUES        #
##################################
#
User_Email = "YOUR EMAIL HERE" #PUT YOUR E-MAIL HERE or NCBI's SERVER WILL COMPLAIN
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
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
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
    #logging.debug(result)
    for each_record in result:
        #because I found not all protein sequences have entries in the nuccore
        # database I needed to add checking for that first
        if each_record["LinkSetDb"] == []:
            #it seems to return the same number for IdList if it finds nothing
            # and so I can use that it the error information
            record_with_issue = each_record["IdList"][0]
            sys.stderr.write("\n*********************NOTICE**********************")
            sys.stderr.write("\n************************************************")
            sys.stderr.write("\nFASTA entry with GI number " + record_with_issue +
             " failed to return mRNA from the nuccore database.")
            sys.stderr.write("\n************************************************")
            sys.stderr.write("\n************END OF THIS NOTICE******************")
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

    Returns those records as a list.
    '''
    #EFetch
    batch_size = 20
    the_records = ""
    the_records_list = []
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
    the_records_list = filter(None, "!>".join(the_records.split('>')).split('!')) # from http://www.gossamer-threads.com/lists/python/python/71964 Splitting on a regex w/o consuming delimite
    # and http://stackoverflow.com/questions/3845423/remove-empty-strings-from-a-list-of-strings; the filter removes the blank one from the beginning of '  print ",NCCO".join(str.split('NCCO')).split(',')   '
    # WOW, I WISH I HAD FOUND THE ABOVE COMBINATION OF SOLUTION EARLIER BECAUSE I HAD OTHER PLACES I COULD HAVE USED THESE, ESPECIALLY 'FILTER'

    return the_records_list



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


def list_to_string(a_list):
    '''
    Prints contents of a list to a string.
    '''
    for item in a_list:
        print item


def list_for_user_not_obtained_mRNAs(not_obtained_uids):
    '''
    List in a summary at the end all those all protein gi numbers where all
    # three methods failed to get an mRNA sequence.
    '''
    problem_uid_feedback_string = ""
    sys.stderr.write("\nThe following are GI_numbers for which an mRNA was not obtained:\n")
    for uid in not_obtained_uids:
        problem_uid_feedback_string = problem_uid_feedback_string + uid + ', '
    sys.stderr.write(problem_uid_feedback_string[0:-2] + ".\n")




def get_nuccore_uid_for_a_protein_GI_number(GI_number):
    '''
    Takes a GI number for a protein and returns info on nuccore link.
    For example if GI_number is 693903775, I should get output equivalent to
    http://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=protein&db=nuccore&id=693903775
    Also returns a list of any that failed to have a nuc elink result. (Even
    though there's only ever be at the most one in the list, I kept as list
    because easy to deal with methods exists for lists in Python.)
    '''
    nuccore_uid = ""
    difficult_GI_mumbers = []
    #pause program to avoid slamming NCBI Entrez server, see
    # http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc109 and
    # http://www.ncbi.nlm.nih.gov/books/NBK25497/#chapter2.Usage_Guidelines_and_Requiremen
    # about the NCBI Entrez use policy.
    time.sleep(0.4) # this way no chance of more than three calls per second.
    #ELink
    e_link_handle = Entrez.elink(dbfrom="protein", db="nuccore", id=GI_number)
    e_link_result = Entrez.read(e_link_handle)
    e_link_handle.close()
    logging.debug(e_link_result)
    for each_record in e_link_result:
        #because I may stumble upon protein sequences lacking corresponding nuccore
        # enttries I should build in some error handliing
        if each_record["LinkSetDb"] == []:
            #it seems to return the same number for IdList if it finds nothing
            # and so I can use that it the error information
            record_with_issue = each_record["IdList"][0]
            sys.stderr.write("\n*********************ERROR**********************")
            sys.stderr.write("\n************************************************")
            sys.stderr.write("\nFASTA entry with GI number " + record_with_issue +
             " failed to return an entry from the nuccore database.")
            sys.stderr.write("\n************************************************")
            sys.stderr.write("\n************END OF THIS ERROR NOTICE************")
            difficult_GI_mumbers.append(record_with_issue)
        else:
            nuccore_uid = each_record["LinkSetDb"][0]["Link"][-1]["Id"] #protein
            # 357618926 failed when I had 'each_record["LinkSetDb"][0]["Link"][0]["Id"]'
            # earlier because the good link was in last position of the 'Link' list
            # I decided that hopedully in most cases there'd only be one and so last
            # in list would be that one anyway so maybe this is a better selction
            # method to get the uid.
    logging.debug(nuccore_uid)
    return nuccore_uid, difficult_GI_mumbers


def mine_the_cds_sequence_from_nuccore_seq(protein_GI_number, genbank_seq_to_mine):
    '''
    Takes two GI_numbers, one for a protein sequence of interest and one
    for a genbank nucleotide sequence that harbors the coding sequence for that
    protein sequence and returns the coding sequence for the protein in fasta format.
    initially adapted from http://biopython.org/pipermail/biopython/2009-August/005471.html
    '''
    #pause program to avoid slamming NCBI Entrez server, see
    # http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc109 and
    # http://www.ncbi.nlm.nih.gov/books/NBK25497/#chapter2.Usage_Guidelines_and_Requiremen
    # about the NCBI Entrez use policy.
    time.sleep(0.4) # this way no chance of more than three calls per second.
    gb_get_handle = Entrez.efetch(db="nuccore",id=genbank_seq_to_mine,rettype="gb")
    #net_handle = Entrez.efetch(db="nuccore",id="693903281",rettype="gb")
    #e_link_result = Entrez.read(net_handle)  #it seems, this is for xml data so not for gb format
    #e_link_result = net_handle.read() #this was what was used in http://biopython.org/pipermail/biopython/2009-August/005471.html
    # to write and so I thought it would work to get it in Python buffer without
    # needing to write a file.
    record = SeqIO.read(gb_get_handle, "gb") #After worked out above line, I realized
    # SewIO had a read method that needed ot be used to parse so why not use
    # that here and already be ready for next step.
    gb_get_handle.close()


    for feature in record.features:
        if feature.type != "CDS" :
            continue
        gi = feature.qualifiers.get("db_xref", [])
        logging.debug(gi)
        logging.debug([protein_GI_number])
        if gi != ['GI:'+ protein_GI_number]:
            continue
        logging.debug(gi)
        seq = feature.extract(record) #from http://wilke.openwetware.org/Parsing_Genbank_files_with_Biopython.html#Getting_the_actual_coding_sequence
        # 'feature.extract()' handles strandedness automatically
        logging.debug(seq)

        #Need to set ID and description that will be in FASTA to be what I want
        start=feature.location.start.position + 1 # want +1 because this comes from the zero-index python slice but need adjust
        end=feature.location.end.position #help from http://wilke.openwetware.org/Parsing_Genbank_files_with_Biopython.html#Getting_the_actual_coding_sequence
        # getting start and end easily.
        strand = (str(feature.location))[-3:]

        #found gi number 523422945 gives keyerror for 'product' with line
        '''
        seq.id ="gi|"+ genbank_seq_to_mine + ":"+ str(start) + "-" + str(end
            ) + strand + " " + record.description + feature.qualifiers['product'][0]
        '''
        # so rewrote to handle based on http://stackoverflow.com/questions/2247412/python-best-way-to-check-for-existing-key
        # If the 'product' key throws an error, seq.id is just built without
        # using it.
        try:
            seq.id ="gi|"+ genbank_seq_to_mine + ":"+ str(start) + "-" + str(end
            ) + strand + " " + record.description + feature.qualifiers['product'][0]
        except KeyError:
            seq.id ="gi|"+ genbank_seq_to_mine + ":"+ str(start) + "-" + str(end
            ) + strand + " " + record.description
        seq.description = ", corresponds to cds for protein sequence GI_number " + protein_GI_number

        cds_seq_fasta = seq.format("fasta")  #adapted from http://biopython.org/wiki/SeqRecord
        logging.debug("cds sequence")
        logging.debug(cds_seq_fasta[0:100])
        #return cds_seq_fasta


def tsplit(string, delimiters):
    """
    str.split that supports multiple delimiters for splitting on.
    from http://code.activestate.com/recipes/577616-split-strings-w-multiple-separators/

    Examples:
    s = 'thing1,thing2/thing3-thing4'
    tsplit(s, (',', '/', '-'))

    gives:
    ['thing1', 'thing2', 'thing3', 'thing4']


    s = "AP014603.1:1116679..1117326"
    tsplit(s, (':', '..'))
    gives:
    ['AP014603.1', '1116679', '1117326']

    NOTE: In an IPython console on Pythonanywhere I needed to enter by first
    typing 'cpaste'. It worked directly in an Ipython notebook.
    Use '?cpaste' to get info. For some reason 'paste' magic term seemed to not work.

    """
    delimiters = tuple(delimiters)
    stack = [string,]
    for delimiter in delimiters:
        for i, substring in enumerate(stack):
            substack = substring.split(delimiter)
            stack.pop(i)
            for j, _substring in enumerate(substack):
                stack.insert(i+j, _substring)
    return stack

def split_and_grab_gi_num_and_start_and_end_loc(string_to_parse):
    '''
    takes a sting from after coded_by text in a protein genbank format sequence
    and returns gi num of coding sequence and start and end position

    example:
    input:
    AP014603.1:1116679..1117326

    returns:
    AP014603.1,1116679,1117326

    '''
    multi_split_cds_info = tsplit(string_to_parse, (':', '..'))
    gi_number_of_the_nt_seq = multi_split_cds_info[0]
    logging.debug(gi_number_of_the_nt_seq)
    the_start_pos = int(multi_split_cds_info[1]) - 1 #need the minus 1 because Python zero
    # indexed and biopython '.seq' method which I plan to use to get sequence
    # uses zero index, but the Entrez records I retrieved have index of one because
    # first base noted as 1 in sequences and not 0
    logging.debug(the_start_pos)
    the_end_pos= int(multi_split_cds_info[2])
    logging.debug(the_end_pos)
    return gi_number_of_the_nt_seq, the_start_pos, the_end_pos

def fix_cds_info_to_parse_for_rev(a_string):
    '''
    takes a string and removes 'complement(' at start and ')' at end to set up
    text for next step
    '''
    return a_string[11:-1]

def extract_number_of_Ns(the_part):
    '''
    takes something like 'gap(679)' that deals with a gap in a nucleotide sequence
    of a contig such as in entry gi=159122327
    at http://www.ncbi.nlm.nih.gov/nuccore/159122327 (related to protein
    EDP47457.1  GI:159122336) and returns the integer number of Ns in the gap.
    Also handles 'gap(unk100)' style gaps as found in entry gi=225561288,
    # related to protein  EEH09737.1  GI:225561457. All are
    100 so I assume it defaults to that when the size of gap not known. And that
    for sake of listing contig internal positions, it gets treated as 100 bo even though
    size not specifically known yet. Seems so because assuming that, the assembled
    cds for GI:225561457 matched the amino acids.
    '''
    #split on '(' and then take part after that and string through to closing ')'
    parts_list = the_part.split('(')
    info_item_of_parts_list = parts_list[1]
    if info_item_of_parts_list.startswith('unk'):
        #extract number after 'unk'. Examples I found all were 100 but I don't want
        # to assume all use that. I recall original mouse genome had 50 bp.
        return int(info_item_of_parts_list[3:info_item_of_parts_list.index(')')])
    else:
        #extract number in parantheses
        return int(info_item_of_parts_list[:info_item_of_parts_list.index(')')])



def get_nt_seq_using_defined_start_end(record_uid, the_starting_point,
    the_ending_point, want_rev_compl):
    '''
    Takes a gi number (or SEEMS TO WORK WITH accession.version uid as well [ I
     found this out when working on assembling contigs])
    and a start and finish point, plus whether need
    reverse complement, and gets that nucleotide sequence.
    Returning the sequence as a SeqIO object.
    Based on http://biopython.org/pipermail/biopython/2009-August/005471.html
    and
    http://biopython.org/wiki/SeqRecord
    and
    http://www2.warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/genbank/

    See about the different types of uids such as gi number vs. accession.version at
    http://www.ncbi.nlm.nih.gov/Sitemap/sequenceIDs.html

    '''
    #pause program to avoid slamming NCBI Entrez server, see
    # http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc109 and
    # http://www.ncbi.nlm.nih.gov/books/NBK25497/#chapter2.Usage_Guidelines_and_Requiremen
    # about the NCBI Entrez use policy.
    time.sleep(0.4) # this way no chance of more than three calls per second.
    gb_get_handle = Entrez.efetch(db="nuccore",id=record_uid,rettype="gb")
    record = SeqIO.read(gb_get_handle, "gb")
    logging.debug(record)
    gb_get_handle.close()

    if 'contig' in record.annotations:
        #Assemble minimal needed contig. Several rounds of improvements have
        # been added to try increase efficiency, limit memory, and calls to
        # NCBI Entrez server, but this method can still be inefficient and
        # generate A LOT OF calls to NCBI. Use according to policy if you have
        # a lot of entries that need this method, i.e., don't have a 1:1 entry
        # of the mRNA per protein, such as a lot of genomic derived sequences.
        #
        #For now it ends up being RECURSIVE when handling contigs because it calls itself.
        # Looking at it the way I learned recursion in Rice's Principle of Computing course,
        # the base case would be when I do NOT have 'contig' in record.annotations.
        # The way I stumbled into this instance of recrusion I sort of wrote it
        # backwards from way they taught. I would think I would
        # write base case and say if !(contig' in record.annotations):
        # then return sequence otherwise process to get to base case to ultimately
        # return sequence asked for in first place. Since the contig is larger
        # actually then the cds I am trying to assemble, I think this may be a case
        # of generative recursion. It is not that I am feeding smaller and smaller
        # chunks of data.
        # Another way to look at it, is that here to build up the sequence that contains
        # a coding sequence I first have to build up the contig that contains the coding sequence.
        # I say I stumbled into writing a recursive function because as I was
        # implementing handling contigs, I got to a point I needed to get sequence
        # when I had id and end and start and I knew I had a function for that
        # already written and so I plugged that in as the function, and only later
        # did I realize that same function I was calling was the same function I
        # was adding to at the time.
        #
        #step 1 get the info that will be needed to know what is in contig
        # see http://biopython.org/wiki/SeqRecord
        record_contig_raw_info = record.annotations["contig"]
        #step 2 prepare to parse record_contig_raw_info by removing leading 'join('
        # and ending ')'
        record_contig_info = record_contig_raw_info[5:-1]
        #step 3: check if starts with 'complement' and cds_is_on_complement appropriately
        if record_contig_info.startswith("complement"):
            cds_is_on_complement = True
        else:
            cds_is_on_complement = False
        #step 4: split on comma
        multi_contig_parts_to_parse = record_contig_info.split(',')
        logging.debug(multi_contig_parts_to_parse)
        #step 5: step through produced list, if cds_is_on_complement, remove 'complement(' at start and ')' at end, and
        # collect each gi number and position info for 'sub' section of contig
        # and then get each part of the sequence, appending to previous
        #exon_info_list = []
        contig_seq = record.seq[0:1]
        previous_assembled_part_uid = None #initiallizing as nothing so can use
        # to check if last assembled sub part of a sequence comes from same
        # source, so can just re-use if it does.
        for part in multi_contig_parts_to_parse:
            #no need to build contig any more if already through to AT LEAST end position
            if len(contig_seq) > the_ending_point:
                break
            # if the part to parse specifices a gap, make that
            if part.startswith('gap('):
                number_of_Ns = extract_number_of_Ns(part)
                contig_seq = contig_seq + (number_of_Ns * 'N')
            else:
                if cds_is_on_complement:
                    part = fix_cds_info_to_parse_for_rev(part.strip()) #strip method removes leading and trailing white space these will have
                uid_of_nt_seq, sub_start_pos, sub_end_pos = split_and_grab_gi_num_and_start_and_end_loc(part.strip()) #strip method removes leading and trailing white space these will have
                #first determine if the uid of this instance of 'part' is the
                # same as the previous one. If so, just recycle the previously
                # assembled contig. Otherwise assemble the one needed.
                if uid_of_nt_seq == previous_assembled_part_uid:
                    contig_seq = previous_assembled_part
                else:
                    #because I want to be sure contig_seq is type Bio.Seq.Seq and not a string, I am going to use firs instance to define it
                    # QUESTION: wouldn't 'contig_seq = record.seq[0:1]' above have done that? Going to leave though since also lets me set up for adjusting start and end and shouldn't be a huge cost to efficiency unless first one is huge and part needed is very far away from it.
                    if part == multi_contig_parts_to_parse[0]:
                        contig_seq = get_nt_seq_using_defined_start_end(uid_of_nt_seq, sub_start_pos,sub_end_pos, cds_is_on_complement)
                        logging.debug(contig_seq[0:200])
                        # since this signals start it is a nice time to set what will be
                        # adjusted start and end to starting_point and ending_point
                        # These adjustments will be for when the sequence isn't the
                        # fist but it is before the region involved in coding seq.
                        adjusted_start = the_starting_point #initially
                        adjusted_end = the_ending_point  #initially
                    else:
                        #if getting this next sequence will put it in the range of
                        # the sought sequence, then need actual sequence. Otherwise
                        # can skip and lower the the_starting_point and
                        # the_ending_point by that amount. This will be way more
                        # efficient than getting the sequence or even using placeholder
                        # characters since some sequences assembled can easily be
                        # be several million long before getting to coding sequences.
                        #
                        # SEE THE 'ELSE' section below FOR WHY NO' + 1' used on next line - has to due with zero indexing of biopython
                        #
                        # In additional work on efficiency, I added '+ (the_starting_point - adjusted_start)' to line below because
                        # I needed a way account for part not being added between
                        # end of contig and where the sequence would normally be and
                        # '(the_starting_point - adjusted_start)' should reflect that.
                        if (len(contig_seq) + (the_starting_point - adjusted_start
                            ) + (sub_end_pos - sub_start_pos)) >= the_starting_point:
                            part_of_contig_seq = get_nt_seq_using_defined_start_end(uid_of_nt_seq, sub_start_pos,sub_end_pos, cds_is_on_complement)
                            logging.debug(part_of_contig_seq[0:200])
                        else:
                             adjusted_start = adjusted_start - (sub_end_pos - sub_start_pos)
                             adjusted_end = adjusted_end - (sub_end_pos - sub_start_pos) #don't want a 'plus 1' there
                             # that I'd usually put for getting number of elements
                             # between two integers because the way sub_end_pos and sub_start_pos
                             # come back from 'get_nt_seq_using_defined_start_end' function
                             # they are in zero index already.
                             part_of_contig_seq = "" # because do not need to add to
                             # assembled contig, just going to act as if that
                             # sequence wasn't there and shift what is the start
                             # and end accordingly.
                        contig_seq = contig_seq + part_of_contig_seq
             #setting up for possibly recycling if next 'part' in
             # multi_contig_parts_to_parse happens to be the same as previous.
             # for exampe,  parts for coding sequence for protein 355751990 all
             # come from the same contig that involves assembling at least 9 million
             # bp and reusing the one already assembled will be much more
             # efficient.
            previous_assembled_part_uid = uid_of_nt_seq
            previous_assembled_part = contig_seq
        #now since minimal needed contig built, set record.seq to that to set up for extracting portion wanted
        # and adjust the points to factor out the parts not included proceeding those involved
        record.seq = contig_seq
        logging.debug(the_starting_point)
        logging.debug(the_ending_point)
        the_starting_point = adjusted_start
        the_ending_point = adjusted_end
    logging.debug(the_starting_point)
    logging.debug(the_ending_point)
    if want_rev_compl:
        return record.seq[the_starting_point:the_ending_point].reverse_complement()
    else:
        return record.seq[the_starting_point:the_ending_point]


def try_coded_by_info_to_get_cds(protein_involved_gi_num):
    '''
    Takes a GI number for a protein and returns info on cds by using the information
    on the 'coded_by' line of the genbank entry for the protein sequence.
    '''
    #pause program to avoid slamming NCBI Entrez server, see
    # http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc109 and
    # http://www.ncbi.nlm.nih.gov/books/NBK25497/#chapter2.Usage_Guidelines_and_Requiremen
    # about the NCBI Entrez use policy.
    time.sleep(0.4) # this way no chance of more than three calls per second.
    gb_get_handle = Entrez.efetch(db="protein",id=protein_involved_gi_num,rettype="gb")
    #net_handle = Entrez.efetch(db="nuccore",id="693903281",rettype="gb")
    #e_link_txt_result = Entrez.read(net_handle)  #it seems, this is for xml data so not for gb format
    e_link_txt_result = gb_get_handle.read() #this was what was used in http://biopython.org/pipermail/biopython/2009-August/005471.html
    # to write and so I thought it would work to get it in Python buffer without
    # needing to write a file.
    # NOTE THAT CAN ONLY GET RECORD AS TEXT OR READ INTO SeqIO IT SEEMS. CANNOT DO BOTH
    # EVEN IF READ INTO AS TEXT AND TRY TO FEED THAT TO SeqIO
    # protein_record = SeqIO.read(gb_get_handle, "gb") #After worked out above line, I realized
    # SeqIO had a read method that needed ot be used to parse so why not use
    # that here and already be ready for next step.
    gb_get_handle.close()
    #print e_link_txt_result
    #print e_link_txt_result.split('coded_by="')
    coded_by_list = e_link_txt_result.split('coded_by="')
    cds_info_item_of_list = coded_by_list[1]
    cds_info_to_parse = cds_info_item_of_list[:cds_info_item_of_list.index('"')]
    logging.debug(cds_info_to_parse)
    cds_seq = ""
    if cds_info_to_parse.startswith("join"):
        #If 'join', there than need to deal with multiple exons
        ''' examples
        from http://www.ncbi.nlm.nih.gov/protein/159122336
        CDS             1..195
                         /locus_tag="AFUB_093010"
                         /coded_by="join(DS499602.1:44993..45207,
                         DS499602.1:45294..45666)"

        from http://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?noredirect=1&db=protein&val=CAI43140.1

        /coded_by="join(complement(AL034369.1:58366..58379),
                             complement(AL034369.1:56950..57001),
                             complement(AL031177.1:112440..112726),
                             complement(AL031177.1:109978..110238))"

        '''
        sys.stderr.write("\n**********************WARNING***********************")
        sys.stderr.write("\n************************************************")
        sys.stderr.write("\nFASTA entry for protein with GI number " + protein_involved_gi_num + " involves encoding on multiple exons. Extracting this information is in development.")
        sys.stderr.write("\nREVIEWING THE RESULTS FOR THIS ONE IS HIGHLY RECOMMENDED.")
        sys.stderr.write("\n************************************************")
        sys.stderr.write("\n***********END OF THIS WARNING NOTICE************")
        #Step 1: remove "join("  at start and take off ')' at end
        cds_info_to_parse = cds_info_to_parse[5:-1]
        logging.debug("removing 'join(' at start and ')' at end")
        logging.debug(cds_info_to_parse)
        #step 1.5: check if starts with 'complement' and cds_is_on_complement appropriately
        if cds_info_to_parse.startswith("complement"):
            cds_is_on_complement = True
        else:
            cds_is_on_complement = False
        logging.debug(cds_is_on_complement)
        #step 2: split on comma
        multi_exon_parts_to_parse = cds_info_to_parse.split(',')
        logging.debug(multi_exon_parts_to_parse)
        #step 3: step through produced list, if cds_is_on_complement, remove 'complement(' at start and ')' at end, and
        # collect each gi number and position info and then get each part of the sequence, appending to previous
        #exon_info_list = []
        for part in multi_exon_parts_to_parse:
            if cds_is_on_complement:
                part = fix_cds_info_to_parse_for_rev(part.strip()) #strip method removes leading and trailing white space these will have
            uid_of_nt_seq, start_pos, end_pos = split_and_grab_gi_num_and_start_and_end_loc(part.strip()) #strip method removes leading and trailing white space these will have
            part_of_cds_seq = get_nt_seq_using_defined_start_end(uid_of_nt_seq, start_pos,end_pos, cds_is_on_complement)
            logging.debug(part_of_cds_seq[0:200])
            cds_seq = cds_seq + part_of_cds_seq
            #exon_info = [uid_of_nt_seq, start_pos, end_pos]
            #exon_info_list.append(exon_info)
    else:
        if cds_info_to_parse.startswith("complement"):
            cds_is_on_complement = True
            cds_info_to_parse = fix_cds_info_to_parse_for_rev(cds_info_to_parse)
            sys.stderr.write("\n**********************WARNING***********************")
            sys.stderr.write("\n************************************************")
            sys.stderr.write("\nFASTA entry for protein with GI number " + protein_involved_gi_num + " involves encoding on complement. Extracting this information is in development.")
            sys.stderr.write("\nREVIEWING THE RESULTS FOR THIS ONE IS HIGHLY RECOMMENDED.")
            sys.stderr.write("\n************************************************")
            sys.stderr.write("\n***********END OF THIS WARNING NOTICE************")
        else:
            cds_is_on_complement = False
        if cds_info_to_parse.startswith("join"):
            #if join within complement need to handle
            '''
            example:
            from http://www.ncbi.nlm.nih.gov/protein/226287726
            /coded_by="complement(join(KN275970.1:36527..37040,
                         KN275970.1:37237..37550))"
            '''
            sys.stderr.write("\n**********************WARNING***********************")
            sys.stderr.write("\n************************************************")
            sys.stderr.write("\nFASTA entry for protein with GI number " + protein_involved_gi_num + " involves encoding on multiple exons. Extracting this information is in development.")
            sys.stderr.write("\nREVIEWING THE RESULTS FOR THIS ONE IS HIGHLY RECOMMENDED.")
            sys.stderr.write("\n************************************************")
            sys.stderr.write("\n***********END OF THIS WARNING NOTICE************")
            #Step 1: remove "join("  at start and take off ')' at end
            cds_info_to_parse = cds_info_to_parse[5:-1]
            logging.debug("removing 'join(' at start and ')' at end")
            logging.debug(cds_info_to_parse)
            #step 2: split on comma
            multi_exon_parts_to_parse = cds_info_to_parse.split(',')
            logging.debug(multi_exon_parts_to_parse)
            #step 3: step through produced list, and
            # collect each gi number and position info and then get each part of the sequence, appending to previous
            #exon_info_list = []
            for part in multi_exon_parts_to_parse:
                uid_of_nt_seq, start_pos, end_pos = split_and_grab_gi_num_and_start_and_end_loc(part.strip()) #strip method removes leading and trailing white space these will have
                part_of_cds_seq = get_nt_seq_using_defined_start_end(uid_of_nt_seq, start_pos,end_pos, cds_is_on_complement)
                logging.debug(part_of_cds_seq[0:200])
                cds_seq = cds_seq + part_of_cds_seq
                #exon_info = [uid_of_nt_seq, start_pos, end_pos]
                #exon_info_list.append(exon_info)
        else:
            uid_of_nt_seq, start_pos, end_pos = split_and_grab_gi_num_and_start_and_end_loc(cds_info_to_parse)
            #Now should have information needed to fetch sequence and format into a fasta record
            cds_seq = get_nt_seq_using_defined_start_end(uid_of_nt_seq, start_pos,end_pos, cds_is_on_complement)
    logging.debug(cds_seq)
    #Need to set ID and description that will be in FASTA to be what I want
    if cds_is_on_complement:
        strand='(-)'
    else:
        strand='(+)'
    definition_list = e_link_txt_result.split('DEFINITION')
    protein_description = definition_list[1]
    protein_description = (protein_description)[:protein_description.index('ACCESSION')].strip()
    cds_seq.id ="gi|"+ uid_of_nt_seq + ":"+ str(start_pos) + "-" + str(end_pos
        ) + strand + " " + protein_description[:-1] + " " + "mRNA cds"
    cds_seq.description = ", corresponds to cds for protein sequence GI_number " + protein_involved_gi_num
    #turns out '.format' is a SeqRecord function as described here, so need to make a record first
    # before I can use that to format, see http://biopython.org/wiki/SeqRecord
    cds_seq_record = SeqRecord(cds_seq,
                   #IUPAC.ambiguous_dna),
                   id=cds_seq.id, name="cds",
                   description=cds_seq.description)
    cds_seq_fasta = cds_seq_record.format("fasta")  #adapted from http://biopython.org/wiki/SeqRecord
    logging.debug("cds sequence first 200 characters")
    logging.debug(cds_seq_fasta[0:200])
    return cds_seq_fasta


###--------------------------END OF HELPER FUNCTIONS---------------------------###










###-----------------Actual Main function of script---------------------------###
###------------------GET FASTA FILE AND PREPARE TO PARSE---------------------###
#file to be provided as a argument when call program.
#argparser from http://docs.python.org/2/library/argparse.html#module-argparse and http://docs.python.org/2/howto/argparse.html#id1
parser = argparse.ArgumentParser(
    prog='GetmRNAorCDSforProtein.py',description="GetmRNAorCDSforProtein is designed to get mRNA sequence, or at the least a CDS (coding sequence) for each \nentry in a list of protein sequence in FASTA format,\nsuch as from BLAST.\n\n\nWritten by Wayne Decatur --> Fomightez @ Github or Twitter.  \n\nActual example what to enter on command line to run program:\npython GetmRNAorCDSforProtein.py input_file.txt\n ", formatter_class=RawTextHelpFormatter
    )
#learned how to control line breaks in description above from http://stackoverflow.com/questions/3853722/python-argparse-how-to-insert-newline-the-help-text
#DANG THOUGH THE 'RawTextHelpFormatter' setting seems to apply to all the text for argument choices. I don't know yet if that is what really what I wanted.
parser.add_argument("InputFile", help="name of data file containing a list of FASTA protein sequences\nfor which you want mRNA or CDS sequences. REQUIRED.")
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
    #first, this next 'if' condition is so it doesn't try to access NCBI to
    # epost and efetch nothing, which gives an error
    if len(retrieved_mRNA_uids) > 0:
        #Now real 'NEXT', in preparation for fetching via EFetch follow
        # recommended practice and upload the list of uids needed to the NCBI
        # Entrez history server via EPost
        sys.stderr.write("\nSending NCBI the list of mRNA unique ids for the records to fetch...")
        #EPost step
        webenv, query_key = EPost_mRNA_uids(retrieved_mRNA_uids)

        #NEXT...
        #Fetch the mRNA sequences in FASTA format via EFetch, following recommneded
        # practice and using the NCBI Entrez history server.
        sys.stderr.write("\nFetching the mRNA sequence records...")
        #EFetch step
        records_list = EFetch_mRNA_sequences(
            len(retrieved_mRNA_uids), webenv, query_key)
    else:
        #otherwise it won't exist an error will get thrown when try to use it
        # when making the the_records_string below
        records_list = []


    #print the_records #FOR DEBUGGING ONLY

    # Check all is well
    assert len(list_of_GI_numbers) == (len(records_list) + len(not_directly_linkable_uids))



    #NEXT...
    # Put the mRNA sequences into a string in the same order as looked up,
    # except for redundant ones and those not directly linkable to uids. For the
    # redundant ones, the mRNA sequence will be found in the position of the first
    # instance.
    # For those not directly linkable to mRNA uids, put a placeholder that will
    # substituted in the string later when mRNA sequence is determined.
    the_records_string = ""
    idx = 0
    #for index,value in enumerate(list_of_GI_numbers):
    for value in list_of_GI_numbers:
        if str(value) in not_directly_linkable_uids:
            the_records_string = (the_records_string +
                "Placeholder_for_" + str(value) + "\n")
        else:
            the_records_string = the_records_string + records_list[idx]
            idx += 1
    #if not_directly_linkable_uids has anything in it, give a warning that
    # looking these up could take same time.
    if len(not_directly_linkable_uids) > 0:
        sys.stderr.write("\n\nNow attempting to get coding sequences those " + str(len(not_directly_linkable_uids)) + " proteins that don't have direct mRNA correspondences listed in NCBI.")
        sys.stderr.write("\nPLEASE, BE PATIENT. For some protein sequences, for example those derived from genomic sequence, acquiring the coding sequence may require retriving information about millions of basepairs and take some time.")
        sys.stderr.write("\nAs these cases can involve a lot of calls to the NCBI server, definitely only run the program with such cases on off-hours or weekends in accordance with NCBI policy.")
        sys.stderr.write("\nAdditionally, it may require several attempts to get all the sequences without encountering a network issue.")
        sys.stderr.write("\nYou may wish to break the sequences up into chunks of fewer protein sequences or separate out particularly difficult sequences for easier repeated attempts.\n\n")




    #Next go through list of not_directly_linkable_uids and try in turn, the two
    # methods I came up with for those protein sequences lacking a direct unique,
    # solo mRNA entry and rather have to be teased out of some larger sequences.
    # Return those found in a dictionary where the value of the uid is the key.
    gi_nums_resistant_to_even_thorough_methods = []
    difficult_to_get_dictionary = {} #even if ends up empty, need to set as empty so can
    # empty so can have in lines like 'for key in difficult_to_get_dictionary:'
    # and NOT GET AN ERROR THROWN. -> NameError: name 'difficult_to_get_dictionary' is not defined
    for gi_num_id in not_directly_linkable_uids:
        mined_sequence = None #initialize each round so no chance of having last value
        sys.stderr.write("Beginning attempts with "+ str(gi_num_id) + "...")
        logging.debug("step 1 for "+ str(gi_num_id) + "...")
        nuccore_seq_to_mine, failed_to_return_related_nuc_seq = get_nuccore_uid_for_a_protein_GI_number(gi_num_id)
        if len(failed_to_return_related_nuc_seq) == 0:
            logging.debug("nuccore_seq_to_mine is "+ str(nuccore_seq_to_mine))
            logging.debug("step 2 for "+ str(gi_num_id) + "...")
            mined_sequence = mine_the_cds_sequence_from_nuccore_seq(gi_num_id, nuccore_seq_to_mine)
        else:
            logging.debug(("The following failed to return a related nucleotide sequence " +
                str(failed_to_return_related_nuc_seq)))
            logging.debug("thus, cannot try step 2 for "+ str(gi_num_id) + "...")
        #Go to the second option if the 1st option failed because either no
        # nuc sequence was returned by Elink or the 2nd option failed to mine
        # anything from a feature
        if (len(failed_to_return_related_nuc_seq) == 0) or (mined_sequence == None):
            logging.debug("step 3 (THE FINAL STEP) for "+ str(gi_num_id) + "...")
            #try method 2 which relies on 'coded_by' information. IT is the last line
            # of attempts because it several calls to the NCBI Entrez server each time
            # it is run.
            mined_sequence = try_coded_by_info_to_get_cds(gi_num_id)
        logging.debug("mined_sequence is "+ str(mined_sequence))
        difficult_to_get_dictionary[gi_num_id] = mined_sequence

    # Now step through the dictionary containing the uids of those not directly
    # linkable and replace the placeholder text in the_records_string with
    # the fasta entry stroed in the value
    for key in difficult_to_get_dictionary:
        #print key, 'corresponds to', d[key]
        the_records_string = the_records_string.replace ((
            "Placeholder_for_" + key), difficult_to_get_dictionary[key])




    #FINALLY...
    # Save records to an output file and give feeback to user.
    #preparation for saving and for providing feedback.
    TheFileNameMainPart, fileExtension = os.path.splitext(fasta_file_name) #from http://stackoverflow.com/questions/541390/extracting-extension-from-filename-in-python
    if '.' in fasta_file_name:  #I don't know if this is needed with the os.path.splitext method but I had it before so left it
        OutputFileName=  TheFileNameMainPart+"_mRNAseq"+fileExtension
        Issues_FileName=  TheFileNameMainPart+"_failed"+fileExtension
    else:
        OutputFileName= fasta_file_name +"_mRNAseq"
        Issues_FileName=  TheFileNameMainPart+"_failed"

    #right before saving real results, give user some feedback stats and info
    # about issues; the method used here eomploying the re module
    # assumes '>' at beginning of line probably means it is the fasta entry
    # AT first the regular expression was just giving me start of string.
    # See https://docs.python.org/3/howto/regex.html
    # and http://stackoverflow.com/questions/42581/python-re-sub-multiline-caret-match
    # for how I used MULTILINE flag to get it to see all starts of lines .
    # SEe http://stackoverflow.com/questions/1374457/find-out-how-many-times-a-regex-matches-in-a-string-in-python
    # for the basic regular expression use.
    if len(not_directly_linkable_uids) > 0:
        sys.stderr.write(
            "\n"+ str(len(re.findall("^>", the_records_string, flags=re.MULTILINE))
            )+ " mRNA and cds sequence records in FASTA format retrieved from NCBI's Entrez server.")
    else:
        sys.stderr.write(
            "\n"+ str(len(re.findall("^>", the_records_string, flags=re.MULTILINE))
            )+ " mRNA sequence records in FASTA format retrieved from NCBI's Entrez server.")

    #if there were any the program didn't obtain, list those
    if gi_nums_resistant_to_even_thorough_methods != []:
        list_for_user_not_obtained_mRNAs(gi_nums_resistant_to_even_thorough_methods)
        Write2File_orPrint_Lines_Of_Output_List(gi_nums_resistant_to_even_thorough_methods,Issues_FileName)
        sys.stderr.write("\nList of GI_numbers for which no mRNA was obtained written to file '"+ Issues_FileName +"'.\n")


    #Save the good results
    fastafileoutput = open(OutputFileName, "w")
    #the_records_string = '\n'.join(records_list) #this is from when plan was records were a list; from http://www.decalage.info/en/python/print_list
    fastafileoutput.write(the_records_string)
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




