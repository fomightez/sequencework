#!/usr/bin/env python
# check_for_omega_intron.py
__author__ = "Wayne Decatur" #fomightez on GitHub
__license__ = "MIT"
__version__ = "0.1.0"


# check_for_omega_intron.py by Wayne Decatur
# ver 0.1
#
#*******************************************************************************
# Verified compatible with both Python 2.7 and Python 3.7; written initially in 
# Python 3. 
#
#
# PURPOSE: Takes a sequence file (FASTA format) containing a fungal mitochondrial 
# large ribosomal subunit RNA and checks if it contains what seems to be a match 
# to the omega intron of S. cerevisiae. Because it uses S. cerevisiae S288C 
# mitochondrial large ribosomal subunit RNA for this effort, I say it is for 
# fungal mitochondrial large ribosomal subunit RNA but depending on 
# conservation, it may be limited to budding yeast or sub-classes of budding 
# yeast. It should help determine that by using with the Shen et al 2018 data 
# for genomes of 332 budding yeast.
#
# The code for this check was originally built in as part of
# `fix_lsu_rRNA_annotation_in_gff_resulting_from_mfannot.py`, but it also 
# could simply be applied when looking at the sequences when you have a 
# previously identified fungal large ribosomal subunit RNA gene or when you are 
# in the process of trying to find them from sequence alone and might want to 
# check intron-harboring status (the latter now being done in 
# `find_mito_fungal_lsu_rRNA_and_check_for_omega_intron.py`). For those reasons 
# I moved it out to its own script so that other scripts can use it by importing 
# it.
# 
# Script requires BLAST+ to be already intalled in environment. An option is
# available at https://github.com/fomightez/blast-binder; go there and click
# `launch binder`, upload this script, and then use the terminal or notebook to
# execute the script.
# 
#
#
#
# CLOSELY RELATED: 
# Other scripts import the main function of this one as part of their
# processing. Current list:
# - fix_lsu_rRNA_annotation_in_gff_resulting_from_mfannot.py
# -find_mito_fungal_lsu_rRNA_and_check_for_omega_intron.py
#
# For impetus, see 
# `for developing find_mito_fungal_lsu_rRNA_and_check_for_omega_intron script.md`
#
#
#
#
# Dependencies beyond the mostly standard libraries/modules:
# sh
# ** NEEDS BLAST+ ALREADY INSTALLED IN ENVIRONMENT.
#
#
#
# VERSION HISTORY:
# v.0.1. basic working version

#
# To do:
# - 
#
#
#
#
# TO RUN:
# Examples,
# Enter on the command line of your terminal, the line
#-----------------------------------
# python check_for_omega_intron.py mito_seq.fa
#-----------------------------------
# Issue `check_for_omega_intron.py -h` for 
# details.
# 
#
#
#
# To use this after importing/pasting or loading into a cell in a Jupyter 
# notebook, specify at least the file of annotations:
# from check_for_omega_intron import check_for_omega_intron
# check_for_omega_intron("mito_seq.fa")
#
# 
#
'''
CURRENT ACTUAL CODE FOR RUNNING/TESTING IN A NOTEBOOK WHEN IMPORTED/LOADED OR 
PASTED IN ANOTHER CELL:
from check_for_omega_intron import check_for_omega_intron
check_for_omega_intron("mito_seq.fa")
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












# Downloaded from https://www.yeastgenome.org/locus/S000007288 using pull down
# to get 'genomic' (S288C_Q0158_21S_RRNA_genomic.fsa) and then to get 'coding' 
# (S288C_Q0158_21S_RRNA_coding.fsa), the latter being without the omega intron
# harboring I-SceI .           This is S. cereivise S288C rnl.
# 4439 bp
cer_rnl  = '''>21S_RRNA Q0158 SGDID:S000007288, chrMito:58009..62447
GTAAAAAGTAGAATAATAGATTTGAAATATTTATTATATAGATTTAAAGAGATAATCATG
GAGTATAATAATTAAATTTAATAAATTTAATATAACTATTAATAGAATTAGGTTACTAAT
AAATTAATAACAATTAATTTTAAAACCTAAAGGTAAACCTTTATATTAATAATGTTATTT
TTTATTATTTTTATAATAAGAATAATTATTAATAATAATAAACTAAGTGAACTGAAACAT
CTAAGTAACTTAAGGATAAGAAATCAACAGAGATATTATGAGTATTGGTGAGAGAAAATA
ATAAAGGTCTAATAAGTATTATGTGAAAAAAATGTAAGAAAATAGGATAACAAATTCTAA
GACTAAATACTATTAATAAGTATAGTAAGTACCGTAAGGGAAAGTATGAAAATGATTATT
TTATAAGCAATCATGAATATATTATATTATATTAATGATGTACCTTTTGTATAATGGGTC
AGCAAGTAATTAATATTAGTAAAACAATAAGTTATAAATAAATAGAATAATATATATATA
TAAAAAAATATATTAAAATATTTAATTAATATTAATTGACCCGAAAGCAAACGATCTAAC
TATGATAAGATGGATAAACGATCGAACAGGTTGATGTTGCAATATCATCTGATTAATTGT
GGTTAGTAGTGAAAGACAAATCTGGTTTGCAGATAGCTGGTTTTCTATGAAATATATGTA
AGTATAGCCTTTATAAATAATAATTATTATATAATATTATATTAATATTATATAAAGAAT
GGTACAGCAATTAATATATATTAGGGAACTATTAAAGTTTTATTAATAATATTAAATCTC
GAAATATTTAATTATATATAATAAAGAGTCAGATTATGTGCGATAAGGTAAATAATCTAA
AGGGAAACAGCCCAGATTAAGATATAAAGTTCCTAATAAATAATAAGTGAAATAAATATT
AAAATATTATAATATAATCAGTTAATGGGTTTGACAATAACCATTTTTTAATGAACATGT
AACAATGCACTGATTTATAATAAATAAAAAAAAATAATATTTAAAATCAAATATATATAT
ATTTGTTAATAGATAATATACGGATCTTAATAATAAGAATTATTTAATTCCTAATATGGA
ATATTATATTTTTATAATAAAAATATAAATACTGAATATCTAAATATTATTATTACTTTT
TTTTTAATAATAATAATATGGTAATAGAACATTTAATGATAATATATATTAGTTATTAAT
TAATATATGTATTAATTAAATAGAGAATGCTGACATGAGTAACGAAAAAAAGGTATAAAC
CTTTTCACCTAAAACATAAGGTTTAACTATAAAAGTACGGCCCCTAATTAAATTAATAAG
AATATAAATATATTTAAGATGGGATAATCTATATTAATAAAAATTTATCTTAAAATATAT
ATATTATTAATAATTATATTAATTAATTAATAATATATATAATTATATTATATATTATAT
ATTTTTTATATAATATAAACTAATAAAGATCAGGAAATAATTAATGTATACCGTAATGTA
GACCGACTCAGGTATGTAAGTAGAGAATATGAAGGTGAATTAGATAATTAAAGGGAAGGA
ACTCGGCAAAGATAGCTCATAAGTTAGTCAATAAAGAGTAATAAGAACAAAGTTGTACAA
CTGTTTACTAAAAACACCGCACTTTGCAGAAACGATAAGTTTAAGTATAAGGTGTGAACT
CTGCTCCATGCTTAATATATAAATAAAATTATTTAACGATAATTTAATTAAATTTAGGTA
AATAGCAGCCTTATTATGAGGGTTATAATGTAGCGAAATTCCTTGGCCTATAATTGAGGT
CCCGCATGAATGACGTAATGATACAACAACTGTCTCCCCTTTAAGCTAAGTGAAATTGAA
ATCGTAGTGAAGATGCTATGTACCTTCAGCAAGACGGAAAGACCCTATGCAGCTTTACTG
TAATTAGATAGATCGAATTATTGTTTATTATATTCAGCATATTAAGTAATCCTATTATTA
GGTAATCGTTTAGATATTAATGAGATACTTATTATAATATAATGATAATTCTAATCTTAT
AAATAATTATTATTATTATTATTAATAATAATAATATGCTTTCAAGCATAGTGATAAAAC
ATATTTATATGATAATCACTTTACTTAATAGATATAATTCTTAAGTAATATATAATATAT
ATTTTATATATATTATATATAATATAAGAGACAATCTCTAATTGGTAGTTTTGATGGGGC
GTCATTATCAGCAAAAGTATCTGAATAAGTCCATAAATAAATATATAAAATTATTGAATA
AAAAAAAAATAATATATATTATATATATTAATTATAAATTGAAATATGTTTATATAAATT
TATATTTATTGAATATATTTTAGTAATAGATAAAAATATGTACAGTAAAATTGTAAGGAA
AACAATAATAACTTTCTCCTCTCTCGGTGGGGGTTCACACCTATTTTTAATAGGTGTGAA
CCCCTCTTCGGGGTTCCGGTTCCCTTTCGGGTCCCGGAACTTAAATAAAAATGGAAAGAA
TTAAATTAATATAATGGTATAACTGTGCGATAATTGTAACACAAACGAGTGAAACAAGTA
CGTAAGTATGGCATAATGAACAAATAACACTGATTGTAAAGGTTATTGATAACGAATAAA
AGTTACGCTAGGGATAATTTACCCCCTTGTCCCATTATATTGAAAAATATAATTATTCAA
TTAATTATTTAATTGAAGTAAATTGGGTGAATTGCTTAGATATCCATATAGATAAAAATA
ATGGACAATAAGCAGCGAAGCTTATAACAACTTTCATATATGTATATATACGGTTATAAG
AACGTTCAACGACTAGATGATGAGTGGAGTTAACAATAATTCATCCACGAGCGCCCAATG
TCGAATAAATAAAATATTAAATAAATATCAAAGGATATATAAAGATTTTTAATAAATCAA
AAAATAAAATAAAATGAAAAATATTAAAAAAAATCAAGTAATAAATTTAGGACCTAATTC
TAAATTATTAAAAGAATATAAATCACAATTAATTGAATTAAATATTGAACAATTTGAAGC
AGGTATTGGTTTAATTTTAGGAGATGCTTATATTCGTAGTCGTGATGAAGGTAAACTATA
TTGTATGCAATTTGAGTGAAAAAATAAGGCATACATGGATCATGTATGTTTATTATATGA
TCAATGAGTATTATCACCTCCTCATAAAAAAGAAAGAGTTAATCATTTAGGTAATTTAGT
AATTACCTGAGGAGCTCAAACTTTTAAACATCAAGCTTTTAATAAATTAGCTAACTTATT
TATTGTAAATAATAAAAAACTTATTCCTAATAATTTAGTTGAAAATTATTTAACACCTAT
AAGTTTAGCATATTGATTTATAGATGATGGAGGTAAATGAGATTATAATAAAAATTCTCT
TAATAAAAGTATTGTATTAAATACACAAAGTTTTACTTTTGAAGAAGTAGAATATTTAGT
TAAAGGTTTAAGAAATAAATTTCAATTAAATTGTTATGTTAAAATTAATAAAAATAAACC
AATTATTTATATTGATTCTATAAGTTATTTAATTTTTTATAATTTAATTAAACCTTATTT
AATTCCTCAAATGATATATAAATTACCTAATACTATTTCATCCGAAACTTTTTTAAAATA
ATATTCTTATTTTTATTTTATGATATATTTCATAAATATTTATTTATATTAAATTTTATT
TGATAATGATATAGTCTGAACAATATAGTAATATATTGAAGTAATTATTTAAATGTAATT
ACGATAACAAAAAATTTGAACAGGGTAATATAGCGAAAGAGTAGATATTGTAAGCTATGT
TTGCCACCTCGATGTCGACTCAACATTTCCTCTTGGTTGTAAAAGCTAAGAAGGGTTTGA
CTGTTCGTCAATTAAAATGTTACGTGAGTTGGGTTAAATACGATGTGAATCAGTATGGTT
CCTATCTGCTGAAGGAAATATTATCAAATTAAATCTCATTATTAGTACGCAAGGACCATA
ATGAATCAACCCATGGTGTATCTATTGATAATAATATAATATATTTAATAAAAATAATAC
TTTATTAATATATTATCTATATTAGTTTATATTTTAATTATATATTATCATAGTAGATAA
GCTAAGTTGATAATAAATAAATATTGAATACATATTAAATATGAAGTTGTTTTAATAAGA
TAATTAATCTGATAATTTTATACTAAAATTAATAATTATAGGTTTTATATATTATTTATA
AATAAATATATTATAATAATAATAATTATTATTATTAATAAAAAATATTAATTATAATAT
TAATAAAATACTAATTTATCAGTTATCTATATAATATCTAATCTATTATTCTATATACT
'''

# VERSION WITHOUT INTRON, I.E. 'coding'. I added 'coding' to description line
# 3296 bps
cer_rnl_coding  = '''>21S_RRNAcoding Q0158 SGDID:S000007288
GTAAAAAGTAGAATAATAGATTTGAAATATTTATTATATAGATTTAAAGAGATAATCATG
GAGTATAATAATTAAATTTAATAAATTTAATATAACTATTAATAGAATTAGGTTACTAAT
AAATTAATAACAATTAATTTTAAAACCTAAAGGTAAACCTTTATATTAATAATGTTATTT
TTTATTATTTTTATAATAAGAATAATTATTAATAATAATAAACTAAGTGAACTGAAACAT
CTAAGTAACTTAAGGATAAGAAATCAACAGAGATATTATGAGTATTGGTGAGAGAAAATA
ATAAAGGTCTAATAAGTATTATGTGAAAAAAATGTAAGAAAATAGGATAACAAATTCTAA
GACTAAATACTATTAATAAGTATAGTAAGTACCGTAAGGGAAAGTATGAAAATGATTATT
TTATAAGCAATCATGAATATATTATATTATATTAATGATGTACCTTTTGTATAATGGGTC
AGCAAGTAATTAATATTAGTAAAACAATAAGTTATAAATAAATAGAATAATATATATATA
TAAAAAAATATATTAAAATATTTAATTAATATTAATTGACCCGAAAGCAAACGATCTAAC
TATGATAAGATGGATAAACGATCGAACAGGTTGATGTTGCAATATCATCTGATTAATTGT
GGTTAGTAGTGAAAGACAAATCTGGTTTGCAGATAGCTGGTTTTCTATGAAATATATGTA
AGTATAGCCTTTATAAATAATAATTATTATATAATATTATATTAATATTATATAAAGAAT
GGTACAGCAATTAATATATATTAGGGAACTATTAAAGTTTTATTAATAATATTAAATCTC
GAAATATTTAATTATATATAATAAAGAGTCAGATTATGTGCGATAAGGTAAATAATCTAA
AGGGAAACAGCCCAGATTAAGATATAAAGTTCCTAATAAATAATAAGTGAAATAAATATT
AAAATATTATAATATAATCAGTTAATGGGTTTGACAATAACCATTTTTTAATGAACATGT
AACAATGCACTGATTTATAATAAATAAAAAAAAATAATATTTAAAATCAAATATATATAT
ATTTGTTAATAGATAATATACGGATCTTAATAATAAGAATTATTTAATTCCTAATATGGA
ATATTATATTTTTATAATAAAAATATAAATACTGAATATCTAAATATTATTATTACTTTT
TTTTTAATAATAATAATATGGTAATAGAACATTTAATGATAATATATATTAGTTATTAAT
TAATATATGTATTAATTAAATAGAGAATGCTGACATGAGTAACGAAAAAAAGGTATAAAC
CTTTTCACCTAAAACATAAGGTTTAACTATAAAAGTACGGCCCCTAATTAAATTAATAAG
AATATAAATATATTTAAGATGGGATAATCTATATTAATAAAAATTTATCTTAAAATATAT
ATATTATTAATAATTATATTAATTAATTAATAATATATATAATTATATTATATATTATAT
ATTTTTTATATAATATAAACTAATAAAGATCAGGAAATAATTAATGTATACCGTAATGTA
GACCGACTCAGGTATGTAAGTAGAGAATATGAAGGTGAATTAGATAATTAAAGGGAAGGA
ACTCGGCAAAGATAGCTCATAAGTTAGTCAATAAAGAGTAATAAGAACAAAGTTGTACAA
CTGTTTACTAAAAACACCGCACTTTGCAGAAACGATAAGTTTAAGTATAAGGTGTGAACT
CTGCTCCATGCTTAATATATAAATAAAATTATTTAACGATAATTTAATTAAATTTAGGTA
AATAGCAGCCTTATTATGAGGGTTATAATGTAGCGAAATTCCTTGGCCTATAATTGAGGT
CCCGCATGAATGACGTAATGATACAACAACTGTCTCCCCTTTAAGCTAAGTGAAATTGAA
ATCGTAGTGAAGATGCTATGTACCTTCAGCAAGACGGAAAGACCCTATGCAGCTTTACTG
TAATTAGATAGATCGAATTATTGTTTATTATATTCAGCATATTAAGTAATCCTATTATTA
GGTAATCGTTTAGATATTAATGAGATACTTATTATAATATAATGATAATTCTAATCTTAT
AAATAATTATTATTATTATTATTAATAATAATAATATGCTTTCAAGCATAGTGATAAAAC
ATATTTATATGATAATCACTTTACTTAATAGATATAATTCTTAAGTAATATATAATATAT
ATTTTATATATATTATATATAATATAAGAGACAATCTCTAATTGGTAGTTTTGATGGGGC
GTCATTATCAGCAAAAGTATCTGAATAAGTCCATAAATAAATATATAAAATTATTGAATA
AAAAAAAAATAATATATATTATATATATTAATTATAAATTGAAATATGTTTATATAAATT
TATATTTATTGAATATATTTTAGTAATAGATAAAAATATGTACAGTAAAATTGTAAGGAA
AACAATAATAACTTTCTCCTCTCTCGGTGGGGGTTCACACCTATTTTTAATAGGTGTGAA
CCCCTCTTCGGGGTTCCGGTTCCCTTTCGGGTCCCGGAACTTAAATAAAAATGGAAAGAA
TTAAATTAATATAATGGTATAACTGTGCGATAATTGTAACACAAACGAGTGAAACAAGTA
CGTAAGTATGGCATAATGAACAAATAACACTGATTGTAAAGGTTATTGATAACGAATAAA
AGTTACGCTAGGGATAACAGGGTAATATAGCGAAAGAGTAGATATTGTAAGCTATGTTTG
CCACCTCGATGTCGACTCAACATTTCCTCTTGGTTGTAAAAGCTAAGAAGGGTTTGACTG
TTCGTCAATTAAAATGTTACGTGAGTTGGGTTAAATACGATGTGAATCAGTATGGTTCCT
ATCTGCTGAAGGAAATATTATCAAATTAAATCTCATTATTAGTACGCAAGGACCATAATG
AATCAACCCATGGTGTATCTATTGATAATAATATAATATATTTAATAAAAATAATACTTT
ATTAATATATTATCTATATTAGTTTATATTTTAATTATATATTATCATAGTAGATAAGCT
AAGTTGATAATAAATAAATATTGAATACATATTAAATATGAAGTTGTTTTAATAAGATAA
TTAATCTGATAATTTTATACTAAAATTAATAATTATAGGTTTTATATATTATTTATAAAT
AAATATATTATAATAATAATAATTATTATTATTAATAAAAAATATTAATTATAATATTAA
TAAAATACTAATTTATCAGTTATCTATATAATATCTAATCTATTATTCTATATACT
'''
def stringFASTA2seq(s):
    '''
    Takes a FASTA file contents that are currently as a Python string and 
    converts it to the sequence string. In other words, it removes the 
    description line and the line breaks and just makes a sequence string.
    '''
    l = s.split("\n")
    return "".join(l[1:])
length_cer_rnl_coding = len(stringFASTA2seq(cer_rnl_coding)) #= 3296<--length 
# in bp of 21S_rRNA without introns; with introns in S288C, it is 4439 bp; 
# S288C omega intron size 1143 bp (Can't just use length of the `cer_rnl_coding` 
# because not read in yet as FASTA and so string includes the description line 
# and line endings.)
near_junction = 2706 #not using last actual base of first exon, 2716, because 
#observed with S. paradoxus CBS432 that because of the limited variety of DNA
# it had 2716 include as the start of a match to the last section of the coding
# sequence, specifically `2716  3296`. 
beyond_insertion_point = 2720 # added a way to check cases where coding 
# (without intron) produces better match than anything to genomic case and 
# checking if greater than a point beyond the junction point (intron insertion 
# point) will make it more stringent. So need such a site defined too.

#*******************************************************************************
#*******************************************************************************
###DO NOT EDIT BELOW HERE - ENTER VALUES ABOVE###

import sys
import os
import subprocess
import pandas as pd


###---------------------------HELPER FUNCTIONS-------------------------------###

    
def determine_omega_presence(seq_file, df = None, bitscore_cutoff = 99):
    '''
    Takes the following:
    `seq_file`, name of sequence file
    `df`, optional dataframe resulting from BLAST query of the coding sequence

    Returns:
    True or false

    Takes a sequence file name, and optionally a dataframe resulting from a 
    BLAST query of the coding sequence of S. cerevisiae S228C 
    Q0158 / 21S_RRNA, and determines whether the omega intron seems present or 
    not.

    The ability to take a dataframe is meant for when using this script imported
    into another where related queries already have taken place and makes no
    sense to repeat.
    '''
    if df is None:
        # Need to do BLAST query with coding sequence of S. cerevisiae S228C 
        # first to make the dataframe other scripts may already have made.
        #
        # Need to store 'cer_rnl' and 'cer_rnl_coding' as files so BLAST can use 
        # them.
        cer_rnl_fn = "cer_rnl_coding.fa"
        with open(cer_rnl_fn, "w") as q_file:
            q_file.write(cer_rnl_coding)
        sys.stderr.write("Checking hits against just coding portion of 21S rRNA"
            "...\n")
        cmd="makeblastdb -in {} -dbtype nucl".format(seq_file)
        try:
            subprocess.run(cmd, shell=True) # based on 
                # https://docs.python.org/3/library/subprocess.html#using-the-subprocess-module
                # and https://stackoverflow.com/a/18739828/8508004
                # NOTE THAT subprocess.run is in Python 3 and not 2, so this 
                # `try` and `except` is to handle Python 2. see 
                # https://stackoverflow.com/a/40590445/8508004
        except AttributeError:
            subprocess.run = py2_run
            subprocess.run(cmd, shell=True) # based on 
            # https://stackoverflow.com/a/40590445/8508004 and 
            # https://github.com/jotyGill/openpyn-nordvpn/issues/82#issue-286134642
        cmd = ('blastn -query {} -db {} -outfmt "6 qseqid sseqid '
            'stitle pident qcovs length mismatch gapopen qstart qend sstart '
            'send qframe sframe frames evalue bitscore qseq '
            'sseq" -task blastn'.format(cer_rnl_fn, seq_file))
        result = subprocess.check_output(cmd, shell=True) # based on 
            # https://stackoverflow.com/a/18739828/8508004
        result = result.decode("utf-8") # so it isn't bytes
        from blast_to_df import blast_to_df
        df = blast_to_df(result)
        df = df[df.bitscore > bitscore_cutoff]
    
    # Determine if it contains omega intron. Doing this by now checking with the 
    # full genomic sequence of cerevisiae 21S rRNA region, exlcuding matches that
    # start above a site near the 5'-exon-into junction, and then examining if 
    # the hits containing that site increase substantially in size. For defining
    # a substantial size increase I am arbitrarily using 20% of 
    # the intron size (1143 bp) or 
    # `length(cer_rnl SEQUENCE)-length(cer_rnl_coding SEQUENCE) * 0.2` ). 
    # Using 5'-end because that part is the largest segment that matches other 
    # omega introns even if homing endonuclease is absent. The part after 
    # I-SceI is much smaller, and thus harder to detect well.
    omega_present = False
    cer_rnl_fn = "cer_rnl.fa"
    with open(cer_rnl_fn, "w") as q_file:
        q_file.write(cer_rnl)
    sys.stderr.write("\nChecking for omega intron presence"
        "...\n")
    cmd="makeblastdb -in {} -dbtype nucl".format(seq_file)
    subprocess.run(cmd, shell=True) # based on 
        # https://docs.python.org/3/library/subprocess.html#using-the-subprocess-module
        # and https://stackoverflow.com/a/18739828/8508004
    cmd = ('blastn -query {} -db {} -outfmt "6 qseqid sseqid '
        'stitle pident qcovs length mismatch gapopen qstart qend sstart '
        'send qframe sframe frames evalue bitscore qseq '
        'sseq" -task blastn'.format(cer_rnl_fn, seq_file))
    full_result = subprocess.check_output(cmd, shell=True) # based on 
        # https://stackoverflow.com/a/18739828/8508004
    full_result = full_result.decode("utf-8") # so it isn't bytes
    from blast_to_df import blast_to_df
    # do some shuffling so the new pickled dataframe with the data from BLAST 
    # query with full genomic sequence of 21S rRNA WITH INTRON doesn't clobber 
    # original
    from sh import mv
    mv("BLAST_pickled_df.pkl", "origBLAST_pickled_df.pkl")
    blast_df_forfull = blast_to_df(full_result)
    # finish shuffle now that have new pickled file and restore original
    mv("BLAST_pickled_df.pkl", "fullBLAST_pickled_df.pkl")
    mv("origBLAST_pickled_df.pkl", "BLAST_pickled_df.pkl")
    sys.stderr.write("\n**DESPITE TEXT ABOVE SAYING `Returning a dataframe "
        "with the information as well`,\nNO DATAFRAME RETURNED YET.**.")
    sys.stderr.write("\nFirst pickled dataframe file remains '{}' and the new "
        "one from\nquery with full, intron-containing cerevisiae 21S rRNA was "
        "named to '{}'.\n".format("BLAST_pickled_df.pkl",
        "fullBLAST_pickled_df.pkl"))
    blast_df_forfull = (
        blast_df_forfull[blast_df_forfull.bitscore > bitscore_cutoff])
    # In order to avoid issues when checking sequences with lots of sequences 
    # that may cause spurious signals, restrict dataframes to contig/sequences, 
    # which are specified by `sseqid`, of best scoring matches.  
    best_scoring_sequence_id = df.loc[df.bitscore.idxmax()].sseqid # best match
    # to the 21S sequence without the intron will specify 21S-containing id. 
    blast_df_forfull = (
        blast_df_forfull[blast_df_forfull.sseqid == best_scoring_sequence_id])
    df = (df[df.sseqid == best_scoring_sequence_id])
    # Before restricting dataframes to sequences spanning in front of the 
    # insertion site, store rawer dfs can use them in case no sequences span the 
    # intron insertion site and need to look at the hits to the intron sequences
    unsubset_blast_df_forfull = blast_df_forfull.copy()
    unsubset_df = df.copy()
    # now to remove matches that start above the `near_junction` site because
    # I want the information for the matches that would span into what would
    # be the 5'-end of the intron, if intron present.
    blast_df_forfull =blast_df_forfull[blast_df_forfull.qstart < near_junction]
    df = df[df.qstart < near_junction]
    # because I chose the segments that begin less than but closest to the 
    # potential 5'-junction to check, it is possible the segments I am checking
    # don't span into the intron actually and I tried to build some testing
    # for that into below.
    # This is also in part why I restricted to the best scoring sequence earlier
    # as well. When I searched against entire genome assembly I was getting hits
    # closer to the junction with scores barely above the cut-off but closer
    # than those in the best scoring contigs/sequences and causing incorrect 
    # calls.
    row_of_interest_for_full = abs(
        blast_df_forfull['qstart'] - near_junction).idxmin()
    row_of_interest_for_coding = abs(df['qstart'] - near_junction).idxmin()
    omega_intron_size = len(stringFASTA2seq(cer_rnl)) - length_cer_rnl_coding
    substantial_increase_in_matches_cutoff = omega_intron_size * 0.2
    # I was seeing an example from the 1011 cerevisiae collection for CAV_4 
    # where it seemed not to have omega intron but it was showing an apparent
    # increase in the size of the so-called hits of interest near the junction 
    # because not comparing same region in the rows of interest.  Instead, the 
    # set of hits (FILTERED by bitscore) from that query (see my file 
    # `Looking at Mexican agave clade of 1011 set relative omega intron.md`) 
    # showed no match spanning the intron insertion site due to the intron 
    # being there, and not matching well enough for long enough for the region 
    # spanning the junction site to produce a hit with a reasonable bitscore. 
    # And the hit of interest from versus the coding sequence was indeed 
    # spanning the intron insertion site,and thus informing that there was a 
    # good match TO NOT HAVING AN INTRON. So this next conditional is to allow 
    # for the situation where the strong match to the sequence without the 
    # intron and none among the filtered set to the genomic S288C seqeunce with 
    # the omega intron.
    full_match_spans_near_junction_val = (
        blast_df_forfull.loc[row_of_interest_for_full].qend > 
        beyond_insertion_point) # probably could just use `> near_junction` but 
    # making it a little larger will make check even more robust/stringent.
    coding_match_spans_near_junction_val = (
        df.loc[row_of_interest_for_coding].qend > beyond_insertion_point)
    if coding_match_spans_near_junction_val and (
        not full_match_spans_near_junction_val):
        omega_present = False
    # This next `elif` is to handle cases like CRI_2 of 1011 cerevisiae where 
    # there are no hits spanning insertion site coding (without intron) sequence
    # but there is a hit spanning the intron insertion site for the genomic
    # sequence with full omega intron. So instead of comparing length of hit for
    # the two rows of interest since not a valid comparison. I just see if
    # the hit spanning the insertion site is at least greater than 18% of the 
    # intron size. (The original problem with CRI_2 may be, at least in part, 
    # due to large chunks of ambiguous nucleotides in the intron. This came up 
    # as an issue when I tried to use `check_for_omega_HEG.py` with this one.)
    elif (not coding_match_spans_near_junction_val) and (
        full_match_spans_near_junction_val):
        if (blast_df_forfull.loc[row_of_interest_for_full].length >= 
            omega_intron_size * 0.18):
            omega_present = True
    # I was seeing some rare cases (YDJ, BSR_4, and GAT of the 1011 collection) 
    # where no well-scoring hits spanning intron insertion site but very well 
    # scoring hits that include intron segment in genomic as part of a hit 
    # spanning the 3'end of the intron. This is to deal with them by looking for 
    # cases where there are no spanning hits but where if you apply a
    # more stringent bitscore cutoff, the size of hits against genomic that 
    # aren't mirrored by coding increase. And the size of the extra hit(s) is 
    # substantial. (This is a variation of my simplisitc previous approach.)
    elif (not coding_match_spans_near_junction_val) and (
        not full_match_spans_near_junction_val):
        blast_df_forfull = unsubset_blast_df_forfull.copy()
        df = unsubset_df.copy()
        blast_df_forfull = (
            blast_df_forfull[blast_df_forfull.bitscore > 300])
        df = (df[df.bitscore > 300])
        full_unique_match_sizes = [x for x in blast_df_forfull.length.to_list(
            ) if x not in df.length.to_list()] # based on 
        # https://stackoverflow.com/a/30040183/8508004
        if full_unique_match_sizes and (max(full_unique_match_sizes) >= 
            substantial_increase_in_matches_cutoff):
            omega_present = True
    # the PRIMARY check for omega intron presence/absence state
    elif (blast_df_forfull.loc[row_of_interest_for_full].length >= 
        df.loc[row_of_interest_for_coding].length + 
        substantial_increase_in_matches_cutoff):
        omega_present = True

    #Previous approach
    # (a variant of this is now been incorporated into being used when no hits
    # span the intron insertion site)
    '''
    # Determine if it contains omega intron. Doing this by now checking with the 
    # full genomic sequence of cerevisiae 21S rRNA region, exlcude lengths of 
    # sequences of the same size matched in both, and see if the length of any 
    # remaining segment increased by a substantial size (arbitrarily use 20% of 
    # the intron size (1143 bp) or `len(cer_rnl)-len(cer_rnl_coding) * 0.2` ). 
    ...
    # now to remove matches seen in both and then see if a substantial length
    # increase has happened
    ...
    full_unique_match_sizes = [x for x in blast_df_forfull.length.to_list(
        ) if x not in df.length.to_list()] # based on 
    # https://stackoverflow.com/a/30040183/8508004
    coding_unique_match_sizes = [x for x in df.length.to_list(
        ) if x not in blast_df_forfull.length.to_list()] 
    if full_unique_match_sizes and coding_unique_match_sizes:
        if max(full_unique_match_sizes) >= max(
            coding_unique_match_sizes) + substantial_increase_in_matches_cutoff:
            omega_present = True

    '''


    return omega_present



def py2_run(*popenargs, **kwargs):
    '''
    backport of `subprocess.run()` to Python 2
    see https://stackoverflow.com/a/40590445/8508004

    Despite being clearly listed at https://stackoverflow.com/a/40590445/8508004 
    and https://github.com/jotyGill/openpyn-nordvpn/issues/82 , the presence of
    `input=None, check=False,` caused issues, and so I removed those
    and corresponding code in the function and then it seemed to work when I
    started a fresh Python console in Jupyter session.
    '''
    process = subprocess.Popen(*popenargs, **kwargs)
    try:
        stdout, stderr = process.communicate(input)
    except:
        process.kill()
        process.wait()
        raise
    retcode = process.poll()
    return retcode, stdout, stderr


###--------------------------END OF HELPER FUNCTIONS-------------------------###
###--------------------------END OF HELPER FUNCTIONS-------------------------###

#*******************************************************************************
###------------------------'main' function of script--------------------------##

def check_for_omega_intron(seq_file, df = None, bitscore_cutoff = 99):
    '''
    Takes the following:
    `seq_file`, name of sequence file
    `df`, optional dataframe resulting from BLAST query of the coding sequence

    Returns:
    True or false


    Main function of script. 
    Takes a name of a sequence file (FASTA format) of a fungal mitochondrial 
    large ribosomal subunit rRNA and checks if it contains what seems to be a 
    match to the omega intron of S. cerevisiae. 

    Has the option to provide a dataframe for results where BLAST query was 
    performed with coding sequence of cerevisiae S288C 21S rRNA large ribosomal 
    RNA of the mitochondria. This is meant to be used in conjunction with 
    another script where that was already done as part of identification, for 
    example the script 
    `fix_lsu_rRNA_annotation_in_gff_resulting_from_mfannot.py` so that the 
    query doesnt need repeating.


    Needs BLAST so make sure already installed from any script calling this one.
    Not repeating that here since will need for most others where this is used
    as a function. Or will be using this as a function imported
    '''

    # Verify gene sequence file present:
    #---------------------------------------------------------------------------
    # Don't continue if not.
    # before BLAST command, make sure `seq_file` is a file present b/c 
    # the error when specified file isn't there seemed unclear because BLAST
    # just returns an error like:
    # 'CalledProcessError: Command 'blastn ....' returned non-zero exit status 2.'
    assert os.path.isfile(seq_file), ("Specified sequence file '{}' "
            "not found.".format(seq_file))
    '''
    try:
        len(seq) > 2 and isinstance(
            seq, str if sys.version_info[0] >= 3 else basestring)
        # the isinstance check comes from comment under 
        # https://stackoverflow.com/a/20612311/8508004 and is meant to make it 
        # work with both Python 2 and 3.
    except NameError:
        github_link = ("https://github.com/fomightez/sequencework/blob/"
            "master/blast-utilities/blast_to_df.py")
        sys.stderr.write("\nSequence string not found."
        ".\n**EXITING !!**.\n".format(github_link))
        sys.exit(1)
    '''

    
    # Verify BLAST installed in environment.
    #--------------------------------------------------------------------
    # 
    #sys.stderr.write("Checking blastn (BLAST+) installed""...\n")
    try:
        cmd="blastn -version"
        result = subprocess.check_output(cmd, shell=True) # based on 
        # https://stackoverflow.com/a/18739828/8508004 
        result = result.decode("utf-8") # so it isn't bytes
        if "blastn:" in result:
            pass
            #sys.stderr.write("Detected '{}'.".format(result.replace("\n","")))
    except CalledProcessError:
        sys.stderr.write("\nblastn not detected. Please install BLAST+ or "
        "run in an\nenvironment launched from "
        "https://github.com/fomightez/blast-binder.\n**EXITING !!**.\n")
        sys.exit(1)
    

    # Get the script containing the function needed for sending 
    # BLAST data to pandas (if not already obtained) & import the function
    #--------------------------------------------------------------------
    #
    # It may have already been brought in because several scripts that call this
    # function already have need for it. So check first.
    file_needed = "blast_to_df.py"
    if not os.path.isfile(file_needed):
        sys.stderr.write("Obtaining script containing function to pass "
        "BLAST results to python"
        "...\n")
        # based on http://amoffat.github.io/sh/
        from sh import curl
        curl("-OL",
            "https://raw.githubusercontent.com/fomightez/sequencework/"
            "master/blast-utilities/blast_to_df.py")
        # verify that worked & ask for it to be done manually if fails
        if not os.path.isfile(file_needed):
            github_link = ("https://github.com/fomightez/sequencework/blob/"
                "master/blast-utilities/blast_to_df.py")
            sys.stderr.write("\n'blast_to_df.py' not found. Please add it to "
            "your current working\ndirectory from {}"
            ".\n**EXITING !!**.\n".format(github_link))
            sys.exit(1)
    from blast_to_df import blast_to_df



    # Check for omega intron in provided sequence
    #---------------------------------------------------------------------
    #
    omega_present = determine_omega_presence(seq_file,df,bitscore_cutoff)



    # Report on omega
    #---------------------------------------------------------------------

    # Report on omega via feedback if called from command line
    if omega_present and __name__ == "__main__":
        sys.stderr.write("\n**Seems omega intron is present. Further "
            "investigation is suggested.**")
    elif (not omega_present) and __name__ == "__main__":
        sys.stderr.write("\n**Simple examination for presence of omega "
            "intron didn't indicate\nobvious presence. Further "
            "investigation is suggested.**")

    if omega_present:
        return True
    else:
        return False




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
    #kwargs['return_df'] = False #probably don't want dataframe returned if 
    # calling script from command line
    check_for_omega_intron(seq_file,**kwargs)
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
    ###-----------------for parsing command line arguments-------------------###
    import argparse
    parser = argparse.ArgumentParser(prog='check_for_omega_intron.py',
        description="check_for_omega_intron.py \
        takes the file name of a sequence containing mitochondrial large \
        ribosomal subunit sequence and returns if omega intron seems present \
        or not.\
        **** Script by Wayne Decatur   \
        (fomightez @ github) ***")

    parser.add_argument("sequence_file", help="Name of file containing the \
        fungal mitochondrial large ribosomal subunit RNA sequence (FASTA \
        format).\
        ", metavar="SEQ_FILE")




    #I would also like trigger help to display if no arguments provided because 
    # need at least one for url
    if len(sys.argv)==1:    #from http://stackoverflow.com/questions/4042452/display-help-message-with-python-argparse-when-script-is-called-without-any-argu
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    seq_file = args.sequence_file


    main()

#*******************************************************************************
###-***********************END MAIN PORTION OF SCRIPT***********************-###
#*******************************************************************************
