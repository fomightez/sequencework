#!/usr/bin/env python
# find_mito_fungal_lsu_rRNA_and_check_for_omega_intron.py
__author__ = "Wayne Decatur" #fomightez on GitHub
__license__ = "MIT"
__version__ = "0.1.0"


# find_mito_fungal_lsu_rRNA_and_check_for_omega_intron.py by Wayne Decatur
# ver 0.1
#
#*******************************************************************************
# Verified compatible with both Python 2.7 and Python 3.7; written initially in 
# Python 3. 
#
#
# PURPOSE: Takes a sequence file of a fungal genome / genome assembly and 
# determines if there is an entry containing the large ribosomal subunit RNA 
# (rnl; 21S rRNA) present. If there is an entry for the large ribosomal subunit 
# RNA (rnl), 
# it notes details and checks if it contains what seems to be a match to the 
# omega intron of S. cerevisiae. Because it uses S. cerevisiae S288C 
# mitochondrial large ribosomal subunit RNA for this effort, I say it is for 
# fungal mitochondrial genomes but depending on conservation, it may be limited 
# to a budding yeast or sub-classes of budding yeast. Use of the script should 
# help determine that by using with the Shen et al 2018 data for genomes of 332 
# budding yeast.
#
# Script requires BLAST+ to be already intalled in environment. An option is
# available at https://github.com/fomightez/blast-binder; go there and click
# `launch binder`, upload this script, and then use the terminal or notebook to
# execute the script.
# 
#
# For impetus, see 
# `for developing find_mito_fungal_lsu_rRNA_and_check_for_omega_intron script.md
# 
# This owes a lot to `fix_lsu_rRNA_annotation_in_gff_resulting_from_mfannot.py`
# which I developed first. The main differences are scope and what are done with 
# the mined details. Here the scope is meant to be all genomic sequences and the 
# details on location of mitochondrial large ribosomal subunit are collected. 
# For `fix_lsu_rRNA_annotation_in_gff_resulting_from_mfannot.py`, the scope is
# a mitochondrial genome sequence that has been annotated by MFannot and the
# details are used to update the annotation to add data on the large ribosomal
# subunit RNA coding region. Because it deals with a sequence and annotation 
# file, the determinations and checks built in are more substantial.
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
# v.0.1. basic working version and fixed to work in Python 2.7 via backport of
# `subprocess.run()` 
#
#
# To do:
# - add passing of adjusted bitscore cutoff when calling `check_for_omega_intron()`
# so that in the cases where the start and end seemed to be on different strands
# and so I adjusted bitscore, that adjusted bitscore used when comparing rows
# generated with BLAST with genomic and coding
#
#
#
#
# TO RUN:
# Examples,
# Enter on the command line of your terminal, the line
#-----------------------------------
# python find_mito_fungal_lsu_rRNA_and_check_for_omega_intron.py mito_seq.fa
#-----------------------------------
# Issue `find_mito_fungal_lsu_rRNA_and_check_for_omega_intron.py -h` for 
# details.
# 
#
#
#
# To use this after importing/pasting or loading into a cell in a Jupyter 
# notebook, specify at least the file of annotations:
# from find_mito_fungal_lsu_rRNA_and_check_for_omega_intron import find_mito_fungal_lsu_rRNA_and_check_for_omega_intron
# find_mito_fungal_lsu_rRNA_and_check_for_omega_intron("mito_seq.fa")
#
# 
#
'''
CURRENT ACTUAL CODE FOR RUNNING/TESTING IN A NOTEBOOK WHEN IMPORTED/LOADED OR 
PASTED IN ANOTHER CELL:
from find_mito_fungal_lsu_rRNA_and_check_for_omega_intron import find_mito_fungal_lsu_rRNA_and_check_for_omega_intron
find_mito_fungal_lsu_rRNA_and_check_for_omega_intron("mito_seq.fa")
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

suffix_for_file_name = "_rnl_info"

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
# 3295 bps
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
length_cer_rnl_coding = 3295 #length in bp of 21S_rRNA without introns; with
# introns in S288C, it is 4439 bp


#*******************************************************************************
#*******************************************************************************
###DO NOT EDIT BELOW HERE - ENTER VALUES ABOVE###

import sys
import os
import subprocess
import pandas as pd




###---------------------------HELPER FUNCTIONS---------------------------------###

def generate_output_file_name_prefix(file_name, suffix_for_file_name):
    '''
    Takes a file name and suffix text as arguments and returns string for the 
    name of the output file. Hence, the generated name is based on the input.

    suffix_for_file_name = "_rnl_info"


    Specific example
    =================
    Calling function with
        ("ogataea_zsoltii_160519-NODE_85_length_19973_cov_251.57_ID_29.fa", "_rnl_info")
    returns
        "ogataea_zsoltii_160519-NODE_85_length_19973_cov_251.57_ID_29_rnl_info"
    '''
    main_part_of_name, file_extension = os.path.splitext(
        file_name) #from http://stackoverflow.com/questions/541390/extracting-extension-from-filename-in-python
    return main_part_of_name + suffix_for_file_name




def extract_hit_strand(frames_data):
    '''
    Takes a string of either "1/1" or "1/-1" and returns
    corresponding strand. I only will be supplying query as positive and so I 
    assume it will always be "1", but best to assert that.

    Returns:
    a string of "+" or "-" to indicate plus/forward or minus/revers strand for 
    the feature. This is the convention used to represent this in GFF3 according 
    to https://useast.ensembl.org/info/website/upload/gff3.html .
    '''
    assert frames_data[0] == "1", ("Strand of query should be a 1. So expect "
        "either\n'1/1' or '1/-1'. Seeing '{}' with {}. **Problem???**".format(
        frames_data, seq_file_name))
    if frames_data.split("/")[1].strip() == "1":
        return "+"
    elif frames_data.split("/")[1].strip() == "-1":
        return "-"
    else:
        sys.stderr.write("\n**ERROR** Cannot extract strand of hit from '{}'"
            ".\n**EXITING !!**.\n".format(frames_data))
        sys.exit(1)





    
def determine_rnl_details(df, seq_file_name, check_omega=True):
    '''
    Takes a blast dataframe from a query of the S. cerevisiae S228C 
    Q0158 / 21S_RRNA against a provided mitochondrial genome and 
    extracts the location and other details of the 21S rRNA gene in the provided 
    sequence. Also because importing another function for handling BLAST query 
    on sequence to check for omega, needs to 'take' sequence file name too.

    Returns start position, and end postion of gene in provided sequence as 
    integers.
    Returns identifying information on specific sequence entry in MULTI-sequence
    FASTA file where the 21S rRNA gene is found. (The idea is that this can be
    used later to compare to other data, such as best matches to other genes, or
    to extract sequence or nearby sequence.)
    Also returns an `+` or `-` to indicate strand.
    And a rough estimate of the quality rated zero to 100.
    And whether the omega intron seems present or not. (Because function was
    inherited from elsewhere has ability to turn off that part but here it is
    part of purpose of script and so don't opt out.)
    '''
    bitscore_cutoff = 99
    df = df[df.bitscore > bitscore_cutoff]
    start_pos = df['sstart'].iloc[df.qstart.idxmin]
    end_pos = df['send'].iloc[df.qend.idxmax]

    # some checks need to be done. Some should trigger errors(?), some warnings, 
    # some trigger fixes(?), and some just mentions. Ideas:
    # - the rows where start_pos and end_pos are harvested should have the same
    # value. Should this trigger an attempt at a fix? It would nice to have an 
    # exampleof it happening to determine how to address. And so maybe for now 
    # just trigger an error so that can find and build in a fix that works for 
    # at least some cases?
    # - start_pos and end_pos shouldn't be separated by more than double the 
    # the size of the rnl coding sequence. Trigger warning? or Error?
    # 

    
    # Get strand info. Note: considered just using `sframe` info but checking
    # `frames` lets me check more things consistent and isn't much overhead.
    strand_for_start_alignment = extract_hit_strand(
        df['frames'].iloc[df.qstart.idxmin])
    strand_for_end_alignment = extract_hit_strand(
        df['frames'].iloc[df.qend.idxmax])
    # I was seeing some where `bitscore > 99` alone was producing cases where
    # a still poorly scoring segment was being assigned as `end` and result
    # seemed to suggest that start and end on different strands and so trying
    # larger cutoff in such instances to see if go away.
    if strand_for_start_alignment != strand_for_end_alignment:
        bitscore_cutoff = 300
        df = df[df.bitscore > bitscore_cutoff]
        start_pos = df['sstart'].iloc[df.qstart.idxmin]
        end_pos = df['send'].iloc[df.qend.idxmax]
        strand_for_start_alignment = extract_hit_strand(
            df['frames'].iloc[df.qstart.idxmin])
        strand_for_end_alignment = extract_hit_strand(
            df['frames'].iloc[df.qend.idxmax])
    # Beyond that attempt to fix, for now triggering an error but would be nice 
    # to fix more if I see a pattern to new cases. (I could imagine this could 
    # happen if rnl spanned start and end of arbitraily placed mito chromosome
    # and then could be reported back so the chromosomes and gff3 could be 
    # adjusted to new start sites, see my circular permuation of mito genomes 
    # stuff, and then rerun.)
    assert strand_for_start_alignment == strand_for_end_alignment, ("Strand "
        "information for fragments containing start and end should be same. "
        "Here seeing '{}' for start and '{}' for end with {}. "
        "**Problem???**".format(
        strand_for_start_alignment, strand_for_end_alignment, seq_file_name))
    strand = strand_for_start_alignment

    # React if start and end coordinates identified are farther apart then
    # double the size of cerevisiae coding. Putting this AFTER checking strand
    # because the strand stuff will handle the cases where this happened because
    # on separate strands.
    if (end_pos - start_pos) > (2 * length_cer_rnl_coding):
        sys.stderr.write("\n**WARNING** The candidate 21S rRNA seems to be "
            "{} bp and this is larger than double the size of S. cerevisiae "
            "21S rRNA. Lots of introns or another issue?"
            ".\n**END OFWARNING**.\n".format(end_pos - start_pos))

    # Collect identifying information about sequence entry where 21S rRNA 
    # located
    id_for_start_alignment = df['sseqid'].iloc[df.qstart.idxmin]
    id_for_end_alignment = df['sseqid'].iloc[df.qend.idxmax]
    # for now triggering an error but would be nice to make so it can be fixed
    # if there is a pattern to most of the cases. (I could imagine this could 
    # happen if rnl spanned start and end of arbitraily placed mito chromosome
    # and then could be reported back so the chromosomes and gff3 could be 
    # adjusted to new start sites, see my circular permuation of mito genomes 
    # stuff, and then rerun.)
    assert id_for_start_alignment == id_for_end_alignment, ("ID "
        "information for fragments containing start and end should be same. "
        "Here seeing '{}' for start and '{}' for end with {}. "
        "**Problem???**".format(
        id_for_start_alignment, id_for_end_alignment,seq_file_name))
    entry_id = id_for_start_alignment

    # calculate a score for quality
    qual_score = 0.0
    for row in df.itertuples():
        qual_score += row.pident * (row.length/float(length_cer_rnl_coding))

    # fix quality if it is over 100% because it is meant to be 100 or less. In
    # theory of the quality is really high and there is some overlap/ repeats it
    # could techically get over 100% if added length actually exceeded the 
    # length of the cerevisiae rnl coding. So just fix by substracting what 
    # amount it is over since redundancy could indicate something so it is 
    # reasonable to lower score for it.
    if qual_score > 100.00:
        qual_score = 100.00 - (qual_score-100.0)

    if check_omega:
        # Determine if it contains omega intron. 
        # I have a separate script with a function to do this because other 
        # scripts have to do this too.
        # Get that script and import the main function.
        file_needed = "check_for_omega_intron.py"
        if not os.path.isfile(file_needed):
            sys.stderr.write("\nObtaining script containing function to check "
                "for omega intron"
                "...\n")
            # based on http://amoffat.github.io/sh/
            from sh import curl
            curl("-OL",
                "https://raw.githubusercontent.com/fomightez/sequencework/"
                "master/omega-presence/check_for_omega_intron.py")
            # verify that worked & ask for it to be done manually if fails
            if not os.path.isfile(file_needed):
                github_link = ("https://github.com/fomightez/sequencework/blob/"
                    "master/omega-presence/")
                sys.stderr.write("\n'check_for_omega_intron.py' not found. "
                    "Please add it to your current working\ndirectory from {}"
                ".\n**EXITING !!**.\n".format(github_link))
                sys.exit(1)
        from check_for_omega_intron import check_for_omega_intron
        # Use the main function to check
        omega_present = check_for_omega_intron(
            seq_file_name, df,bitscore_cutoff)


    
    if check_omega:
        return start_pos,end_pos, entry_id, strand, qual_score, omega_present
    else:
        return start_pos,end_pos, entry_id, strand, qual_score



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


###--------------------------END OF HELPER FUNCTIONS---------------------------###
###--------------------------END OF HELPER FUNCTIONS---------------------------###

#*******************************************************************************
###------------------------'main' function of script---------------------------##

def find_mito_fungal_lsu_rRNA_and_check_for_omega_intron(seq_file_name, 
    return_df=False):
    '''
    Takes the following:
    `seq_file_name`, name of sequence file

    Returns:
    Optionally returns a dataframe of info (meant for when in IPython or Jupyter)

    Generates:
    Saves a tabular text file of the information.


    Main function of script. 
    Takes a sequence file of a fungal genome / genome assembly and determines if 
    there is an entry containing the large ribosomal subunit RNA (rnl; 21S_rRNA) 
    present. If there is an entry for the large ribosomal subunit RNA (rnl), it 
    notes details and checks if it contains what seems to be a match to the 
    omega intron of S. cerevisiae.


    Has options to return some items. These are meant for when calling this
    main function after import into Jupyter/IPython:
    - `return_df` - set to true if you want to return a dataframe with the 
    details on the best hit.


    '''
    # Verify gene sequence file present:
    #---------------------------------------------------------------------------
    # Don't continue if not.
    # before BLAST command, make sure `seq_file_name` is a file present b/c 
    # the error when specified file isn't there seemed unclear because BLASt
    # just returns an error like:
    # 'CalledProcessError: Command 'blastn ....' returned non-zero exit status 2.'
    assert os.path.isfile(seq_file_name), ("Specified sequence file '{}' "
            "not found.".format(seq_file_name))
  

    # Verify BLAST installed in environment.
    #--------------------------------------------------------------------
    # 
    sys.stderr.write("Checking blastn (BLAST+) installed"
        "...\n")
    try:
        cmd="blastn -version"
        result = subprocess.check_output(cmd, shell=True) # based on 
        # https://stackoverflow.com/a/18739828/8508004 
        result = result.decode("utf-8") # so it isn't bytes
        if "blastn:" in result:
            sys.stderr.write("Detected '{}'...\n".format(
                result.replace("\n","")))
    except CalledProcessError:
        sys.stderr.write("\nblastn not detected. Please install BLAST+ or "
        "run in an\nenvironment launched from "
        "https://github.com/fomightez/blast-binder.\n**EXITING !!**.\n")
        sys.exit(1)

    # Get the script containing the function needed for sending BLAST data
    # to pandas and import the function
    #--------------------------------------------------------------------
    #
    # Because I want to build in a way to allow this to be done manually
    # check if already there first and then otherwise try to get.
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


    # Determine location of the mitochondrial large ribosomal subunit in the 
    # provided sequence
    #---------------------------------------------------------------------
    #
    # need to store 'cer_rnl' as a file so BLAST can use it.
    cer_rnl_fn = "cer_rnl_coding.fa"
    with open(cer_rnl_fn, "w") as q_file:
        q_file.write(cer_rnl_coding)
    sys.stderr.write("Locating LSU rRNA (rnl) in the provided "
        "sequence file"
        "...\n")
    cmd="makeblastdb -in {} -dbtype nucl".format(seq_file_name)
    try:
        subprocess.run(cmd, shell=True) # based on 
            # https://docs.python.org/3/library/subprocess.html#using-the-subprocess-module
            # and https://stackoverflow.com/a/18739828/8508004
            # NOTE THAT subprocess.run is in Python 3 and not 2, so this `try` and
            # `except` is to handle Python 2. see 
            # https://stackoverflow.com/a/40590445/8508004
    except AttributeError:
        subprocess.run = py2_run
        subprocess.run(cmd, shell=True) # based on 
        # https://stackoverflow.com/a/40590445/8508004 and 
        # https://github.com/jotyGill/openpyn-nordvpn/issues/82#issue-286134642
    # before BLAST command, make sure seq_file_name is a file present b/c 
    # the error when specified file isn't there seemed unclear because BLAST
    # just returns an error like:
    # 'CalledProcessError: Command 'blastn ....' returned non-zero exit status 2.'
    assert os.path.isfile(seq_file_name), ("Specified sequence file '{}' "
        "not found.".format(seq_file_name))
    cmd = ('blastn -query {} -db {} -outfmt "6 qseqid sseqid '
        'stitle pident qcovs length mismatch gapopen qstart qend sstart '
        'send qframe sframe frames evalue bitscore qseq '
        'sseq" -task blastn'.format(cer_rnl_fn, seq_file_name))
    result = subprocess.check_output(cmd, shell=True) # based on 
        # https://stackoverflow.com/a/18739828/8508004
    result = result.decode("utf-8") # so it isn't bytes
    from blast_to_df import blast_to_df
    blast_df = blast_to_df(result)
    rnl_start, rnl_end, id_, strand, quality, omega_present = (
        determine_rnl_details(blast_df,seq_file_name))
    
    # Report on omega
    if omega_present:
        sys.stderr.write("\n**Seems omega intron is present. Further "
            "investigation is suggested.**")
    else:
        sys.stderr.write("\n**Simple examination for presence of omega "
            "intron didn't indicate\nobvious presence. Further "
            "investigation is suggested.**")

    # Make a dataframe of the details
    #--------------------------------------------------------------------
    details = [(seq_file_name,id_,rnl_start, rnl_end, 
        strand, quality, omega_present)]
    col_labels = (['fasta_source','entry_id', 'start', 'end', 
        'strand','quality','omega_presence'])
    df = pd.DataFrame.from_records(details, columns=col_labels)


    # Reporting and Saving
    #--------------------------------------------------------------------
    main_fn = generate_output_file_name_prefix(seq_file_name, 
        suffix_for_file_name)
    tsv_fn = main_fn + ".tsv"
    pkl_fn = main_fn + ".pkl"
    json_fn = main_fn + ".json"
    df.to_csv(tsv_fn, sep='\t',index = False)
    # feedback
    sys.stderr.write(
        "\n\nDONE.\nDetails stored in tabular text form in the "
        "file '{}'.\n".format(tsv_fn))
    df.to_json(json_fn)
    sys.stderr.write(
        "\n\nDetails stored as json in '{}'.\n".format(json_fn))
    df.to_pickle(pkl_fn)
    sys.stderr.write(
        "\n\nDetails pickled for later use in Pandas in the "
        "file '{}'.\n".format(pkl_fn))

    if return_df:
        sys.stderr.write("Returning a dataframe.\n")
        return df
    else:
        sys.stderr.write("\n**DESPITE TEXT ABOVE SAYING,'Returning a dataframe'"
            " with the information as well`,\nNO DATAFRAME RETURNED.**")




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
    find_mito_fungal_lsu_rRNA_and_check_for_omega_intron(seq_file_name,**kwargs)
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
    parser = argparse.ArgumentParser(prog='find_mito_fungal_lsu_rRNA_and_check_for_omega_intron.py',
        description="find_mito_fungal_lsu_rRNA_and_check_for_omega_intron.py \
        takes a sequence file of a fungal genome / genome assembly and \
        determines if there is an entry containing the large ribosomal subunit \
        RNA (rnl; 21S rRNA) present. If there is an entry for the large \
        ribosomal \
        subunit RNA (rnl), it notes details and checks if it contains what \
        seems to be a match to the omega intron of S. cerevisiae.\
        **** Script by Wayne Decatur   \
        (fomightez @ github) ***")

    parser.add_argument("sequence_file", help="Name of file containing the \
        fungal genomic sequences (FASTA format) harboring mitochondrial large \
        ribosomal subunit RNA somewhere.\
        ", metavar="SEQ_FILE")




    #I would also like trigger help to display if no arguments provided because 
    # need at least one for url
    if len(sys.argv)==1:    #from http://stackoverflow.com/questions/4042452/display-help-message-with-python-argparse-when-script-is-called-without-any-argu
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    seq_file_name = args.sequence_file


    main()

#*******************************************************************************
###-***********************END MAIN PORTION OF SCRIPT***********************-###
#*******************************************************************************
