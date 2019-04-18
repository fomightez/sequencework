#!/usr/bin/env python
# fix_lsu_rRNA_annotation_in_gff_resulting_from_mfannot.py
__author__ = "Wayne Decatur" #fomightez on GitHub
__license__ = "MIT"
__version__ = "0.1.0"


# fix_lsu_rRNA_annotation_in_gff_resulting_from_mfannot.py by Wayne Decatur
# ver 0.1
#
#*******************************************************************************
# Verified compatible with both Python 2.7 and Python 3.6; written initially in 
# Python 3. 
#
#
# PURPOSE: Takes takes a file of annotations (in gff form) of a mitochondrial 
# genome produced by MFAnnot and followed by conversion to gff via a script, the
# corresponding mitochondrial sequence in FASTA format, and 
# determines if there is an entry for the large ribosomal subunit rRNA (rnl) 
# present. If there is no entry for the large ribosomal subunit rRNA (rnl), it 
# adds the annotation to the gff.
# Note: Going from the MFAnnot master files (`.new` extesion) preferably use 
# `mfannot2gff3.pl` from
# https://github.com/yjx1217/LRSDAY/blob/master/scripts/mfannot2gff3.pl . I have 
# also used in the past `mfannot2gff.pl`from 
# https://github.com/kbseah/mitonotate/blob/master/mfannot2gff.pl` but noted it 
# suffered from an apparent one-off error for the RNAs at the start, see 
# https://github.com/kbseah/mitonotate/issues/5. And so I tried the LRSday one 
# (https://github.com/yjx1217/LRSDAY/blob/master/scripts/mfannot2gff3.pl) most
# recently. NOTE that there is a difference in output from these two. 
# mfannot2gff.pl includes notes on the translation table used in the mfannot run 
# on each feature. I have added in checking if the note about the translation 
# table is present in each feature and hande the output to try and match 
# conversion script output.
# Script requires BLAST+ to be allready intalled in environment. An option is
# available at https://github.com/fomightez/blast-binder; go there and click
# `launch binder`, upload this script, and then use the terminal or notebook to
# execute the script.
# 
#
#
#
#
# Developed as a follow-up to 
# `Trying some extracted mitos from 332 in MFannot April 5 2019.ipynb` and 
# general dissatification with MFannot annotation of the large ribosome subunit 
# in budding yeast mitochondrial sequences going back to around August 2018.
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
# - add outputing quality assessment in demo (I think I make it now but don't have anywhere used yet)
# - Mke a separate script that just checks of omega present or not, which would be useful for 
# my extracted cerevisiae mitos where I don't need to fix the gffs at this time 
# (or after I fix the cerevisiae ones to run along with those that didn't need fixing), I should
# move the handling in this script to that one and just import that main 
# function into here like I did with the blast_to_df. Include in notes on that
#https://www.ncbi.nlm.nih.gov/pubmed/6361491 and as noted in `Counting putative 
# promoters among 332 Saccharomycotina.ipynb`. Also these notes just below should be moved 
# there are the associatd markdown file:
'''
yueomyces_sinensis_mito 20145-20664 looks to be the ortholog of the cerevisiae promoter at the start of the 21S rRNA ortholog; however, it seems yueomyces_sinensis maybe(?) lacks the Omega intron because I am seeing the parts of the yueomyces_sinensis be split between best matches covering 58020-58313 (corresponding to 17-295 of the 520 bp yueomyces_sinensis sequence) and 58318-58482 of cerevisiae (the latter part corresponding to 343-520 of 520 bp yueomyces_sinensis sequence).

to help with that FROM SGD:
noncoding_exon  1..2716 chrmt:58009..60724  2000-05-19  2000-05-19
intron  2717..3859  chrmt:60725..61867  2000-05-19  2000-05-19
noncoding_exon  3860..4439  chrmt:61868..62447  2000-05-19  2000-05-19
'''
#
#
#
#
# TO RUN:
# Examples,
# Enter on the command line of your terminal, the line
#-----------------------------------
# python fix_lsu_rRNA_annotation_in_gff_resulting_from_mfannot.py mito_annotations.gff3 mito_seq.fa
#-----------------------------------
# Issue `fix_lsu_rRNA_annotation_in_gff_resulting_from_mfannot.py -h` for 
# details.
# 
#
#
#
# To use this after importing/pasting or loading into a cell in a Jupyter 
# notebook, specify at least the file of annotations:
# from fix_lsu_rRNA_annotation_in_gff_resulting_from_mfannot import fix_lsu_rRNA_annotation_in_gff_resulting_from_mfannot
# fix_lsu_rRNA_annotation_in_gff_resulting_from_mfannot("mito_annotations.gff3", "mito_seq.fa")
#
# 
#
'''
CURRENT ACTUAL CODE FOR RUNNING/TESTING IN A NOTEBOOK WHEN IMPORTED/LOADED OR 
PASTED IN ANOTHER CELL:
from fix_lsu_rRNA_annotation_in_gff_resulting_from_mfannot import fix_lsu_rRNA_annotation_in_gff_resulting_from_mfannot
fix_lsu_rRNA_annotation_in_gff_resulting_from_mfannot("mito_annotations.gff3", "mito_seq.fa")
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

text_to_add_to_altered_file = "_rnlFIXED"

#
#*******************************************************************************
#**********************END USER ADJUSTABLE VARIABLES****************************












# Downloaded from https://www.yeastgenome.org/locus/S000007288 using pull down
# to get 'genomic' (S288C_Q0158_21S_RRNA_genomic.fsa) and then to get 'coding' 
# (S288C_Q0158_21S_RRNA_coding.fsa), the latter being without the omega intron
# harboring I-SceI .           This is S. cereivise S288C rnl.
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
length_cer_rnl_coding = 3295 #length in bp of 21S_rRNA without introns


#*******************************************************************************
#*******************************************************************************
###DO NOT EDIT BELOW HERE - ENTER VALUES ABOVE###

import sys
import os
import subprocess
import pandas as pd




###---------------------------HELPER FUNCTIONS---------------------------------###

def generate_output_file_name(file_name, text_to_add_to_altered_file):
    '''
    Takes a file name as an argument and returns string for the name of the
    output file. The generated name is tagged to indicate difference.


    Specific example
    =================
    Calling function with
        ("mito_annotation.gff3", "_rnlFIXED")
    returns
        "mito_annotation_rnlFIXED.gff3"
    '''
    main_part_of_name, file_extension = os.path.splitext(
        file_name) #from http://stackoverflow.com/questions/541390/extracting-extension-from-filename-in-python
    if '.' in file_name:  #I don't know if this is needed with the os.path.splitext method but I had it before so left it
        return main_part_of_name + text_to_add_to_altered_file  + file_extension
    else:
        return file_name + text_to_add_to_altered_file + ".gff3"

def get_start_end_gff_line(current_line):
    '''
    takes a line in gff3 format and returns `start` and `end` columns

    returns the values as integers

    see https://useast.ensembl.org/info/website/upload/gff3.html
    '''
    parts = current_line.split("\t")
    start = parts[3]
    end = parts[4]
    return int(start),int(end)

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
        "either\n'1/1' or '1/-1'. Seeing '{}'. **Problem???**".format(
        frames_data))
    if frames_data.split("/")[1].strip() == "1":
        return "+"
    elif frames_data.split("/")[1].strip() == "-1":
        return "-"
    else:
        sys.stderr.write("\n**ERROR** Cannot extract strand of hit from '{}'"
            ".\n**EXITING !!**.\n".format(frames_data))
        sys.exit(1)

def determine_if_annotations_include_rnl(file_name):
    '''
    takes a file of annotations (in gff form) of a mitochondrial genome 
    produced by MFAnnot followed by conversion via a script and determines 
    if there is an entry for the large ribosomal subunit rRNA (rnl) present.

    RETURNS True or false.

    Going from the MFAnnot master files (`.new` extesion) preferably use 
    `mfannot2gff3.pl` from 
    https://github.com/yjx1217/LRSDAY/blob/master/scripts/mfannot2gff3.pl . I 
    have also used in the past `mfannot2gff.pl`from 
    https://github.com/kbseah/mitonotate/blob/master/mfannot2gff.pl` but noted
    it suffered from an apparent one-off error for the RNAs at the start
    , see https://github.com/kbseah/mitonotate/issues/5 and so I tried the
    LRSday one 
    (https://github.com/yjx1217/LRSDAY/blob/master/scripts/mfannot2gff3.pl) most
    recently.


    It is able to discern `trnl` from `rnl`.
    Case-insensitive. (Note I didn't have to make that way and then indeed it
    would have been easier to discern `trnL` from `rnl`. However, making it
    insensitive will make it more robust in long run because won't need to rely
    on someone being sure to capitlize `L` at of `trnL`,)
    RETURNS True or false.
    '''
    rnl_present = False
    with open(file_name, 'r') as input:
        for line in input:
            if "rnl" in line.lower() and "trnl" not in line.lower():
                rnl_present = True
                return rnl_present # no point in continuing to read in rest
    return rnl_present

    
def determine_rnl_start_end(df, check_omega=True):
    '''
    Takes a blast dataframe from a query of the S. cerevisiae S228C 
    Q0158 / 21S_RRNA against a provided mitochondrial genome and 
    extracts the location of the rnl gene in the provided sequence.

    Returns start position, and end postion of gene in provided sequence as 
    integers.
    Also returns an `+` or `-` to indicate strand.
    And a rough estimate of the quality rated zero to 100.
    And, optionally, whether the omega intron seems present or not.
    '''
    df = df[df.bitscore > 99]
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
        df['frames'].iloc[df.qstart.idxmin])
    # for now triggering an error but would be nice to make so it can be fixed
    # if there is a pattern to most of the cases. (I could imagine this could 
    # happen if rnl spanned start and end of arbitraily placed mito chromosome
    # and then could be reported back so the chromosomes and gff3 could be 
    # adjusted to new start sites, see my circular permuation of mito genomes 
    # stuff, and then rerun.)
    assert strand_for_start_alignment == strand_for_end_alignment, ("Strand "
        "information for fragments containing start and end should be same. "
        "Here seeing '{}' for start and '{}' for end. **Problem???**".format(
        strand_for_start_alignment, strand_for_end_alignment))
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

    # calculate a score for quality
    qual_score = 0.0
    for row in df.itertuples():
        qual_score += row.pident * (row.length/length_cer_rnl_coding)

    # fix quality if it is over 100% because it is meant to be 100 or less. In
    # theory of the quality is really high and there is some overlap/ repeats it
    # could techically get over 100% if added length actually exceeded the 
    # length of the cerevisiae rnl coding. So just fix by substracting what 
    # amount it is over since redundancy could indicate something so it is 
    # reasonable to lower score for it.
    if qual_score > 100.00:
        qual_score = 100.00 - (qual_score-100.0)

    # Determine if it contains omega intron. Do this by now checking with the 
    # full genomic sequence of cerevisiae 21SrRNA region and see if the matches
    # are less in number
    omega_present = False
    cer_rnl_fn = "cer_rnl.fa"
    with open(cer_rnl_fn, "w") as q_file:
        q_file.write(cer_rnl)
    sys.stderr.write("\nChecking for omega intron presence"
        "...\n")
    cmd="makeblastdb -in {} -dbtype nucl".format(seq_file_name)
    subprocess.run(cmd, shell=True) # based on 
        # https://docs.python.org/3/library/subprocess.html#using-the-subprocess-module
        # and https://stackoverflow.com/a/18739828/8508004
    cmd = ('blastn -query {} -db {} -outfmt "6 qseqid sseqid '
        'stitle pident qcovs length mismatch gapopen qstart qend sstart '
        'send qframe sframe frames evalue bitscore qseq '
        'sseq" -task blastn'.format(cer_rnl_fn, seq_file_name))
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
    sys.stderr.write("\nFirst pickled dataframe file remains '{}' and the new "
        "one from\nquery with full, intron-containing cerevisiae 21S rRNA was "
        "named to '{}'.\n".format("BLAST_pickled_df.pkl","fullBLAST_pickled_df.pkl"))
    blast_df_forfull = blast_df_forfull[blast_df_forfull.bitscore > 99]
    if len(blast_df_forfull) < len(df):
        omega_present = True


    strand = "+"
    qual_score = 50

    if check_omega:
        return start_pos,end_pos, strand, qual_score, omega_present
    else:
        return start_pos,end_pos, strand, qual_score


def make_updated_txt(gff_text_lines, line_to_insert_after,
    line_to_insert_before,new_entry):
    '''
    Takes the following:
    - a list of text lines for the gff3 file 
    - line_to_insert_after, which is an integer
    - line_to_insert_before, which is either an integer or "NOT APPLICABLE"
    - new_entry , which is the text to add

    Returns:
    updated text including the new entry
    '''
    if line_to_insert_before == "NOT APPLICABLE":
        return "\n".join(gff_text_lines + [new_entry])
    return "\n".join(gff_text_lines[:line_to_insert_after + 1] + [new_entry] + 
        gff_text_lines[line_to_insert_before:])


###--------------------------END OF HELPER FUNCTIONS---------------------------###
###--------------------------END OF HELPER FUNCTIONS---------------------------###

#*******************************************************************************
###------------------------'main' function of script---------------------------##

def fix_lsu_rRNA_annotation_in_gff_resulting_from_mfannot(gff_file_name, 
        seq_file_name, return_status=False, return_quality_estimate=False):
    '''
    Main function of script. 
    Takes takes a file of annotations (in gff form) of a mitochondrial 
    genome produced by MFAnnot and followed by conversion to gff via a script,
    the corresponding fungal mitochondrial genome (FASTA format), 
    and determines if there is an entry for the large ribosomal subunit rRNA 
    (rnl) present. If there is no entry for the large ribosomal subunit rRNA 
    (rnl), it adds the annotation to the gff.

    Has options to return some items. These are meant for when calling this
    main function after import into Jupyter/IPython:
    - `return_status` - set to true if you want to return whether `rnl` gene
    found in original or a fix was implemented
    - `return_quality_estimate` - set to true if you want to return quality 
    estimate with zero to 100 score.

    Note: Going from the MFAnnot master files (`.new` extesion) preferably use 
    `mfannot2gff3.pl` from 
    https://github.com/yjx1217/LRSDAY/blob/master/scripts/mfannot2gff3.pl . I 
    have also used in the past `mfannot2gff.pl`from 
    https://github.com/kbseah/mitonotate/blob/master/mfannot2gff.pl` but noted 
    it suffered from an apparent one-off error for the RNAs at the start, see 
    https://github.com/kbseah/mitonotate/issues/5. And so I tried the LRSday one 
    (https://github.com/yjx1217/LRSDAY/blob/master/scripts/mfannot2gff3.pl) most
    recently.
    '''
    # See if annotation present or not:
    #---------------------------------------------------------------------------
    # Don't continue if already present
    if determine_if_annotations_include_rnl(gff_file_name):
        sys.stderr.write("\n***************************************************"
            "******\n"
            "Annotation of large ribosomal subunit by MFannot detected in\n'{}'"
            " and so no further processing necessary."
            "\n********************************************************"
            "*\nDONE.\n\n".format(gff_file_name))
        if return_status:
            status = "present"
        if return_quality_estimate:
            quality = "100.0" #because must be good if mfannot identified it
    else:
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
        sys.stderr.write("Obtaining script containing function to pass "
            "BLAST results to python"
            "...\n")
        # based on http://amoffat.github.io/sh/
        from sh import curl
        curl("-OL",
            "https://raw.githubusercontent.com/fomightez/sequencework/"
            "master/blast-utilities/blast_to_df.py")
        # verify that worked
        file_needed = "blast-utilities/blast_to_df.py"
        try:
            os.path.isfile(file_needed)
        except NameError:
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
        sys.stderr.write("Locating LSU rRNA (rnl) in the provided fungal "
            "mitochondrial genome sequence file"
            "...\n")
        cmd="makeblastdb -in {} -dbtype nucl".format(seq_file_name)
        subprocess.run(cmd, shell=True) # based on 
            # https://docs.python.org/3/library/subprocess.html#using-the-subprocess-module
            # and https://stackoverflow.com/a/18739828/8508004
        # before BLAST command, make sure seq_file_name is a file present b/c 
        # the error when specified file isn't there seemed unclear because BLASt
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
        rnl_start, rnl_end, strand, quality, omega_present = (
            determine_rnl_start_end(blast_df))
        
        # Report on omega
        if omega_present:
            sys.stderr.write("\n**Seems omega intron is present. Further "
                "investigation is suggested.**")
        else:
            sys.stderr.write("\n**Simple examination for presence of omega "
                "intron didn't indicate obvious presence. Further "
                "investigation is suggested.**")


        # Get the file to be fixed in as an object, determination point to add
        # new line and add it.
        #--------------------------------------------------------------------
        # When ran `determine_if_annotations_include_rnl()`, earlier just read 
        # one at a time and didn't store anything. The idea being there that 
        # they might already have annotation and so keep as simple as possible
        # stopping as soon as dtermined.
        # Easy here to read in whole thing so I can insert at proper point
        lines_to_insert_between = (0,0)
        with open(gff_file_name, 'r') as gff_handle:
            gff_text=gff_handle.read()
        gff_text_lines = gff_text.split("\n")
        gff_text_lines = [x for x in gff_text_lines if x] #remove blank lines 
        # so any blanks at end of dile not an issue with iterating.
        # The first line should be a header line indicating gff3 accroding to 
        # https://useast.ensembl.org/info/website/upload/gff3.html and that 
        # also means there is no chance of needing to insert the rnl data before
        # the first line and so except for the possibility of the last line 
        # there is no need to have to handle edge case where there isn't two 
        # lines to insert between.
        assert gff_text_lines[0].startswith("#") , ("The first line of a GFF3 "
            "file must be a comment that\nidentifies the version, see "
            "https://useast.ensembl.org/info/website/upload/gff3.html .")
        line_to_insert_after = 0
        line_to_insert_before = 0
        # determine `line_to_insert_after` and `line_to_insert_before`. Also
        # check if `transl_table=` used or not?
        transl_table_present = False # this is used to tell if source of gff3
        # from 
        # https://github.com/yjx1217/LRSDAY/blob/master/scripts/mfannot2gff3.pl,
        # will note have translation table, -or- 
        # https://github.com/kbseah/mitonotate/blob/master/mfannot2gff.pl , will
        # have translation table used in mfannot execution noted on each line.
        for indx, current_line in enumerate(gff_text_lines):
            # first line will be comment with gff3 noted so skip. Also skip
            # if line starts with comment indicating symbol of `#`.
            if indx == 0  or current_line.startswith("#"):
                #current_line_num += 1 #no longer used since have indx
                continue
            else:
                if transl_table_present == False and (
                    ";transl_table=" in current_line):
                    transl_table_present = True
                line_start,line_end = get_start_end_gff_line(current_line)
                #print (line_start,line_end)
                if rnl_end < line_start:
                    line_to_insert_after = indx - 1
                    line_to_insert_before = indx
                    break # so will stop when first hits where the end of the 
                    # rnl is now less than the start of feature on current line
                elif indx == len(gff_text_lines) - 1 and line_end < rnl_start:
                    line_to_insert_after = indx
                    line_to_insert_before = "NOT APPLICABLE"


        '''
        # split off any header lines for putting back later. This will make 
        # handling first true data line of gff3 easier.
        end_of_header = False
        header = []
        main_lines_of_gff = []
        for line in (gff_text_lines):
            if line.startswith('#') and not end_of_header:
                header.append(line + "\n")
            elif not line.startswith('#') and not end_of_header:
                end_of_header = True
                main_lines_of_gff.append(line)
            elif end_of_header:
                main_lines_of_gff.append(line)
        #first exclude first line because it will be an edge case where there
        # isn't two lines to insert between
        first_feature_start = 
        #now handle the last line as it will also be an edge case where there
        # isn't two lines to insert between
        pass
        # Now handle the standard cases where there is a where there's two lines 
        #to insert between
        lines_to_insert_between = (line_to_insert_after,line_to_insert_before)
        '''



        # Update GFF3 with the new entry
        #REAL ENTRY (WITHOUT TABS SHOWING) AS MODEL ON NEXT LINE:
        #CMF_6-12611    mfannot rRNA    53945   57205   .   +   .   ID=rnl;Name=rnl;transl_table=3;gene=rnl
        # Note that "Fields must be tab-separated", according to 
        # https://useast.ensembl.org/info/website/upload/gff3.html .
        # Extract id_ at beginning of line and translation table value from
        # previous line, unless there is none
        if gff_text_lines[line_to_insert_after].startswith(
            "#") and not line_to_insert_before == "NOT APPLICABLE":
            id_ = gff_text_lines[line_to_insert_before].split()[0]
            if transl_table_present:
                transl_tbl = gff_text_lines[line_to_insert_before].split(
                    "transl_table=")[1].split(";")[0]
        elif not gff_text_lines[line_to_insert_after].startswith("#"):
            id_ = gff_text_lines[line_to_insert_after].split()[0]
            if transl_table_present:
                transl_tbl = gff_text_lines[line_to_insert_after].split(
                    "transl_table=")[1].split(";")[0]
        else:
            id_ = gff_text_lines[l-1].split()[0]
            if transl_table_present:
                transl_tbl = gff_text_lines[-1].split(
                    "transl_table=")[1].split(";")[0]
        if transl_table_present:
            new_entry = ("{}\tmfannot\trRNA\t{}\t{}\t.\t{}\t.\t"
                "ID=rnl;Name=rnl;transl_table={};gene=rnl".format(
                id_,rnl_start, rnl_end,strand, transl_tbl))
        else:
            new_entry = ("{}\tmfannot\trRNA\t{}\t{}\t.\t{}\t.\t"
                "ID=rnl;Name=rnl".format(
                id_,rnl_start, rnl_end,strand))
        updated_txt = ""
        #updated_txt = new_entry # FOR DEVELOPMENT ONLY
        updated_txt = make_updated_txt(gff_text_lines, line_to_insert_after,
            line_to_insert_before,new_entry)


        # Reporting and Saving (beginning)
        #--------------------------------------------------------------------
        # feedback
        fixed_fn = generate_output_file_name(gff_file_name, 
            text_to_add_to_altered_file)
        sys.stderr.write(
            "\n\n**rnl added and file saved as {}.**\n".format(fixed_fn))
        with open(fixed_fn, 'w') as output:
            output.write(updated_txt)
        status = "fixed"

    # Reporting and Saving continued (optional portion)
    #---------------------------------------------------------------------------

    # Do these no matter outcome, if needed. These are for when main function
    # called via import
    if return_status and return_quality_estimate:
        return status, quality
        sys.stderr.write("...status and quality score returned.\n")
    elif return_quality_estimate and not return_status:
            return quality
            sys.stderr.write("...quality score returned.\n")
    elif return_status and not return_quality_estimate:
            return status
            sys.stderr.write("...status returned.\n")
    # This was a dubm way to do this. I should have just built what to return 
    # as a list and then returned that! Fix when add more and then change how
    # feedback notifications done.
        



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
    fix_lsu_rRNA_annotation_in_gff_resulting_from_mfannot(gff_file_name, 
        seq_file_name,**kwargs)
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
    parser = argparse.ArgumentParser(prog='fix_lsu_rRNA_annotation_in_gff_resulting_from_mfannot.py',
        description="fix_lsu_rRNA_annotation_in_gff_resulting_from_mfannot.py \
        takes output from MFannot that has been converted to a gff file by a\
        subsequent script, a corresponding fungal mitochondrial sequence file, \
        and adds annotation for large ribosomal RNA if it is not already there.\
        **** Script by Wayne Decatur   \
        (fomightez @ github) ***")

    parser.add_argument("gff_file", help="Name of annotaion results \
        file (gff format) from a mitochondrial genome to parse and possibly \
        fix.\
        ", metavar="GFF_FILE")
    parser.add_argument("sequence_file", help="Name of file containing the \
        fungal mitochondrial sequence (FASTA format) corresponding to the \
        annotation file.\
        ", metavar="SEQ_FILE")




    #I would also like trigger help to display if no arguments provided because 
    # need at least one for url
    if len(sys.argv)==1:    #from http://stackoverflow.com/questions/4042452/display-help-message-with-python-argparse-when-script-is-called-without-any-argu
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    gff_file_name = args.gff_file
    seq_file_name = args.sequence_file


    main()

#*******************************************************************************
###-***********************END MAIN PORTION OF SCRIPT***********************-###
#*******************************************************************************
