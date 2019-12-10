#!/usr/bin/env python
# report_diff_between_two_seq_strings.py
__author__ = "Wayne Decatur" #fomightez on GitHub
__license__ = "MIT"
__version__ = "0.1.0"


# report_diff_between_two_seq_strings.py by Wayne Decatur
# ver 0.1
#
#*******************************************************************************
# Verified compatible with both Python 2.7 and Python 3.6; written initially in 
# Python 3. 
#
#
# PURPOSE: Takes two Python strings that correspond to biological sequences and
# reports differences between them as calculated from a global pairwise 
# alignment and optionally reports the length difference of the two unaligned 
# strings. 
# The idea for now is that this reported data can be used to gauge amount of 
# changes between matches to a sequence pattern in order to summarize a lot of
# matches to the pattern in many genomes. However, I am writing it as a general
# function and so it may be useable in several other contexts.
# 
#
#
# Based on `score_differences_between_sequences_by_pairwise_alignment.py` 
# directly, which was based on
# `roughly_score_relationships_to_subject_seq_pairwise_premsa.py`.
# Core processing of pairwaise alignmets based on 
# https://github.com/berrisfordjohn/adding_stats_to_mmcif/blob/master/adding_stats_to_mmcif/pairwise_align.py.
#
# As of writing, this script is intended for short sequences, but important to 
# keep in mind the following that is inherited from software used to handle the
# alignments/scoring:
# IF THE SEQUENCES ARE LONGER THAN THE 
# DOUBLE THE LIMIT THAT CAUSES THINGS TO SILENTLY FAIL, IT WILL SAMPLE THE 
# ALIGNMENT OF THAT LENGTH AT THE START AND THA SAME SIZE BACK FROM THE END AND 
# COMBINE THE SCORES. IF THE SEQUENCE . IF THE SEQUENCES ARE BETWEEN THE SIZE
# OF THE LIMIT AND DOUBLE THE LIMIT, it will only sample the alignment of the
# length of the limit so that there is no chance of double sampling. If the size
# of the sequences is less than the limit, it will align and score the entire
# sequence. It creates a dataframe ranking the sequences from most similar to 
# most different relative the first one in the supplied file, using whatever
# amount of sequence it was able to compare.
# TROUBLESHOOTING STEP#1: If you are working with moderately- to large-sized 
# sequences (greater than 5kb) and you are seeing the program hang for more 
# than a minute, make the size of the squence block, `block_len`, being 
# compared very small, say `100` to see if anything works.
#
# HUGE CAVEAT THAT CAUSED THIS SCRIPT TO BE ADJUSTED TO SIMPLY SAMPLE PART OF 
# SEQUENCES FOR LONG SEQUENCES (i.e., greater than 9151 bp): Won't work on long 
# sequences (fails silently and just hangs), it seems, due to memory use and 
# speed limitations. Even with just getting the score. See 
# https://github.com/biopython/biopython/pull/1655 ACtually, I think it is the
# reading in of all the sequences upfront that is compounding problems. However,
# I don't think it is failing at the reading in step because when tell it to 
# process less of sequences after that step, it works. So reading isn't the 
# limiting step. Maybe divide up and then read in only the two being processed 
# for each round? Other evidence that reading step isn't necessarily part of 
# issue is that removing all but two of the sequences, still has the issue 
# where the size of 9151 works and 9155 from the start doesn't work.
# Some info on runs on Binder:
# With limiting to first 5000 bp, runs instantly and works with thirteen 
# sequences of chromosomes of about 80 kb. With limiting to 
# first 9151 bp, takes a minute or so on Binder but still works. With limiting 
# to first 9155 or more though it fails on Binder. 
# Ughh. Next day anything greater than 9144 is failing to run?!?! So not 
# exactly consistent.
# But maybe this is good enough to get a general gauge of relatedness for a
# small chromosome? Actually when sequence is longer than double the limit, I 
# edited the code to do two passes, one at the start of all the sequences 
# spanning up to the size limit and one back from the end the same size. 
#
# 
#
#
#
#

#
#
# Dependencies beyond the mostly standard libraries/modules:
# Biopython
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
# python report_diff_between_two_seq_strings.py GAATTC GATTTC --report_len_diff
#-----------------------------------
# Issue `report_diff_between_two_seq_strings.py.py -h` for details.
# 
# More examples from running from the command line are at the links below: 
# https://git.io/????????  <-- for when make a demo
#
#
#
# To use this after importing/pasting or loading into a cell in a Jupyter 
# notebook, specify at least the results file (or results as a string) in the 
# call to the main function similar to below:
# from report_diff_between_two_seq_strings import report_diff_between_two_seq_strings
# diffs = report_diff_between_two_seq_strings("GAATTC","GATTTC ")
# 
#
# to adjust sequence block size used, 
# diffs = report_diff_between_two_seq_strings(
#   "GAATTC","GATTTC ", block_len = 15000) # NOT RECOMMENDED AS MIGHT CAUSE SILENT FAILURE (just hangs); ONLY TRY FOR LARGE MACHINES
#
# A more in-depth series of examples of using this script within a notebook 
# is found at:
# https://git.io/???????? <-- for when make a demo
#
# 
#
'''
CURRENT ACTUAL CODE FOR RUNNING/TESTING IN A NOTEBOOK WHEN IMPORTED/LOADED OR 
PASTED IN ANOTHER CELL:
from report_diff_between_two_seq_strings import report_diff_between_two_seq_strings
df = report_diff_between_two_seq_strings("seqs.fa")
df
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




default_block_length = 9141 # I found that in Jupyter sessions hosted on Binder
# that memory was an issue that causes things to silently fail if I tried using
# sequences 9141 bp or longer. FEEL FREE TO ADJUST THIS IF YOU ARE ON A LARGER 
# MACHINE. There are ways to change it for specific jobs on the command line, or 
# when you call the main function. Make this value tiny if things seem to hang for 
# more than five minutes so you can verify everything else is fine.

#
#*******************************************************************************
#**********************END USER ADJUSTABLE VARIABLES****************************













#*******************************************************************************
#*******************************************************************************
###DO NOT EDIT BELOW HERE - ENTER VALUES ABOVE###

import sys
import os
from Bio import Seq
from Bio import SeqIO
from Bio import Align
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
from collections import defaultdict



###---------------------------HELPER FUNCTIONS---------------------------------###


###--------------------------END OF HELPER FUNCTIONS---------------------------###
###--------------------------END OF HELPER FUNCTIONS---------------------------###


###---------------------------WORK HORSE CLASS-------------------------------###


# this comes straight from https://github.com/berrisfordjohn/adding_stats_to_mmcif/blob/master/adding_stats_to_mmcif/pairwise_align.py

class SequenceAlign:
    # this is modified to get mininum number of residue differences for the 
    # pairwise alignmnet from
    # https://github.com/berrisfordjohn/adding_stats_to_mmcif/blob/master/adding_stats_to_mmcif/pairwise_align.py

    def __init__(self, sequence1, sequence2):
        self.sequence1 = sequence1
        self.sequence2 = sequence2
        self.score = None
        self.differences = None
        #logging.debug(self.sequence1)
        #logging.debug(self.sequence2)

    def rna(self, seq):
        return set(seq).issubset(set("AUGC"))

    def dna(self, seq):
        return set(seq).issubset(set("ATGC"))

    # def dnarna(self, seq):
    #    return self.rna(seq=seq) or self.dna(seq=seq)

    def both_sequences_same_type(self):
        if self.dna(self.sequence1) == self.dna(self.sequence2) and self.rna(self.sequence1) == self.rna(
                self.sequence2):
            return True, ''
        return False, 'sequences not the same type'
        # if self.dnarna(self.sequence1) != self.dnarna(self.sequence2):
        #    return False, 'sequences not the same type'
        # return True, ''

    def remove_gaps(self, sequence):
        return str(sequence).replace("\n", "").replace(" ", "")

    def prepare_sequences(self):
        self.sequence1 = self.remove_gaps(self.sequence1)
        self.sequence2 = self.remove_gaps(self.sequence2)
        # if len(self.sequence1) > 2000 or len(self.sequence2) > 2000:
        #    return False, 'sequences too long. Please install emboss needle'
        return self.both_sequences_same_type()

    def pairwise2(self):

        matrix = matlist.blosum62
        gap_open = -10
        gap_extend = -0.5
        alns = pairwise2.align.globalds(self.sequence1, self.sequence2, matrix, gap_open, gap_extend)[0]
        #logging.info(pairwise2.format_alignment(*alns))
        #logging.info(alns)
        self.score = alns[2]

    def pairwise_aligner(self):

        aligner = Align.PairwiseAligner()
        aligner.open_gap_score = -10
        aligner.extend_gap_score = -0.5
        aligner.target_end_gap_score = 0.0
        aligner.query_end_gap_score = 0.0
        # aligner.match = 2
        # aligner.mismatch = -1
        alignments = aligner.align(self.sequence1, self.sequence2)
        differences_list = []
        alignments_by_differences = defaultdict(list) #keys are number of 
        # differences and values are a list of all alignmnets matching that 
        # number of differences
        for alignment in sorted(alignments):
            split_alignment = str(alignment).split()
            assert len(split_alignment) % 3 == 0, "should be divisible by three"
            #symbolic_alignment_summary = split_alignment[int(len(split_alignment)/3):(int(len(split_alignment)/3)*2)]
            # turns out each part of pairwise alignment, that being 
            # first aligned sequence, symbolic summary of alignment, and 
            # other aligned seqeunce is on separate line already, & so can just 
            # do:
            symbolic_alignment_summary = split_alignment[1]
            unaligned = symbolic_alignment_summary.count("X") #by limiting
            # to the symbolic summary of the alignment, I don't have to
            # worry about 'X's in sequence making count wrong. And I also
            # don't have to worry about dividing gap indicators out of 
            # double counting on the next line either.
            #print(unaligned)
            gaps = symbolic_alignment_summary.count("-")
            differences = unaligned + gaps
            differences_list.append(differences)
            alignments_per_differences = (
                alignments_by_differences[differences].append(alignment))
        align_score = aligner.score(self.sequence1, self.sequence2)
        self.score = align_score
        self.differences = min(differences_list)
        self.alignments_w_min_differences = (
            alignments_by_differences[self.differences])
    """
    def do_alignment_emboss(self):
        needle_cline.asequence = "asis:" + self.sequence1
        needle_cline.bsequence = "asis:" + self.sequence2
        stdout, stderr = needle_cline()
        result = [AlignIO.read(StringIO.StringIO(stdout), "emboss")]
        self.score = self.get_emboss_score(result[0].seq, result[1].seq)
    def get_emboss_score(self, seq1, seq2):
        return sum(aa1 == aa2 for aa1, aa2 in zip(seq1, seq2))
    """

    def get_alignment_differences(self):
        return self.differences

    def do_sequences_align(self):
        if self.score > 0:
            return True
        return False

    def do_sequence_alignment(self):
        sequences_ok, error = self.prepare_sequences()
        if not sequences_ok:
            return False, error, max(
            len(self.sequence1),len(self.sequence2))
        self.pairwise_aligner()
        if self.do_sequences_align():
            return True, '', self.differences, self.alignments_w_min_differences
        return False, 'sequences do not align', max(
            len(self.sequence1),len(self.sequence2))
###--------------------------END OF WORK HORSE CLASS-------------------------###
###--------------------------END OF WORK HORSE CLASS-------------------------###

#*******************************************************************************
###------------------------'main' function of script---------------------------##

def report_diff_between_two_seq_strings(
    seq_string1, seq_string2, report_len_diff = False, report_alignment = False,
    block_len=default_block_length):
    '''
    Main function of script. 
    Takes two Python strings that correspond to biological sequences and reports 
    differences between them as calculated from a global pairwise alignment and 
    optionally reports the length difference of the two unaligned strings. 

    Setting `report_len_diff` to True is to be used when you want to report the 
    length difference of the two unaligned strings as well. The size difference
    will be returned second.

    Setting `report_alignment` to True is to be used when you want to report the 
    alignments corresponding to the reported differences.

    `block_len` is for dealing with large sequences and hopefully isn't 
    necessary if sticking with small ones as intended right now. (See 
    `roughly_score_relationships_to_subject_seq_pairwise_premsa.py` for more 
    about that)

    Returns either one,two, or three items:
    The first item returned is always the amount of differences.
    If `report_len_diff` is True, the unaligned size difference will be returned
    following the number of differences.
    If `report_report_alignment` is True, formatted text strings of the 
    alignment(s) associated with the minimal number of differences will be 
    returned last.
    '''
    # Compare unaligned if `report_len_diff` True
    if report_len_diff:
        len_diff = max([len(seq_string1), len(seq_string2)]) - min([len(seq_string1), len(seq_string2)])


    # Limit to handle if exceeds memory limits for aligning
    max_seq_size = max([len(seq_string1), len(seq_string2)])
    if max_seq_size > block_len:
        records = [x[:block_len] for x in records]
        sys.stderr.write("\nAnticipated memory issues with long sequence and\n"
            "so only block of {} bps from the start "
            "compared.\n...".format(len(records[0])))
    

    #Pairwise comparison for the sequence pair
    sa = SequenceAlign(sequence1=seq_string1, sequence2=seq_string2)
    aligned, error, differences, alignments_w_min_differences = (
        sa.do_sequence_alignment())
    #print(differences)  # for debugging/ development
    #print(alignments_per_differences)  # for debugging/ development

    '''
    # Set up for reporting alignment(s), if needed
    if report_alignment:
        formatted_alignments = format_alignments(differences,
            alignments_w_min_differences)
    '''


    # Reporting and saving(?)/feedback
    #---------------------------------------------------------------------------
    if report_len_diff and (not __name__ == "__main__"):
        to_return = (differences,len_diff)
    elif not __name__ == "__main__":
        to_return = (differences,) # see 
        # https://wiki.python.org/moin/TupleSyntax for creating single-item 
        # tuple because "the essential element here is the trailing comma". 
        # Otherwise `to_return = (differences)` results in `to_return` being an 
        # int that cannot be iterated over for unpacking to possibly add more 
        # later or return. (was considering what was better for building what to 
        # return and read "I'd suggest that a keyed-access return value starts 
        # making more sense than a tuple when there are more than about three 
        # member" in https://softwareengineering.stackexchange.com/a/327761), 
        # which would be good to keep in mind. (This all came up when I was
        # refactoring to add returning POSSIBLY an additional item after coding 
        # around only returning one or two items.)
    elif report_len_diff and __name__ == "__main__":
        # save a text file of the diff_score and len_diff?
        sys.stderr.write("\n\nReporting amount of differences and difference "
            "in length of the two input sequences...")
        print (differences,len_diff)
    else:
        # save a text file of the diff_score?
        sys.stderr.write("\n\nReporting amount of differences...")
        print (differences)
    if report_alignment and (not __name__ == "__main__"):
        to_return = (
            *to_return, [str(x) for x in alignments_w_min_differences])
    elif report_alignment:
        aln_txt_end = "ments"
        if len(alignments_w_min_differences) == 1:
            aln_txt_end = "ment"
        sys.stderr.write(
            f"\nReporting align{aln_txt_end} with {differences} differences...")
        print(f"Align{aln_txt_end} displaying {differences} differences:\n")
        print(*alignments_w_min_differences,sep='\n') # based on 
        # https://stackoverflow.com/a/52097312/8508004

    # once what to return is built entirely, return it
    if (not __name__ == "__main__"):
        return to_return


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
    kwargs['report_len_diff'] = report_len_diff
    kwargs['report_alignment'] = report_alignment
    report_diff_between_two_seq_strings(args.seq_string1, 
        args.seq_string2, **kwargs)
    # using https://www.saltycrane.com/blog/2008/01/how-to-use-args-and-kwargs-in-python/#calling-a-function
    # to build keyword arguments to pass to the function above
    # (see https://stackoverflow.com/a/28986876/8508004 and
    # https://stackoverflow.com/a/1496355/8508004 
    # (maybe https://stackoverflow.com/a/7437238/8508004 might help too) for 
    # related help). Makes it easy to add more later.





if __name__ == "__main__":
    ###-----------------for parsing command line arguments-----------------------###
    import argparse
    parser = argparse.ArgumentParser(prog='report_diff_between_two_seq_strings.py',
        description="report_diff_between_two_seq_strings.py \
        Takes two Python strings that correspond to biological sequences and \
        reports differences between them as calculated from a global pairwise \
        alignment and optionally reports the length difference of the two \
        unaligned strings. \
        **** Script by Wayne Decatur   \
        (fomightez @ github) ***")

    parser.add_argument("seq_string1", help="text representing first \
        biological sequence.\
        ", metavar="SEQ1")
    parser.add_argument("seq_string2", help="text representing other \
        biological sequence.\
        ", metavar="SEQ2")
    parser.add_argument("-rld", "--report_len_diff",help=
    "Add this flag to report the length difference of the unaligned sequences \
    additionally following the amount of differences.",
    action="store_true")
    parser.add_argument("-ra", "--report_alignment",help=
    "Add this flag to report the alignment(s) associated with the minimal \
    number of differences.", action="store_true")



    #I would also like trigger help to display if no arguments provided because 
    # need at least one for url
    if len(sys.argv)==1:    #from http://stackoverflow.com/questions/4042452/display-help-message-with-python-argparse-when-script-is-called-without-any-argu
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    report_len_diff = args.report_len_diff
    report_alignment = args.report_alignment


    main()

#*******************************************************************************
###-***********************END MAIN PORTION OF SCRIPT***********************-###
#*******************************************************************************
