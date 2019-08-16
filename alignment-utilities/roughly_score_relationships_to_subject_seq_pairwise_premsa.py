#!/usr/bin/env python
# roughly_score_relationships_to_subject_seq_pairwise_premsa.py
__author__ = "Wayne Decatur" #fomightez on GitHub
__license__ = "MIT"
__version__ = "0.1.0"


# roughly_score_relationships_to_subject_seq_pairwise_premsa.py by Wayne Decatur
# ver 0.1
#
#*******************************************************************************
# Verified compatible with both Python 2.7 and Python 3.6; written initially in 
# Python 3. 
#
#
# PURPOSE: Takes a file of multiple sequences in FASTA format and aligns each
# of them in turn to the first sequence. IF THE SEQUENCES ARE LONGER THAN THE 
# DOUBLE THE LIMIT THAT CAUSES THINGS TO SILENTLY FAIL, IT WILL SAMPLE THE 
# ALIGNMENT OF THAT LENGTH AT THE START AND THA SAME SIZE BACK FROM THE END AND 
# COMBINE THE SCORES. IF THE SEQUENCE . IF THE SEQUENCES ARE BETWEEN THE SIZE
# OF THE LIMIT AND DOUBLE THE LIMIT, it will only sample the alignment of the
# length of the limit so that there is no chance of double sampling. If the size
# of the sequences is less than the limit, it will align and score the entire
# sequence. It creates a dataframe ranking the sequences from most similar to 
# most different relative the first one in the supplied file, using whatever
# amount of sequence it was able to compare.
# The dataframe is saved as a tabular text file when used on the command line. 
# Optionally, it can also return that dataframe for use inside a Jupyter 
# notebook.
# Meant to be used before aligning in order to possibly divide up jobs or screen 
# for which sequences to leave out as repeats or overly highly related 
# sequences that add nothing more to the alignment.
# 
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
# NOTE ON OTHER APPROACH TO LIMIT MEMORY REQUIREMENTS:
# To save time, no alignment data is actually produced. 
# "This improves memory usage and speed., see 
# https://github.com/berrisfordjohn/adding_stats_to_mmcif/blob/master/adding_stats_to_mmcif/pairwise_align.py
# and July 11, 2018 Biopython User's list email.
#
#
#
#
#
# Core processing of pairwaise alignmets based on 
# https://github.com/berrisfordjohn/adding_stats_to_mmcif/blob/master/adding_stats_to_mmcif/pairwise_align.py.
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
# python roughly_score_relationships_to_subject_seq_pairwise_premsa.py -SEQS_FILE
#-----------------------------------
# Issue `roughly_score_relationships_to_subject_seq_pairwise_premsa.py.py -h` for details.
# 
# More examples from running from the command line are at the links below: 
# https://git.io/???????? <-- for when make a demo
#
#
#
# To use this after importing/pasting or loading into a cell in a Jupyter 
# notebook, specify at least the results file (or results as a string) in the 
# call to the main function similar to below:
# from roughly_score_relationships_to_subject_seq_pairwise_premsa import roughly_score_relationships_to_subject_seq_pairwise_premsa
# df = roughly_score_relationships_to_subject_seq_pairwise_premsa("seqs.fa")
# df
#
# to adjust sequence block size used, 
# df = roughly_score_relationships_to_subject_seq_pairwise_premsa(
#   "seqs.fa", block_len = 15000) # NOT RECOMMENDED AS MIGHT CAUSE SILENT FAILURE (just hangs); ONLY TRY FOR LARGE MACHINES
#
# A more in-depth series of examples of using this script within a notebook 
# without need to save file intermediates is found at:
# https://git.io/???????? <-- for when make a demo
#
# 
#
'''
CURRENT ACTUAL CODE FOR RUNNING/TESTING IN A NOTEBOOK WHEN IMPORTED/LOADED OR 
PASTED IN ANOTHER CELL:
from roughly_score_relationships_to_subject_seq_pairwise_premsa import roughly_score_relationships_to_subject_seq_pairwise_premsa
df = roughly_score_relationships_to_subject_seq_pairwise_premsa("seqs.fa")
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

## Settings and options for output plot 
df_save_as_name = 'ranked_seqs_df.tsv' # name for saving dataframe as text table

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



###---------------------------HELPER FUNCTIONS---------------------------------###


###--------------------------END OF HELPER FUNCTIONS---------------------------###
###--------------------------END OF HELPER FUNCTIONS---------------------------###


###---------------------------WORK HORSE CLASS-------------------------------###


# this comes straight from https://github.com/berrisfordjohn/adding_stats_to_mmcif/blob/master/adding_stats_to_mmcif/pairwise_align.py

class SequenceAlign:
    # this comes straight from 
    # https://github.com/berrisfordjohn/adding_stats_to_mmcif/blob/master/adding_stats_to_mmcif/pairwise_align.py

    def __init__(self, sequence1, sequence2):
        self.sequence1 = sequence1
        self.sequence2 = sequence2
        self.score = None
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
        # only need to run aligner.score. This improves memory usage and speed. 
        # alignments = aligner.align(self.sequence1, self.sequence2)
        # for alignment in sorted(alignments):
        #    logging.debug(alignment)
        align_score = aligner.score(self.sequence1, self.sequence2)
        #logging.info(align_score)

        self.score = align_score

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

    def get_alignment_score(self):
        return self.score

    def do_sequences_align(self):
        if self.score > 0:
            return True
        return False

    def do_sequence_alignment(self):
        sequences_ok, error = self.prepare_sequences()
        if not sequences_ok:
            return False, error, 0
        self.pairwise_aligner()
        if self.do_sequences_align():
            return True, '', self.score
        return False, 'sequences do not align', 0
###--------------------------END OF WORK HORSE CLASS-------------------------###
###--------------------------END OF WORK HORSE CLASS-------------------------###

#*******************************************************************************
###------------------------'main' function of script---------------------------##

def roughly_score_relationships_to_subject_seq_pairwise_premsa(
    seqs, block_len=default_block_length, return_df = True, save_text_of_df=True):
    '''
    Main function of script. 
    sequences scored of how similar the sequences are to first one in provided
    sequence file.
    based on using 
    https://github.com/berrisfordjohn/adding_stats_to_mmcif/blob/master/adding_stats_to_mmcif/pairwise_align.py

    Optionally also returns a dataframe of the results data, meant for use in a
    Jupyter notebook.
    '''
    # Read in the sequences
    records = []
    for record in SeqIO.parse(seqs, "fasta"):
        records.append(record)
    # feedback
    sys.stderr.write("Sequences read in...")
    max_seq_size = max([len(x.seq) for x in records])
    sys.stderr.write("\nLongest sequence in input detected as "
        "{}.\n...".format(max_seq_size))


    # Prepare to pairwise comparisons to the first one in the list of supplied sequences
    sys.stderr.write("calculating scores of pairwise alignments...")
    subj = records[0]
    sequences_to_compare = records[1:]

    # Limit to handle if exceeds memory limits for aligning
    if max_seq_size > block_len:
        subj = subj[:block_len]
        sequences_to_compare = [x[:block_len] for x in sequences_to_compare]
        sys.stderr.write("\nAnticipated memory issues with long sequence and\n"
            "so only block of {} bps from the start "
            "compared.\n...".format(len(sequences_to_compare[0])))
    
    # clearing records help any on performance? Answer seemed to be 'no'.
    #records = [] # didn't seem to hep as still failed at same point in test
    # with or without that step. Plus this will mean cannot add in another
    # block of aligned sequence for assessment below.


    score_results = {}
    for sequence in sequences_to_compare:
        sa = SequenceAlign(sequence1=subj.seq, sequence2=sequence.seq)
        aligned, error, score = sa.do_sequence_alignment()
        score_results[sequence.id] =  score


    # for those longer than 2x the size of the block the memory resources seem
    # to handle, also assess a block of sequence of the safe length back from 
    # the end to get an even larger representation of the diversity of the 
    # sequence (not going to divide up and assess further because beyond ends
    # might not match up well without alignment anyway and so no point. THIS IS
    # WHOL ENDEAVOR IS MEANT AS QUCIK AND DIRTY ASESSMENT prior to real 
    # alignment of at least some of the sequences and SO NO NEED TO BE
    # THOROUGH.)
    if max_seq_size > (block_len*2):
        subj = records[0]
        sequences_to_compare = records[1:]
        subj = subj[-block_len:]
        sequences_to_compare = [x[-block_len:] for x in sequences_to_compare]
        sys.stderr.write("\nSuper long sequence detected: Comparing a 2nd "
            "block back from the\nsequence 'end' as well and "
            "combining scores.\n...")
        for sequence in sequences_to_compare:
            sa = SequenceAlign(sequence1=subj.seq, sequence2=sequence.seq)
            aligned, error, score = sa.do_sequence_alignment()
            score_results[sequence.id] +=  score


    # Make dataframe of the results and sort on score
    sys.stderr.write("summarizing scores...")
    import pandas as pd
    df = pd.DataFrame(list(score_results.items()), columns = ['id', 'score_vs_'+subj.id])
    df = df.sort_values('score_vs_'+subj.id, ascending=False)
        

    # feedback
    sys.stderr.write("\nResults converted to a dataframe...")


    # Reporting and Saving
    #---------------------------------------------------------------------------
    #print(df)#originally for debugging during development,added..
    # Document the full set of data collected in the terminal or 
    # Jupyter notebook display in some manner. 
    # Using `df.to_string()` because more universal than `print(df)` 
    # or Jupyter's `display(df)`.
    #sys.stderr.write("\nFor documenting purposes, the following lists the "
    #    "parsed data:\n")
    #with pd.option_context('display.max_rows', None, 'display.max_columns', None):
    #    display(df)
    #sys.stderr.write(df.to_string())

    # Handle saving the dataframe
    if save_text_of_df == False:
        sys.stderr.write("\n\nTabular text of the data "
        "was not stored for use\nelsewhere "
        "because `no_table` was specified.")
    else:
        df.to_csv(df_save_as_name, sep='\t',index = False)
        # Let user know
        sys.stderr.write("\n\nA table of the data "
        "has been saved as a text file (tab-delimited).\n"
        "DATA is stored as ==> '{}'".format(df_save_as_name ))

    
    # Return dataframe (optional)
    #---------------------------------------------------------------------------
    if return_df:
        sys.stderr.write("\n\nReturning a dataframe with the information "
                "as well.")
        return df

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
    kwargs['block_len'] = block_len 
    if df_save_as_name == 'no_table':
        kwargs['save_text_of_df'] = False
    kwargs['return_df'] = False #probably don't want dataframe returned if 
    # calling script from command line
    roughly_score_relationships_to_subject_seq_pairwise_premsa(seqs,**kwargs)
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
    ###-----------------for parsing command line arguments-----------------------###
    import argparse
    parser = argparse.ArgumentParser(prog='roughly_score_relationships_to_subject_seq_pairwise_premsa.py',
        description="roughly_score_relationships_to_subject_seq_pairwise_premsa.py \
        Takes a file of multiple sequences in FASTA format and aligns each of \
        them in turn to the first sequence. However, if the sequence happens \
        to be moderate- or large-sized (> 5 kb), by default it only samples \
        part of the sequence due to memory limitations. It scores the \
        alignments and \
        produces a dataframe ranking the sequences from most similar to most \
        different relative the first one in the supplied file. The dataframe \
        is saved as a tabular text file when used on the command line. \
        Optionally, it can also return that dataframe for use \
        inside a Jupyter notebook. \
        **** Script by Wayne Decatur   \
        (fomightez @ github) ***")

    parser.add_argument("seqs_file", help="Name of file of sequences \
        (all FASTA-formatted) to compare to the first in that file .\
        ", metavar="SEQS_FILE")
    parser.add_argument('-bl', '--block_len', action='store', type=int, 
    default= default_block_length, help="**FOR ADVANCED USE.*** Allows for \
    setting the sequence blocks compared for long sequences. The default of \
    {} was worked out for uisng in a session launched form MyBinder.org 2018. \
    You can free to make larger if you have more computational resources. For \
    example, `-bl 15000`; however, the failing condition is just a \
    silent(hanging) state.".format(default_block_length)) 
    parser.add_argument('-dfo', '--df_output', action='store', type=str, 
    default= df_save_as_name, help="OPTIONAL: Set file name for saving tabular \
    text (tab-delimited) derived from the produced dataframe. If none \
    provided, '{}' will be used. To force no table to \
    be saved, enter `-dfo no_table` without quotes or ticks as output file \
    (ATYPICAL).".format(df_save_as_name))



    #I would also like trigger help to display if no arguments provided because 
    # need at least one for url
    if len(sys.argv)==1:    #from http://stackoverflow.com/questions/4042452/display-help-message-with-python-argparse-when-script-is-called-without-any-argu
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    seqs = args.seqs_file
    block_len = args.block_len
    df_save_as_name = args.df_output


    main()

#*******************************************************************************
###-***********************END MAIN PORTION OF SCRIPT***********************-###
#*******************************************************************************
