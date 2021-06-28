# make_report_on_equivalents_of_interacting_residues.py
__author__ = "Wayne Decatur" #fomightez on GitHub
__license__ = "MIT"
__version__ = "0.1.0"


# make_report_on_equivalents_of_interacting_residues.py by Wayne Decatur
# ver 0.1.0
#
#*******************************************************************************
# 
# PURPOSE: For a protein that has a structure that you've used to find a remote
# homolog using the HH-suite 3 software, this script reports if residues of that 
# protein (chain #1) that interact with another specified chain have equivalents
# in the identified remote homolog.
# Needs to work in conjunction with the notebook 
# `Report if residues interacting with a specific chain have equivalent residues in an hhsuite-generated alignment.ipynb` 
# that is presently in https://github.com/fomightez/hhsuite3-binder . In fact, 
# the easiest way to use this is to launch sessions by clicking on the 
# `launch binder` badge at that repo. In the session that comes up, everything 
# will already be installed and available for working through the notebook 
# `Report if residues interacting with a specific chain have equivalent residues in an hhsuite-generated alignment.ipynb` 
# that does this comparison for a demonstration protein with a structure & 
# identified homolog. Users can then change the PDB codes and chain designations
# and add their HH-suite3 to analyze their own structures, homologs, and 
# potential chain interactions of interest.
# 
#
#
#
#
#
# to do: 
# - decide if need dataframe makeing stuff in 
# `GET NECESSARY COMPANION SCRIPTS AND IMPORT FUNCTIONS` and if not, 
# comment out BUT KEEP HERE BECAUSE MAY NEED FOR SNAKEMAKE worflow!!!

import os
import sys
import glob
from shutil import copyfile
import subprocess
import pandas as pd
import numpy as np
#from halo import HaloNotebook as Halo
from IPython.utils import io





################################################################################
#######----------------------HELPER FUNCTIONS-----------------------------######
from contextlib import contextmanager,redirect_stderr,redirect_stdout
from os import devnull

@contextmanager
def suppress_stdout_stderr():
    """
    A context manager that redirects stdout and stderr to devnull.
    From https://stackoverflow.com/a/52442331/8508004
    """
    with open(devnull, 'w') as fnull:
        with redirect_stderr(fnull) as err, redirect_stdout(fnull) as out:
            yield (err, out)

def res_tuple_simple(items):
    '''
    takes columns with Atom#1 residue number and chain designation along with
    Atom#2 residue number and chain designation
    and returns a tuple of residue number 1 position number first followed by 
    residue number 2 position number.
    '''
    # check for 'empty' dataframe using `np.isnan(items[0])`, based on https://stackoverflow.com/a/29528160/8508004 ; allows rest of script to run
    # gracefully if there were no interactions between the two chains examined
    # for one of the sturctures and thus 'empty' dataframe produced from PDBsum
    # data
    if np.isnan(items[0]):
        return ("0","0")
    return ("{}".format(items[0]),"{}".format(items[2]))
def res_tuple(items):
    '''
    takes columns with Atom#1 residue number and chain designation along with
    Atom#2 residue number and chain designation
    and returns a tuple of residue#1 information first followed by 
    residue#2 information.
    Uses Jmol/Jsmol convention where `161:B` means residue #161 of chain B.
    '''
    # check for 'empty' dataframe using `np.isnan(items[0])`, based on https://stackoverflow.com/a/29528160/8508004 ; allows rest of script to run
    # gracefully if there were no interactions between the two chains examined
    # for one of the sturctures and thus 'empty' dataframe produced from PDBsum
    # data
    if np.isnan(items[0]):
        return ("0:NA","0:NA")
    return ("{}:{}".format(items[0],items[1]),"{}:{}".format(items[2],items[3]))

def write_string_to_file(s, fn):
    '''
    Takes a string, `s`, and a name for a file & writes the string to the file.
    '''
    with open(fn, 'w') as output_file:
        output_file.write(s)

def out2_stderr_n_log(s,log_file_text):
    '''
    Takes a string as input and sends it to the stderr as well as to a building
    string that will everntually get saved as a Log file.
    Also needs the Log file to be sent in because gets assigned within the
    function in order to add to it. Returns the modified `log_file_text`.
    '''
    sys.stderr.write(s)
    log_file_text += s
    return log_file_text


import time
from IPython.display import display, Javascript
import hashlib
def save_notebook(file_path):
    '''
    Function to save a notebook from 
    https://stackoverflow.com/a/57814673/8508004

    IMPORTANTLY, this won't work in the JupyterLab interface for notebooks!
    See https://github.com/jupyterlab/jupyterlab/issues/7627
    '''
    start_md5 = hashlib.md5(open(file_path,'rb').read()).hexdigest()
    display(Javascript('IPython.notebook.save_checkpoint();'))
    current_md5 = start_md5
    while start_md5 == current_md5:
        time.sleep(1)
        current_md5 = hashlib.md5(open(file_path,'rb').read()).hexdigest()



def chunk_string(string, chunk_size):
    """Return a list of n-sized chunks from string of letters."""
    return [string[i:i+chunk_size] for i in range(0, len(string),chunk_size)] 


def strip_off_first_line(fn,set_name,character_to_mark_set_name_end):
    '''
    This takes a name of a file & then uses the shell to remove the first line.
    In order to leave the input file intact, a new multi-sequence FASTA file
    is made and that is used in place of the one where the label was the first
    line. The set sample name extracted gets added to the file name.
    Removing first line based on 
    https://unix.stackexchange.com/questions/96226/delete-first-line-of-a-file
    '''
    name_for_f_without_first_line = (
        f"{set_name}{character_to_mark_set_name_end}set.fa")
    #!tail -n +2 {fn} >{name_for_f_without_first_line} 
    os.system(f"tail -n +2 {fn} >{name_for_f_without_first_line}")
    return name_for_f_without_first_line


def percent_GCcalc(items):
    '''
    takes a list of three and calculates percentage of sum of first
    two itemswithin total (second item)

    Taken from 
    `GSD Adding_percentGC_to_nt_counts_for_mito_genomes_from_1011_collection.ipynb`
    '''
    return (items[0] + items[1])/items[2]

def calculate_end_position(start_pos,aligned_sequence):
    '''
    Takes a `start position` corresponding to the first residue number of the 
    `aligned_sequence` that also gets provided and calculates the position 
    represented by the final residue in the aligned_sequence, accounting for
    gaps that may have been added to the aligned sequence by removing them to
    get length when not included.
    '''
    return start_pos + len(aligned_sequence.replace("-",""))


def get_aln_index_and_real_pos(sequence):
    '''
    a generator that takes a sequence and then returns next position in 
    alignment (equivalent to index, i.e., zero-indexed) and the actual position 
    in the contiguous sequence to which that corresponds. So the second-value
    is the ungapped position in common terms.
    For example, if sequence is `a-b`,the first values returned are `0,1`.
    The second values returned are `1,1`. And third returned values are `2,2`.
    Appears previously in several of my scripts:
    `extract_regions_from_clustal_alignment`
    `calculate_cons_for_clustal_nucleic.py`
    `calculate_cons_for_clustal_protein.py`
    `score_columns_in_clustal_msa.py`
    `score_sequences_in_clustal_msa_favoring_top_line.py`
    `score_sequences_in_clustal_msa.py `
    However, only `extract_regions_from_clustal_alignment_v0.1 VERY SLOW IF ALIGNMENT OVER 10K.py`
    seemed to use it. Looks like I left in others, in case needed. But the slow
    version was tossed out and I used Biopython to make the ungapped. Don't 
    recall if absolutely needed biopython over this generation to speed up for
    larger sequences or was part of overall reworking that made it faster.
    Sequences here won't be overly long since protein sequences don't have as
    many residues as a chromosome may have nucleotides.
    '''
    indx = 0
    while True:
        yield indx, len(sequence[:indx+1].replace("-","")) #because second 
        # value after colon means up to but not including, the `+1` there allows
        # getting first character when index is set to zero
        indx+=1


#######------------------END OF HELPER FUNCTIONS--------------------------######
################################################################################









#### SETTINGS FOR THE MAIN PART OF THE SCRIPT-----------------------------------
num_flanking_to_get_for_details = 5 #number of residues to show on each side in report 
# where context of equivalents are shown

report_details_for_those_with_equivalents = True # show context with a little
# flanking sequence. Can override this default from notebook calling script by
# providing `details_for_those_with_equivalents`.
report_details_for_those_without_equivalents = False # show something for the 
# context of those where there are gaps. Can override this default from notebook 
# calling script by providing `details_for_those_without_equivalents`.


# Conservation calling-related settings lifted from 
# `calculate_cons_for_clustal_protein.py`:
strongly_similar_aa_tuples = [("S","T","A"), ("N","E","Q","K"), 
("N","H","Q","K"), ("N","D","E","Q"), ("Q","H","R","K"), ("M","I","L","V"),
("M","I","L","F"),("H","Y"),("F","Y","W")] #used to decide if
# conservative substition (strongly similar) in case where not all identical.
# Based on
# https://www.ebi.ac.uk/seqdb/confluence/display/JDSAT/Bioinformatics+Tools+FAQ#BioinformaticsToolsFAQ-WhatdoconsensussymbolsrepresentinaMultipleSequenceAlignment?

weakly_similar_aa_tuples = [("C","C","A"), ("A","T","V"), ("S","A","G"), 
("S","T","N","K"), ("S","T","P","A"), ("S","G","N","D"), 
("S", "N", "D", "E","Q","K"), ("N", "D", "E", "Q","H","K"),
("N", "E", "H", "Q","R","K"), ("F", "V", "I", "L","M"), ("H","F","Y")] #used to 
# decide if conservative substition (weakly similar) in case where not identical.
# Based on
# https://www.ebi.ac.uk/seqdb/confluence/display/JDSAT/Bioinformatics+Tools+FAQ#BioinformaticsToolsFAQ-WhatdoconsensussymbolsrepresentinaMultipleSequenceAlignment?

#### END OF SETTINGS FOR THE MAIN PART OF THE SCRIPT----------------------------











################################################################################
#######------------------------MAIN SECTION-------------------------------######

#spinner = Halo(text='Processing...', spinner='dots',color = 'magenta')
#spinner.start()

# GET NECESSARY COMPANION SCRIPTS AND IMPORT FUNCTIONS:
#------------------------------------------------------------------------------#
### WAIT! DATAversionFRAMES SHOULD BE MADE IN NOTEBOOK UNLESS I CHANGE SO NOT NEEDED?!?!
file_needed = "pdbsum_prot_interactions_list_to_df.py"
if not os.path.isfile(file_needed):
    sys.stderr.write("\nObtaining script containing a function to use to parse "
        "the data files from PDBsum "
        "...\n")
  # based on http://amoffat.github.io/sh/
    from sh import curl
    curl("-OL",
        "https://raw.githubusercontent.com/fomightez/structurework/master"
        "/pdbsum-utilities/"+file_needed)
    # verify that worked & ask for it to be done manually if fails
    if not os.path.isfile(file_needed):
        github_link = ("https://github.com/fomightez/structurework/tree"
            "/master/pdbsum-utilities")
        sys.stderr.write("\n'+file_needed+' not found. "
            "Please add it to your current working\ndirectory from {}"
            ".\n**EXITING !!**.\n".format(github_link))
        sys.exit(1)
from pdbsum_prot_interactions_list_to_df import pdbsum_prot_interactions_list_to_df

file_needed = "hhsuite3_results_to_df.py"
if not os.path.isfile(file_needed):
    sys.stderr.write("\nObtaining script containing a function to use to parse "
        "the data files from .hhr files made by HH-suite3 software "
        "...\n")
  # based on http://amoffat.github.io/sh/
    from sh import curl
    curl("-OL",
        "https://raw.githubusercontent.com/fomightez/sequencework/master"
        "/hhsuite3-utilities/"+file_needed)
    # verify that worked & ask for it to be done manually if fails
    if not os.path.isfile(file_needed):
        github_link = ("https://github.com/fomightez/sequencework/tree"
            "/master/hhsuite3-utilities")
        sys.stderr.write("\n'+file_needed+' not found. "
            "Please add it to your current working\ndirectory from {}"
            ".\n**EXITING !!**.\n".format(github_link))
        sys.exit(1)
from hhsuite3_results_to_df import hhsuite3_results_to_df


# IDENTIFY THE RESIDUE NUMBERS OF CHAIN #1 THAT INTERACT WITH OTHER CHAINS:
#------------------------------------------------------------------------------#
sys.stderr.write("Parsing data from PDBsum to identify chain {} residues "
    "interacting with chain {} residues ...\n".format(
    structure_chain1,structure_chain2))
interacting_res_nums = i_df["Atom1 Res no."].tolist()
# Remove any duplicate residue numbers.
interacting_res_nums = set(interacting_res_nums)


# FILTER OUT RESIDUES THAT AREN'T REPRESENTED IN THE ALIGNED SEQUENCES
#------------------------------------------------------------------------------#
# filter out residues that aren't represented in the aligned sequences because
# cannot say anything about those residues relative equivalents.
# To make it easier up front, I only asked for the start position of the 
# sequence to be provided and not the end position. However, because as part of
# this the end position can also be easily determined by determining how many
# actual residues are in the alignment. Cannot simply use the length of the 
# aligned sequence because there may be gaps, and so those have to be accounted
# for in this calculation. Easiest thing is just to remove the `-` symbols and
# then count length.
# Currently the notebook used to run this, provides the following:
# start_pos_query_seq
# start_pos_hit_seq
# aligned_query_seq
# aligned_hit_seq
sys.stderr.write("Identify any residues that are not in the span that is"
    " aligned ...\n")
# determine the end position of the query and the hit sequence in the alignment 
# using the sequence and the start position. Account for gaps.
qend = calculate_end_position(start_pos_query_seq,aligned_query_seq)
hend = calculate_end_position(start_pos_hit_seq,aligned_hit_seq)
# using the start and now-determined end, collect those that before the start 
# position and after the end position. Also collect those in the range that
# spans the portion represented in aligned query sequence.
interacting_res_nums_with_data = []
interacting_res_nums_without_data = []
for x in interacting_res_nums:
    if start_pos_query_seq <= x <= qend:
        interacting_res_nums_with_data.append(x)
    else:
        interacting_res_nums_without_data.append(x)


# SEE IF INTERACTING RESIDUE NUMBERS HAVE EQUIVALENTS IN THE aligned_subj_seq
#------------------------------------------------------------------------------#
# Going forward, only concerned with residues where the alignment provides data.
# For now just identify which have equivalent and which don't. Because maybe 
# none have equivalents and then everything could finish. Doesn't make sense to
# do anything else for this iteration over the list if no equivalents.
sys.stderr.write("Checking for equivalents of the interacting residues in "
    "\nin the other aligned sequence ...\n")
interacting_res_nums_without_equivs = []
interacting_res_nums_with_equivs = []
interacting_res_nums_without_equivs_details = {}
interacting_res_nums_with_equivs_details = {}
for res_num in interacting_res_nums_with_data:
    # Need to determine location corresponding to position of the residue 
    # contiguous sequence in the alignment, disregarding gaps. Will use a
    # generator `get_aln_index_and_real_pos()` I've used in the past that gives 
    # both the current index in thesequence being scanned and the real position 
    # of the sequence at that point if you disregard the gaps.
    # In seeing if that real position matches the interacting residue, also need
    # to take into account the sequences don't necessarily start at the first 
    # position depending on how much is spanned by the aligned sequence.
    # Note that need the -1 in defining `real_pos_in_aligned_seq` below  b/c
    # indexing going on is considering Python zero based indexing and numbers 
    # being provided from PDBsum are typical, common numbering where first 
    # residue in sequence is numbered one.
    for i,num_into_aligned_seg in get_aln_index_and_real_pos(aligned_query_seq):
        real_pos_in_aligned_seq = num_into_aligned_seg + start_pos_query_seq - 1
        if real_pos_in_aligned_seq == res_num:
            break # don't want to go any further b/c have position needed and 
            # can get the `equivalent_at_position` in the aligned hit sequence
            # simply using the index in the aligned sequence now.
    equivalent_at_position = aligned_hit_seq[i]
    # Collect details while have the position information for query and hit 
    # sequences
    # For sequences near the start, the subtraction of the 
    # `num_flanking_to_get_for_details` may result in negative numbers which 
    # cause back-counting when used in a slice and  so adjust for that. At the 
    # right-side end, the value added to the other may just cause a number
    # higher than the length if close to the far end, but Python just gives 
    # through to the end of the string if the second number in the slide is 
    # larger & so it's not an for that side of the strings.
    # Then define the same parts to make up the string on each side of target 
    # residue.
    nplus = num_flanking_to_get_for_details
    q_ind = "q " # indicator to put at left side;make sure same len as h_ind
    h_ind = "h " # indicator to put at left side;make sure same len as g_ind
    if (i-nplus < 0):
        left_side_start = 0
    else:
        left_side_start = i-nplus
    right_side_end = i+nplus
    q_first_part = aligned_query_seq[left_side_start:i]
    q_d = q_first_part + aligned_query_seq[i:right_side_end+1] #plus one because
    # I want to include all the way through the index at right side of the colon
    # and they way Python normally works with slices is to mean up to but NOT
    # including the index at the right side of the colon.
    h_d = (aligned_hit_seq[left_side_start:i] + 
        aligned_hit_seq[i:right_side_end+1])
    id_spacer = len(q_ind)* " "
    indicator_line = id_spacer+ len(q_first_part)*" " + "|"
    if equivalent_at_position == "-":
        interacting_res_nums_without_equivs.append(res_num)
        # collect details while have the position information for query and hit 
        # sequences
        dictionary_to_add_to = interacting_res_nums_without_equivs_details
    else:
        interacting_res_nums_with_equivs.append(res_num)
        # collect details while have the position information for query and hit 
        # sequences
        dictionary_to_add_to = interacting_res_nums_with_equivs_details
    dictionary_to_add_to[res_num] = {
        "residue at pos": aligned_query_seq[i],
        "equivalent_at_position":equivalent_at_position,
        "seq_details": indicator_line + "\n" + q_ind + q_d + "\n" + h_ind + h_d
                                    }
percent_without_equivalents = len(
    interacting_res_nums_without_equivs)/float(
    len(interacting_res_nums_with_data))
percent_with_equivalents = len(
    interacting_res_nums_with_equivs)/float(len(interacting_res_nums_with_data))
assert round(percent_without_equivalents + percent_with_equivalents) ==1,("The "
    "percent_without_equivalents and percent_with_equivalents for the residues "
    "with information in the current alignment should sum to 100% and they "
    "don't. ERROR?")


# CATEGORGIZE CONSERVATION FOR EACH INTERACTING RESIDUE WITH AN EQUIVALENT
#------------------------------------------------------------------------------#
# If there are any interacting residues with equivalents, make a note of those
# that don't change and those that do. Further categorize those that do change
# into 'strongly similar' and 'weakly similar' based on similar approach as 
# used in script `calculate_cons_for_clustal_protein.py`, which bases those 
# calls on
# https://www.ebi.ac.uk/seqdb/confluence/display/JDSAT/Bioinformatics+Tools+FAQ#BioinformaticsToolsFAQ-WhatdoconsensussymbolsrepresentinaMultipleSequenceAlignment?
identical_equivalents = []
different_equivs_strongly_conserved = []
different_equivs_weakly_conserved = []
different_equivs_not_conserved = []
# Note that for the assignment of the 'conservation call' dictionary entry,
# I used Black to reformat an extracted version of this script to see what 
# should be done since the assignment was going beyond 70 characters per line. 
# (Oddly, I put in a test of a a single long text line and it didn't fix that, 
# although it did move the parantheses to other lines to shorten it some. I 
# thought Black was like a miracle worker the way people talk about it.) After,
# I got the pattern, I made it by hand in Sublime because Sumblime didn't seem 
# to like doing it itself.
if interacting_res_nums_with_equivs_details:
    for res in interacting_res_nums_with_equivs_details:
        current_q_res = (
            interacting_res_nums_with_equivs_details[res]['residue at pos'])
        current_h_res = (interacting_res_nums_with_equivs_details[res]
            ['equivalent_at_position'])
        current_pos_list = [current_q_res,current_h_res]
        if all(ele == current_pos_list[0] for ele in current_pos_list):
            identical_equivalents.append(res)
            interacting_res_nums_with_equivs_details[res][
                "conservation_call"
            ] = "identical"
        elif any([set(current_pos_list).issubset(
            x) for x in strongly_similar_aa_tuples]):
            different_equivs_strongly_conserved.append(res)
            interacting_res_nums_with_equivs_details[res][
                "conservation_call"
            ] = "strongly conserved"
        elif any([set(current_pos_list).issubset(
            x) for x in weakly_similar_aa_tuples]):
            different_equivs_weakly_conserved.append(res)
            interacting_res_nums_with_equivs_details[res][
                'conservation_call'] = 'weakly conserved'
        else:
            different_equivs_not_conserved.append(res)
            interacting_res_nums_with_equivs_details[res][
                'conservation_call'] = 'not similar'

assert len(interacting_res_nums_with_equivs_details) == len(
    interacting_res_nums_with_equivs)
percent_with_equivalents_identical = len(identical_equivalents)/float(
    len(interacting_res_nums_with_equivs))
percent_with_equivalents_strongly_conserved = len(
    different_equivs_strongly_conserved)/float(
    len(interacting_res_nums_with_equivs))
percent_with_equivalents_weakly_conserved = len(
    different_equivs_weakly_conserved)/float(
    len(interacting_res_nums_with_equivs))
percent_with_equivalents_not_conserved = len(
    different_equivs_not_conserved)/float(
    len(interacting_res_nums_with_equivs))
'''
for k in interacting_res_nums_with_equivs_details:
    print(interacting_res_nums_with_equivs_details[k])
    print(interacting_res_nums_with_equivs_details[k]['seq_details'])
'''

# SUMMARY REPORT
#------------------------------------------------------------------------------#
separator_string = ("-"*80)
paratheses_note = ("\nThe number in parentheses following the percent is the "
    "specific # of residues.")
sys.stderr.write("\nDetermination of EQUIVALENTS Completed.\n\n"
    "*************************************RESULTS****************************"
    "********")
sys.stderr.write("\nA total of {} residues of chain {} interact with chain {} "
    "in the structure\nwith PDB id of {}.".format(
    len(interacting_res_nums), structure_chain1, structure_chain2, structure))
if interacting_res_nums_without_data:
    sys.stderr.write("\n\n*****NOTE*********NOTE*********NOTE*********NOTE*****"
        "*******")
    sys.stderr.write("\nNo information is provided in the aligned sequences "
    "for the following chain {} residues that interact\nwith chain {}, and "
    "therefore without further work in the form of extending the alignment, "
    "cannot say anything further about "
    "these:".format(structure_chain1,structure_chain2))
    sys.stderr.write("\n{}".format(interacting_res_nums_without_data.join(",")))
    sys.stderr.write("\n\n*END OF NOTE***END OF NOTE**END OF NOTE***END OF NOTE"
        "*******\n")
sys.stderr.write("\nThe following percent of chain {} residues interacting "
    "with chain {} residues\nhave equivalents in the other aligned "
    "sequence:".format(structure_chain1,structure_chain2))
sys.stderr.write("\n{:.1%} ({})".format(percent_with_equivalents, len(
    interacting_res_nums_with_equivs)))
sys.stderr.write(paratheses_note)
if percent_without_equivalents > percent_without_equivalents:
    sys.stderr.write("\nThe following percent of chain {} residues interacting "
        "with chain {} residues lack equivalents in the other aligned "
        "sequence:".format(structure_chain1,structure_chain2))
    sys.stderr.write("\n{:.1%} ({})".format(percent_without_equivalents, len(
        interacting_res_nums_without_equivs)))
    sys.stderr.write(paratheses_note)
sys.stderr.write("\n{}".format(separator_string))
sys.stderr.write("\nPercent breakdown of the conservation categories"
    " for the interacting\nresidues of chain {} that have equivalents in\n"
    "the other aligned sequence:".format(
    structure_chain1))
sys.stderr.write("\nidentical: {:.1%} ({})".format(
    percent_with_equivalents_identical, len(identical_equivalents)))
sys.stderr.write("\nstrongly conserved: {:.1%} ({})".format(
    percent_with_equivalents_strongly_conserved, 
    len(different_equivs_strongly_conserved)))
sys.stderr.write("\nweakly conserved: {:.1%} ({})".format(
    percent_with_equivalents_weakly_conserved,
    len(different_equivs_weakly_conserved)))
sys.stderr.write("\nnot conserved: {:.1%} ({})".format(
    percent_with_equivalents_not_conserved,len(different_equivs_not_conserved)))
sys.stderr.write(paratheses_note)
sys.stderr.write("\n{}".format(separator_string))
sys.stderr.write("\nThe specific residue positions of the conservation "
    "categories for the interacting\nresidues of chain {} that have "
    "equivalents:".format(structure_chain1))
sys.stderr.write("\nidentical: {}".format(
    ", ".join([str(i) for i in sorted(identical_equivalents)])))
sys.stderr.write("\nstrongly conserved: {}".format(
    ", ".join([str(i) for i in sorted(different_equivs_strongly_conserved)])))
sys.stderr.write("\nweakly conserved: {}".format(
    ", ".join([str(i) for i in sorted(different_equivs_weakly_conserved)])))
sys.stderr.write("\nnot conserved: {}".format(
    ", ".join([str(i) for i in sorted(different_equivs_not_conserved)])))
sys.stderr.write("\n{}".format(separator_string))



# Determine if notebook had provided settings for whether to display detailed
# sequence context for residues with and without equivalents. If they, weren't
# set there, use defaults set in this script under the section
# 'SETTINGS FOR THE MAIN PART OF THE SCRIPT' above.
try:
    type(details_for_those_with_equivalents) # triggers check if defined
except NameError:
    # if not defined, set as related default hard-coded in script
    details_for_those_with_equivalents = (
        report_details_for_those_with_equivalents)
try:
    type(details_for_those_without_equivalents) # triggers check if defined
except NameError:
    # if not defined, set as related default hard-coded in script
    details_for_those_without_equivalents = (
        report_details_for_those_without_equivalents)

if details_for_those_with_equivalents:
    sys.stderr.write("\n\nDETAILS for chain {} residues interacting with chain "
    "{} residues\nthat have equivalents in the other aligned "
    "sequence:".format(structure_chain1,structure_chain2))
    for res,details in sorted(interacting_res_nums_with_equivs_details.items()):
        sys.stderr.write("\n\nresidue #{} of chain {}:".format(
            res,structure_chain1))
        sys.stderr.write("\n{}".format(details["seq_details"]))
        sys.stderr.write("\n\nquery aa:{}; hit aa: {}; "
            "\nconservation class: {}\n".format(
            details["residue at pos"],
            details["equivalent_at_position"],
            details["conservation_call"]))
    sys.stderr.write("\n{}".format(separator_string))

if details_for_those_without_equivalents:
    sys.stderr.write("\n\nDETAILS for chain {} residues interacting with chain "
    "{} residues\nthat LACK equivalents in the other aligned "
    "sequence:".format(structure_chain1,structure_chain2))
    for res,details in sorted(interacting_res_nums_without_equivs_details):
        sys.stderr.write("\n\nresidue #{} of chain {}:".format(
            res,structure_chain1))
        sys.stderr.write("\n{}".format(details["seq_details"]))
        sys.stderr.write("\n\nquery aa:{}".format(details["residue at pos"]))
    sys.stderr.write("\n{}".format(separator_string))



#######------------------END OF MAIN SECTION------------------------------######
################################################################################
