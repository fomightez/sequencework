#!/usr/bin/env python
# check_for_omega_HEG.py
__author__ = "Wayne Decatur" #fomightez on GitHub
__license__ = "MIT"
__version__ = "0.1.0"


# check_for_omega_HEG.py by Wayne Decatur
# ver 0.1
#
#*******************************************************************************
# Verified compatible with both Python 2.7 and Python 3.7; written initially in 
# Python 3. 
#
#
# PURPOSE: Takes a sequence file (FASTA format) containing a fungal mitochondrial 
# large ribosomal subunit RNA and checks if it contains what seems to be a match 
# to the Homing Endonuclease Gene (HEG) of the  omega intron of S. cerevisiae. 
# Doesn't make a call about presence of omega intron per se. And so only use
# after positive via `check_for_omega_intron.py`. Better yet, best to use after 
# `find_mito_fungal_lsu_rRNA_and_check_for_omega_intron.py` because you can pass
# the seqeunce id of the 21S containing fragment as well to make detection of
# the homing endonuclease better.
# Because the script uses S. cerevisiae S288C 
# mitochondrial large ribosomal subunit RNA for this effort, I say it is for 
# fungal mitochondrial large ribosomal subunit RNA, but depending on 
# conservation, it may be limited to budding yeast or sub-classes of budding 
# yeast. It should help determine that by using with the Shen et al 2018 data 
# for genomes of 332 budding yeast.
#
# The code is based on heavily `check_for_omega_intron.py`.
# 
# Script requires BLAST+ to be already intalled in environment. An option is
# available at https://github.com/fomightez/blast-binder; go there and click
# `launch binder`, upload this script, and then use the terminal or notebook to
# execute the script.
# 
#
#
#
# Related:
# - check_for_omega_intron.py
# - fix_lsu_rRNA_annotation_in_gff_resulting_from_mfannot.py
# - find_mito_fungal_lsu_rRNA_and_check_for_omega.py
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
# TO RUN:
# Examples,
# Enter on the command line of your terminal, the line
#-----------------------------------
# python check_for_omega_HEG.py mito_seq.fa
#-----------------------------------
# Issue `check_for_omega_HEG.py -h` for 
# details.
# 
#
#
#
# To use this after importing/pasting or loading into a cell in a Jupyter 
# notebook, specify at least the file of annotations:
# from check_for_omega_HEG import check_for_omega_HEG
# check_for_omega_HEG("mito_seq.fa")
#
# 
#
'''
CURRENT ACTUAL CODE FOR RUNNING/TESTING IN A NOTEBOOK WHEN IMPORTED/LOADED OR 
PASTED IN ANOTHER CELL:
from check_for_omega_HEG import check_for_omega_HEG
check_for_omega_HEG("mito_seq.fa")
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












# Downloaded using `GSD get omega intron and omega HEG sequences.ipynb`, see it
# at .
# 1143 bp; HEG ORF is 708 bps
cer_omega  = '''>OMEGA_INTRON chrMT coordinates 60725 to 61867
ATTTACCCCCTTGTCCCATTATATTGAAAAATATAATTATTCAATTAATTATTTAATTGA
AGTAAATTGGGTGAATTGCTTAGATATCCATATAGATAAAAATAATGGACAATAAGCAGC
GAAGCTTATAACAACTTTCATATATGTATATATACGGTTATAAGAACGTTCAACGACTAG
ATGATGAGTGGAGTTAACAATAATTCATCCACGAGCGCCCAATGTCGAATAAATAAAATA
TTAAATAAATATCAAAGGATATATAAAGATTTTTAATAAATCAAAAAATAAAATAAAATG
AAAAATATTAAAAAAAATCAAGTAATAAATTTAGGACCTAATTCTAAATTATTAAAAGAA
TATAAATCACAATTAATTGAATTAAATATTGAACAATTTGAAGCAGGTATTGGTTTAATT
TTAGGAGATGCTTATATTCGTAGTCGTGATGAAGGTAAACTATATTGTATGCAATTTGAG
TGAAAAAATAAGGCATACATGGATCATGTATGTTTATTATATGATCAATGAGTATTATCA
CCTCCTCATAAAAAAGAAAGAGTTAATCATTTAGGTAATTTAGTAATTACCTGAGGAGCT
CAAACTTTTAAACATCAAGCTTTTAATAAATTAGCTAACTTATTTATTGTAAATAATAAA
AAACTTATTCCTAATAATTTAGTTGAAAATTATTTAACACCTATAAGTTTAGCATATTGA
TTTATAGATGATGGAGGTAAATGAGATTATAATAAAAATTCTCTTAATAAAAGTATTGTA
TTAAATACACAAAGTTTTACTTTTGAAGAAGTAGAATATTTAGTTAAAGGTTTAAGAAAT
AAATTTCAATTAAATTGTTATGTTAAAATTAATAAAAATAAACCAATTATTTATATTGAT
TCTATAAGTTATTTAATTTTTTATAATTTAATTAAACCTTATTTAATTCCTCAAATGATA
TATAAATTACCTAATACTATTTCATCCGAAACTTTTTTAAAATAATATTCTTATTTTTAT
TTTATGATATATTTCATAAATATTTATTTATATTAAATTTTATTTGATAATGATATAGTC
TGAACAATATAGTAATATATTGAAGTAATTATTTAAATGTAATTACGATAACAAAAAATT
TGA
'''

# VERSION WITHOUT HOMING ENDONUCLEASE GENE
# predicted to be ~336 bps
cer_omega_sansHEG  = '''>OMEGAminusHEG S.cerevisiae omega with homing endonuclease gene and associated DNA removedARTIFICAL
ATTTACCCCCTTGTCCCATTATATTGAAAAATATAATTATTCAATTAATTATTTAATTGA
AGTAAATTGGGTGAATTGCTTAGATATCCATATAGATAAAAATAATGGACAATAAGCAGC
GAAGCTTATAACAACTTTCATATATGTATATATACGGTTATAAGAACGTTCAACGACTAG
ATGATGAGTGGAGTTAACAATAATTCATCCACGAGCGCCCAATGTCGAATAAATAAAATA
TTAAATAAATTTTTATTTGATAATGATATAGTCTGAACAATATAGTAATATATTGAAGTA
ATTATTTAAATGTAATTACGATAACAAAAAATTTGA
'''
#----------------------------------------------------------------------------
# IDEA IF CURRENT APPROACH SEEMS NOT PERFOMRING WELL:
# OR SHOULD THE VERSION WITHOUT HOMING ENDONUCLEASE GENE TO COMPARE JUST BE THE 
#FIRST 250 of the omega intron???
#----------------------------------------------------------------------------
def stringFASTA2seq(s):
    '''
    Takes a FASTA file contents that are currently as a Python string and 
    converts it to the sequence string. In other words, it removes the 
    description line and the line breaks and just makes a sequence string.
    '''
    l = s.split("\n")
    return "".join(l[1:])
length_cer_omega_sansHEG = len(stringFASTA2seq(cer_omega_sansHEG)) #=336 #<--
# speculative length in bp of S. cerevisiae
# omega intron  without HEG. This is based on paradoxus N44 and CBS432 omega 
# intron sequences which naturally lack the homing endonuclease gene and some 
# flanking DNA. S288C omega intron size with HEG is 1143 bp (Cannot just use 
# length of the `cer_omega_sansHEG` because not read in yet as FASTA and so 
# string includes the description line and line endings.) Could be made better 
# if I find actual S. cerevisiae strains that have omega lacking HEG; however, 
# there are conflicting figure contents in Wu and Hao, 2014 
# (https://www.ncbi.nlm.nih.gov/pubmed/24515269) as whether that exists. 
# Figure 1 suggest not existent. Yet, Figure 2 suggests maybe some have the HEG? 
# And so will be interesting to see what find among the 1011 cerevisiae.
near_junction = 221 #point in S.cerevisiae near the junction where versions of
# intron without the endonuclease start to have sequences associated with 
# upstream of homing endonuclease gene (HEG) 

#*******************************************************************************
#*******************************************************************************
###DO NOT EDIT BELOW HERE - ENTER VALUES ABOVE###

import sys
import os
import subprocess
import pandas as pd


###---------------------------HELPER FUNCTIONS-------------------------------###

    
def determine_omega_HEG_presence(seq_file, 
    mito_LSU_containing_id="NO_seq_PROVIDED", df = None, bitscore_cutoff = 99):
    '''
    Takes the following:
    `seq_file`, name of sequence file
    `df`, optional dataframe resulting from BLAST query of the coding sequence

    Returns:
    True or false

    Takes a sequence file name, and optionally a dataframe resulting from a 
    BLAST query of the coding sequence of S. cerevisiae S228C 
    omega introns, and determines whether the omega intron Homing Endonuclease 
    Gene (HEG) seems present or not.

    Optionally, you can also provide the identitifer specifying the specific 
    sequence entry that has the 21S rRNA-containing sequence. This will improve
    the call when providing a multi-FASTA file from a genome assembly. The 
    identifier will correspond most often to 'entry_id' produced previously by 
    my script `find_mito_fungal_lsu_rRNA_and_check_for_omega_intron.py`

    The ability to take a dataframe is meant for when using this script imported
    into another where related queries already have taken place and makes no
    sense to repeat.
    '''
    if df is None:
        # Need to do BLAST query with GUESS at HEG-less omega intron of S. 
        # cerevisiae S228C first to make the dataframe other scripts may already 
        # have made.
        #
        # Need to store 'cer_omega_sansHEG' and 'cer_omega' as files so BLAST 
        # can use them.
        cer_omega_fn = "cer_omega_sansHEG.fa"
        with open(cer_omega_fn, "w") as q_file:
            q_file.write(cer_omega_sansHEG)
        sys.stderr.write("Checking hits against just HEG-less portion of itron"
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
            'sseq" -task blastn'.format(cer_omega_fn, seq_file))
        result = subprocess.check_output(cmd, shell=True) # based on 
            # https://stackoverflow.com/a/18739828/8508004
        result = result.decode("utf-8") # so it isn't bytes
        from blast_to_df import blast_to_df
        df = blast_to_df(result)
        df = df[df.bitscore > bitscore_cutoff]
    

    # Determine if it contains omega intron HEG. Doing this by now checking with 
    # the full intron sequence of cerevisiae omega intron, exlcuding matches 
    # that start above a site near the 5'-exon-into junction, and then examining 
    # if the hits containing that site increase substantially in size. For 
    # defining a substantial size increase I am arbitrarily using 6% of 
    # the HEG size (807 bp;includes endonuclease plus a little flanking) or 
    # `length(cer_omega SEQUENCE)-length(cer_omega_sansHEG SEQUENCE) * 0.06` ). 
    # Using 5'-end because that part is the largest segment that matches other 
    # omega introns even if homing endonuclease is absent. The part after 
    # I-SceI is much smaller, and thus harder to detect well.
    omegaHEG_present = False
    cer_omega_fn = "cer_omega.fa"
    with open(cer_omega_fn, "w") as q_file:
        q_file.write(cer_omega)
    sys.stderr.write("\nChecking for omega intron HEG presence"
        "...\n")
    cmd="makeblastdb -in {} -dbtype nucl".format(seq_file)
    subprocess.run(cmd, shell=True) # based on 
        # https://docs.python.org/3/library/subprocess.html#using-the-subprocess-module
        # and https://stackoverflow.com/a/18739828/8508004
    cmd = ('blastn -query {} -db {} -outfmt "6 qseqid sseqid '
        'stitle pident qcovs length mismatch gapopen qstart qend sstart '
        'send qframe sframe frames evalue bitscore qseq '
        'sseq" -task blastn'.format(cer_omega_fn, seq_file))
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
    # if `mito_LSU_containing_id` provided, restrict to it.
    if mito_LSU_containing_id != "NO_seq_PROVIDED":
        blast_df_forfull = (
            blast_df_forfull[blast_df_forfull.sseqid == mito_LSU_containing_id])
        df = (df[df.sseqid == mito_LSU_containing_id])
    else:
        # In order to avoid issues when checking sequences with lots of sequences 
        # that may cause spurious signals, restrict dataframes to contig/sequences, 
        # which are specified by `sseqid`, of best scoring matches.  
        best_scoring_sequence_id = df.loc[df.bitscore.idxmax()].sseqid # best match
        # to the omega intron sequence without the HEG will specify 
        # 21S-containing id. 
    # Before restricting dataframes to sequences spanning in front of the HEG
    # 5'-junction, store rawer dfs can use them in case no sequences span the 
    # HEG 5'-junction and need to look at the hits to the HEG sequences more.
    unsubset_blast_df_forfull = blast_df_forfull.copy()
    unsubset_df = df.copy()
    # Now to remove matches that start above the `near_junction` site because
    # I want the information for the matches that would span into what would
    # be the 5'-end of the intron, if intron present.
    blast_df_forfull =blast_df_forfull[blast_df_forfull.qstart < near_junction]
    df = df[df.qstart < near_junction]
    # because I chose the segments that begin less than but closest to the 
    # potential 5'-junction to check, it is possible the segments I am checking
    # don't span into the intron actually and I tried to build some testing
    # for that into below.
    row_of_interest_for_full = abs(
        blast_df_forfull['qstart'] - near_junction).idxmin()
    row_of_interest_for_HEGless = abs(df['qstart'] - near_junction).idxmin()
    omega_HEG_size = len(stringFASTA2seq(cer_omega))-length_cer_omega_sansHEG
    substantial_increase_in_matches_cutoff = omega_HEG_size * 0.06
    # While running the related script,`check_for_omega_intron.py`, I noted one
    # case where I wasn't getting a hit reasonably-scoring hit against genomic 
    # DNA for the region near the intron insertion site but was getting a good 
    # hit spanning the intron insertion site against the non-coding sequence.
    # Thus in this odd case the row examined for the results agains genomic 
    # weren't look at intron insertion site at all & the result against codring
    # was stringly indicating it had no intron. But because of the mismatch 
    # of what was being considered, the sequence was erroneously showing a 
    # size change. So I should probably add eliminating possibility of mismatch
    # about what is being compared and whether match to HEG-less seqeunce
    # is strong and no reasonably scoring data offered for region where HEG is 
    # inserted when queried against the S288C version with HEG intact.
    full_match_spans_near_junction_val = (
        blast_df_forfull.loc[row_of_interest_for_full].qend > near_junction)
    HEGless_match_spans_near_junction_val = (
        df.loc[row_of_interest_for_HEGless].qend > near_junction)
    if HEGless_match_spans_near_junction_val and (
        not full_match_spans_near_junction_val):
        omegaHEG_present = False
    # I was seeing some rare cases (CRI_2 of the 1011 collection so far; plus
    # comparable cases for`check_for_omega_intron.py` [(YDJ, BSR_4, & GAT].) 
    # where no well-scoring hits spanning junction site but very well 
    # scoring hits that include segment shortly following the junction. This is 
    # to deal with them by looking for cases where there are no spanning hits 
    # but where if you apply a more stringent bitscore cutoff, the size of hits 
    # against intact aren't mirrored by increase in other set. And the size of 
    # the extra hit(s) is substantial. (This is a variation of my simplisitc 
    # previous approach for detecting insertions where I looked for number of 
    # fragments to change.)
    elif (not HEGless_match_spans_near_junction_val) and (
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
            omegaHEG_present  = True
    # the primary check for omega HEG presence/absence state
    elif (blast_df_forfull.loc[row_of_interest_for_full].length >= 
        df.loc[row_of_interest_for_HEGless].length + 
        substantial_increase_in_matches_cutoff):
        # ADD CHECK FOR A SUBSTANTIAL ORF IN CANDIDATE REGION?Could be difficult
        # since right now I just have BLAST hit details.
        omegaHEG_present = True


    return omegaHEG_present



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

def check_for_omega_HEG(seq_file, mito_LSU_containing_id="NO_seq_PROVIDED", 
    df = None, bitscore_cutoff = 99):
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

    Optionally, you can also provide the identitifer specifying the specific 
    sequence entry that has the 21S rRNA-containing sequence. This will improve
    the call when providing a multi-FASTA file from a genome assembly. The 
    identifier will correspond most often to 'entry_id' produced previously by 
    my script `find_mito_fungal_lsu_rRNA_and_check_for_omega_intron.py`

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



    # Check for omega intron HEG in provided sequence
    #---------------------------------------------------------------------
    #
    omegaHEG_present = determine_omega_HEG_presence(seq_file,
        mito_LSU_containing_id,df,bitscore_cutoff)



    # Report on omega HEG presense status
    #---------------------------------------------------------------------
    # Report on omega HEG via feedback if called from command line
    if omegaHEG_present and __name__ == "__main__":
        sys.stderr.write("\n**Seems omega intron HEG is present. Further "
            "investigation is suggested.**")
    elif (not omegaHEG_present) and __name__ == "__main__":
        sys.stderr.write("\n**Simple examination for presence of omega "
            "intron HEG didn't indicate\nobvious presence. Further "
            "investigation is suggested.**")


    if omegaHEG_present:
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
    kwargs['mito_LSU_containing_id'] = mito_LSU_containing_id
    #kwargs['return_df'] = False #probably don't want dataframe returned if 
    # calling script from command line
    check_for_omega_HEG(seq_file,**kwargs)
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
    parser = argparse.ArgumentParser(prog='check_for_omega_HEG.py',
        description="check_for_omega_HEG.py \
        takes the file name of a sequence containing mitochondrial large \
        ribosomal subunit sequence and returns if omega intron seems present \
        or not.\
        **** Script by Wayne Decatur   \
        (fomightez @ github) ***")

    parser.add_argument("sequence_file", help="Name of file containing the \
        fungal mitochondrial large ribosomal subunit RNA sequence (FASTA \
        format).\
        ", metavar="SEQ_FILE")
    parser.add_argument("mito_LSU_id", nargs='?', help="**OPTIONAL**Specify \
        the sequence identifier of the 21S-containing sequence if the sequence \
        file supplied is a multi-FASTA sequence. Specifying, the 21S-containing \
        will improve the detection call. The identifier will correspond most \
        often to 'entry_id' produced previously by \
        `find_mito_fungal_lsu_rRNA_and_check_for_omega_intron.py`.", 
        default="NO_seq_PROVIDED", metavar="specific_id")
    # See
    # https://stackoverflow.com/questions/4480075/argparse-optional-positional-arguments 
    # and 
    # https://docs.python.org/2/library/argparse.html#nargs for use of `nargs='?'` 
    # to make input and output file names optional. Note that the square brackets
    # shown in the usage out signify optional according to 
    # https://stackoverflow.com/questions/4480075/argparse-optional-positional-arguments#comment40460395_4480202
    # , but because placed under positional I added clarifying text to help 
    # description.





    #I would also like trigger help to display if no arguments provided because 
    # need at least one for url
    if len(sys.argv)==1:    #from http://stackoverflow.com/questions/4042452/display-help-message-with-python-argparse-when-script-is-called-without-any-argu
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    seq_file = args.sequence_file
    mito_LSU_containing_id=args.mito_LSU_id


    main()

#*******************************************************************************
###-***********************END MAIN PORTION OF SCRIPT***********************-###
#*******************************************************************************