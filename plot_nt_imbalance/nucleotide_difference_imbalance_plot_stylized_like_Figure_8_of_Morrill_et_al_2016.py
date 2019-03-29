#!/usr/bin/env python
# nucleotide_difference_imbalance_plot_stylized_like_Figure_8_of_Morrill_et_al_2016.py
__author__ = "Wayne Decatur" #fomightez on GitHub
__license__ = "MIT"
__version__ = "0.1.0"


# nucleotide_difference_imbalance_plot_stylized_like_Figure_8_of_Morrill_et_al_2016.py by Wayne Decatur
# ver 0.2
#
#*******************************************************************************
# Verified compatible with both Python 2.7 and Python 3.7; written initially in 
# Python 3. 
#
# PURPOSE: Make a plot figure like Figure 8, panel B of Morrill et al 2016 
# (PMID: 27026700), see  
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4882425/figure/F8/ .
# Legend says:
'''
>"GC imbalance marks the CTD coding region and contributes to secondary structure formation. 
A, nucleotide statistics for the core of RPB1 and the CTD region. The CTD coding strand has reduced guanine bases and a high concentration of cytosine (gray).  
B, plot of nucleotide differences for a 20-kilobase region of yeast chromosome IV including RPB1. The data were calculated from 150-base pair windows with 50 base pairs of overlap between adjacent windows and plotted in Excel. The CTD region has a very negative value of G-C. A-T for this same region is less variable."
'''
#
# Note that this script was originally intended to be used after using my script 
# `get_chromosomal_coordinates_as_FASTA.py` from 
# https://github.com/fomightez/yeastmine. Of course, by default the files 
# generated from there have the pertinent infromation about the specific 
# chromosomes and the coordinates spanned both as part of the name of the file 
# and in the description line of the sequence file itself. However, in order to 
# keep it general for using with sequences sourced from elsewhere, I don't 
# assume that as source and instead ask for some of the infromation for what to 
# use in labels. (<--SUBJECT TO CHANGE. Maybe add a flag I can invoke to declare
# source is my `get_chromosomal_coordinates_as_FASTA.py` scipt?)
#
#
# Dependencies beyond the mostly standard libraries/modules:
#  
#
# 
#
#
# VERSION HISTORY:
# v.0.1. basic working version
#
#
# to do:
# - fix 'Purpose' header section above
# - document at github
#
#
# Note for myself:
# For impetus and development see `
# Making a line plot like GC imbalance in Figure 8 of Morrill et al 2016.md` and
# `developing_plot_of_nucelotide_difference_imbalance.ipynb`.
#
#
#
#
# TO RUN:
# Example,
# Enter on the command line of your terminal, the line
#-----------------------------------
# python nucleotide_difference_imbalance_plot_stylized_like_Figure_8_of_Morrill_et_al_2016.py chrV_2000-22000.fa ChrV 2000
#-----------------------------------
#
#
#
# To use this after importing/pasting or loading into a cell in a Jupyter 
# notebook, specify necessary items, which are sequence file, chromosome or  
# contig or scaffold designation, and starting position in provided sequence 
# file, in call to command:
#from nucleotide_difference_imbalance_plot_stylized_like_Figure_8_of_Morrill_et_al_2016 import nucleotide_difference_imbalance_plot_stylized_like_Figure_8_of_Morrill_et_al_2016
#%matplotlib inline
#ax = nucleotide_difference_imbalance_plot_stylized_like_Figure_8_of_Morrill_et_al_2016("chrV_2000-22000.fa", "ChrV", "2000",return_fig=True);
# Using `return_fig=True` it will return a plot figure object, without it an 
# image file will be produced. From the example above, you can redisplay the 
# returned figure in another Jupyter cell later with `ax.figure` without running
# the code again.
#*******************************************************************************
#





#*******************************************************************************
##################################
#  USER ADJUSTABLE VALUES        #

##################################
#


chunk_size = 150 #for 150 bp window, based on Figure 8 of Morrill et al 2016 
overlap_specified = 50 #50 bp overlap of windows, based on Figure 8 of Morrill 
# et al 2016 (PMID: 27026700)
step_size = chunk_size - overlap_specified #100 for chunk_size of 150 & overlap 
# of 50; based on Figure 8 of Morrill et al 2016 (PMID: 27026700)

suffix_for_saving_plot = "_nt_imblance_plot.png"
default_plt_image_size = (13,8)
yaxis_label = "nucleotide imbalance"
xaxis_label_prefix = "position along "

#
#*******************************************************************************
#**********************END USER ADJUSTABLE VARIABLES****************************















#*******************************************************************************
#*******************************************************************************
###DO NOT EDIT BELOW HERE - ENTER VALUES ABOVE###

import sys
import os
import argparse
from pyfaidx import Fasta
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt




###---------------------------HELPER FUNCTIONS---------------------------------###



def generate_output_file_name(file_name, suffix):
    '''
    Takes a file name as an argument and returns string for the name of the
    output file. The generated name is based on the original file
    name.

    It also indicates in resulting file name name specific chromsomes or 
    scaffolds if plotting was limited to those.

    Specific example
    ================
    Calling function with
        ("data1.txt", "_across_chr.png")
    returns
        "data1_across_chr.png"

    '''
    main_part_of_name, file_extension = os.path.splitext(
        file_name) #from https://stackoverflow.com/a/541394/8508004
    return main_part_of_name + suffix

def calc_nt_diff(di_base,seq):
    '''
    Takes a string made of two bases along with a sequence 
    as a string. Calculates the nucleotide different 
    imbalance value for the base combination and returns that.
    It is insensitive to case.
    Nucleotide different imbalance value is meant to match
    Figure 8 of Morrill et al 2016 (PMID: 27026700). I found
    the description of the calculation and what was being 
    plotted rather vague, and so it took some testing of my
    speculations to work this out.
    
    '''
    assert len(di_base) == 2, ("The `calc_nt_diff` function takes "
    "a string of two bases, such as 'GC' or 'AT'.")
    first_base = di_base[0].lower()
    secnd_base = di_base[1].lower()
    return seq.lower().count(secnd_base) - seq.lower().count(first_base)

def gen_chunk_string_with_different_step(a_list, chunk_size, step_amount):
    """Yield successive n-sized chunks from list, stepping /stride by 
    step_amount."""
    for i in range(0, len(a_list), step_amount):
        yield a_list[i:i+chunk_size]

def midpoint(items):
    '''
    takes a iterable of items and returns the midpoint (integer) of the first 
    and second values
    '''
    return int((int(items[0])+int(items[1]))/2)




###--------------------------END OF HELPER FUNCTIONS-------------------------###
###--------------------------END OF HELPER FUNCTIONS-------------------------###






#*******************************************************************************
###------------------------'main' function of script--------------------------##

def nucleotide_difference_imbalance_plot_stylized_like_Figure_8_of_Morrill_et_al_2016(
    sequence_file, chr_, position_corresponding_to_first_nt, save_vg=False, 
    return_fig=False):
    '''
    Main function of script. 

    Make a plot figure like Figure 8, panel B of Morrill et al 2016 
    (PMID: 27026700), see  
    https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4882425/figure/F8/ .

    Inputs:
    - sequence file of spanning desired coordinates in FASTA format with 
    single sequence in it. (Or at least first sequence be the one to use.)
    - text string of designation to use for chromosome/contig/scaffold/region in
    plot labels
    - position corresponding to first bases in provided sequence as an integer.
    To be used to label positions along coordinates in x-axis of plot.

    You can set `return_fig` to return a plot figure object and would be useful
    if calling from a Jupyter notebook or IPython. If you assign it to `ax` in
    your notebook, you can redisplay it in another cell via `ax.figure`.

    Saving as vector graphics is also an option. It is not the default because
    file size can get large with a lot of points.

    Saves an image file of the plot, or optionally, returns a plot object.
    '''

    # Retrieve sequence from sequence file and break it up into the set chunks
    #--------------------------------------------------------------------------- 
    seq_entries = Fasta(sequence_file)
    seq = seq_entries[0] # assume first one is the one to be used
    chunks = (
        list(gen_chunk_string_with_different_step(seq,chunk_size,step_size)))
    #discard any chunks at end less than the size of the set window
    chunks = [x for x in chunks if len(x)== chunk_size]


    # Assign midpoint postions to each chunk based on chunk_size, start position
    # input, and provided sequence length.
    # ALSO, while going through chunks, might as well do calculation too:
    #Calculate G vs. C and A vs. T that I think might give results like Figure 8 
    # of Morrill et al 2016 (PMID: 27026700)
    #---------------------------------------------------------------------------
    #determining midpoints for each chunk (relative length of provided sequence) 
    # and assiging to a list. (used a dictionary in development but a list
    # should be better for python 2.7 compatibility)
    # Had calculation as separate step in development but why interate again.
    # (Doing similar thing where now storing `chunks_diffs` as list and 
    # not diction for  python 2.7 compatibility.)
    import sys
    #chunks_midpoints = {} #key will be index of chunk in chunks
    chunks_midpoints = []
    #chunks_diffs = {} #key will be index of chunk in chunks. value will be tuple 
    # with each diff as an item
    chunks_diffs = []
    for indx,chunk in enumerate(chunks):
        # handle first chunk without much fanfare because easy calculation and
        # doesn't depend on a previous one
        if indx == 0 and len(chunk)== chunk_size:
            #chunks_midpoints[indx] = chunk_size/2
            chunks_midpoints.append(chunk_size/2)
        elif len(chunk)== chunk_size:
            start_curr_chunk = step_size * indx
            end_curr_chunk = step_size * indx + chunk_size
            #chunks_midpoints[indx]= midpoint((start_curr_chunk,end_curr_chunk)) 
            chunks_midpoints.append(midpoint((start_curr_chunk,end_curr_chunk)))
        else:
            sys.stderr.write("\n\nError? Issue with size of sequence chunks "
                "not matching expected/\n") #shouldn't happen
            sys.exit(1)
        GC_diff = calc_nt_diff("GC",str(chunk)) # casting 'Sequence' object to 
        # string with `str(chunk)`
        AT_diff = calc_nt_diff("AT",str(chunk)) # casting 'Sequence' object to 
        # string with `str(chunk)`
        #chunks_diffs[indx] = (GC_diff,AT_diff)
        chunks_diffs.append((GC_diff,AT_diff))

    # Adjust potions of midpoints to account for first base in sequence file 
    # being not first base of chromosome
    #---------------------------------------------------------------------------
    # correct chunks_midpoints to take into account that first postion in 
    # provided sequence might not be first position of that sequence along 
    # chromosome so that labels for x-axis will match situation
    # (This was done with a dictionary comprehension in development where used
    # dictionary)
    start_pos = position_corresponding_to_first_nt
    #chunks_midpoints = {k:v+start_pos for k,v in chunks_midpoints .items()}
    chunks_midpoints = [x+start_pos for x in chunks_midpoints]




    '''Combined into first interation through list of chunks above
    #Calculate G vs. C and A vs. T that I think might give results like Figure 8 
    # of Morrill et al 2016 (PMID: 27026700)
    #---------------------------------------------------------------------------
    chunks_diffs = {} #key will be index of chunk in chunks. value will be tuple 
    # with each diff as an item
    for indx,seq in enumerate(chunks):
        GC_diff = calc_nt_diff("GC",str(seq)) # casting 'Sequence' object to 
        # string with `str(seq)`
        AT_diff = calc_nt_diff("AT",str(seq)) # casting 'Sequence' object to 
        # string with `str(seq)`
        chunks_diffs[indx] = (GC_diff,AT_diff)
    '''


    #Make the plot
    #---------------------------------------------------------------------------
    sns.set()
    plt.figure(figsize=default_plt_image_size)
    #indx = list(chunks_midpoints.values())
    #data = list(chunks_diffs.values())
    #indx and data were originally organized in dictionaries during development
    # but for better compatibility with Python 2.7 lists were used so order
    # maintained
    #df = pd.DataFrame(data, indx, ["GvsC", "AvsT"])
    df = pd.DataFrame(chunks_diffs, chunks_midpoints, ["GvsC", "AvsT"])
    ax = sns.lineplot(data=df)
    ax.set_ylabel(yaxis_label, fontsize = 16);
    ax.set_xlabel(xaxis_label_prefix+chr_, fontsize = 16);
    ax.legend(fontsize= 12);



    # Return plot figure object (meant for calling from Jupyter cell or IPython)
    # or save to file.
    #---------------------------------------------------------------------------
    if return_fig:
        sys.stderr.write("Plot figure object returned.")
        return ax
    else:
        #save image; standard for when called from command line; however, will
        # also be default for when called from Jupyter notebook / IPython unless 
        # called with `return_fig=True` 
        output_file_name = generate_output_file_name(sequence_file, 
            suffix_for_saving_plot)
        if save_vg:
            plt.savefig(output_file_name[:-4]+".svg", 
            orientation='landscape') # FOR VECTOR GRAPHICS; useful if merging 
            # into Adobe Illustrator. Based on 
            # https://neuroscience.telenczuk.pl/?p=331 ; I think ReportLab also 
            # outputs SVG?
            sys.stderr.write("\nPlot image saved to: {}\n".format(
                output_file_name[:-4]+".svg"))
        else:
            # save png
            plt.savefig(output_file_name)
            sys.stderr.write("\nPlot image saved to: {}\n".format(
                output_file_name))

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
    kwargs['return_fig'] = False #probably don't want plot object returned if 
    # calling script from command line & instead will trigger save of plot image
    #kwargs['return_df'] = False #probably don't want dataframe returned if 
    # calling script from command line
    kwargs['save_vg'] = False
    nucleotide_difference_imbalance_plot_stylized_like_Figure_8_of_Morrill_et_al_2016(
        sequence_file, chr_, position_corresponding_to_first_nt,**kwargs)
    # using 
    # https://www.saltycrane.com/blog/2008/01/how-to-use-args-and-kwargs-in-python/#calling-a-function
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
    parser = argparse.ArgumentParser(prog='\
    nucleotide_difference_imbalance_plot_stylized_like_Figure_8_of_Morrill_et_al_2016.py',
    description="\
    nucleotide_difference_imbalance_plot_stylized_like_Figure_8_of_Morrill_et_al_2016.py \
    makes a plot figure like Figure 8, panel B of Morrill et al 2016 \
    (PMID: 27026700), see  \
    https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4882425/figure/F8/.\
    **** Script by Wayne Decatur   \
    (fomightez @ github) ***")

    parser.add_argument("sequence_file", help="Name of sequence file to \
        use as input. Must be FASTA format. Can be a \
        multi-FASTA file, i.e., multiple sequences in FASTA format in one \
        file.", metavar="SEQUENCE_FILE")

    parser.add_argument("chr", type=str, help="Designation of chromosome / \
        contig/ scaffold / region corresponding to that in sequence file. To \
        be used for plot label. No spaces. Feel free to prefix with a species \
        indicator as well\
        .", metavar="CHR")

    parser.add_argument("position_corresponding_to_first_nt", type=int, help="\
        Coordinate position that corresponds to first position of specified \
        along chromosome. To be used for label of x-axis. No spaces.", 
        metavar="START_POS")
    parser.add_argument("-svg", "--save_vg",help=
    "add this flag to save as vector graphics \
    (**RECOMMENDED FOR PUBLICATION***) instead of default png. Not default or \
    saved alongside default because file size can get large due to the large \
    number of points.",
    action="store_true")



    #I would also like trigger help to display if no arguments provided because 
    #need at least one input file
    if len(sys.argv)==1:    #from https://stackoverflow.com/a/4042861/8508004
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    sequence_file = args.sequence_file
    chr_ = args.chr #cannot use `chr` because it is a built-in python function
    position_corresponding_to_first_nt = args.position_corresponding_to_first_nt
    save_vg = args.save_vg


    main()

#*******************************************************************************
###-***********************END MAIN PORTION OF SCRIPT***********************-###
#*******************************************************************************
