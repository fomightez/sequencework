#!/usr/bin/env python
# two_nucleotides_in_proximity_difference_imbalance_plot.py
__author__ = "Wayne Decatur" #fomightez on GitHub
__license__ = "MIT"
__version__ = "0.1.0"


# two_nucleotides_in_proximity_difference_imbalance_plot.py by Wayne Decatur
# ver 0.2
#
#*******************************************************************************
# Verified compatible with both Python 2.7 and Python 3.7; written initially in 
# Python 3. 
#
# PURPOSE: Make a plot similar to figure like Figure 8, panel B of Morrill et al 2016 
# (PMID: 27026700), see  
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4882425/figure/F8/ , but for 
# two different nucleotides occuring near each other and not just one. Was going
# to call it `dinucleotide_difference_imbalance_plot.py` but worried that 
# implied the two nucleotides had to be directly adjacent instead of just 
# nearby in sequence.
#
#
#
#
# Based on my script 
# `nucleotide_difference_imbalance_plot_stylized_like_Figure_8_of_Morrill_et_al_2016.py`
# that handled single nucleotide imbalances, spliitng the four into two sets of 
# two.
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
# python two_nucleotides_in_proximity_difference_imbalance_plot.py chrV_2000-22000.fa ChrV 2000
#-----------------------------------
#
#
#
# To use this after importing/pasting or loading into a cell in a Jupyter 
# notebook, specify necessary items, which are sequence file, chromosome or  
# contig or scaffold designation, and starting position in provided sequence 
# file, in call to command:
#from two_nucleotides_in_proximity_difference_imbalance_plot import two_nucleotides_in_proximity_difference_imbalance_plot
#%matplotlib inline
#ax = two_nucleotides_in_proximity_difference_imbalance_plot("chrV_2000-22000.fa", "ChrV", "2000",return_fig=True);
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

suffix_for_saving_plot = "_dibase_imblance_plot.png"
default_plt_image_size = (13,8)
yaxis_label = "dinucleotide imbalance"
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

def calc_dibase_diff(di_base_txt1,di_base_txt2,seq):
    '''
    Takes two strings both made of two bases along with a sequence 
    as a string. Calculates the nucleotide different 
    imbalance value for the base combination and returns that.
    It is insensitive to case.
    DiNucleotide difference imbalance value is meant to be similar to the 
    nucleotide difference imbalance value of
    Figure 8 of Morrill et al 2016 (PMID: 27026700), except here it will be 
    calculated with the sum of occurences of two nucleotides vs the sum of the
    occurences of the other two. Using `di` just two signal 'two'; they do not
    have to be adjacent, just within the same window, in other words in 
    proximity to each other.
    
    '''
    assert len(di_base_txt1) == 2 and len(di_base_txt2) == 2, ("The "
        "`calc_dint_diff` function takes two strings, both of two bases, such "
        "as 'GC' and 'AT'.")
    first_set_bases = list(di_base_txt1.lower())
    secnd_set_bases  = list(di_base_txt2.lower())
    tally_first_nts = seq.lower().count(first_set_bases[0]) + seq.lower(
        ).count(first_set_bases[1])
    tally_second_nts = seq.lower().count(secnd_set_bases[0]) + seq.lower(
        ).count(secnd_set_bases[1])
    return tally_first_nts - tally_second_nts

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
nt_set = {"G","C","T","A"}
def two_nucleotides_in_proximity_difference_imbalance_plot(
    sequence_file, chr_, position_corresponding_to_first_nt, dibase_text1="GC",
    chunk_size = chunk_size, overlap_specified = overlap_specified,
    save_vg=False, return_fig=False):
    '''
    Main function of script. 

    Make a plot figure like Figure 8, panel B of Morrill et al 2016 
    (PMID: 27026700), see  
    https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4882425/figure/F8/ , 
    except do it for the difference of two bases proximal to each other vs
    the other two.

    Inputs:
    - sequence file of spanning desired coordinates in FASTA format with 
    single sequence in it. (Or at least first sequence be the one to use.)
    - text string of designation to use for chromosome/contig/scaffold/region in
    plot labels
    - position corresponding to first bases in provided sequence as an integer.
    To be used to label positions along coordinates in x-axis of plot.
    - A string of two letters that represent two of the four bases `GATC`, and 
    so `GC`, `AT`, 'GA', 'GT', 'TC', 'AC', or any differently ordered 
    combination of letters is valid. Only ask for two because the remaining two
    are determined from the possibe four.

    Optionally you can provide integer settings for the `chunk_size` and 
    `overlap_specified` for the analyses windows. Without specifying them, by 
    default they are set to mirror Figure 8, panel B of Morrill et al 2016 
    (PMID: 27026700).

    You can set `return_fig` to return a plot figure object and would be useful
    if calling from a Jupyter notebook or IPython. If you assign it to `ax` in
    your notebook, you can redisplay it in another cell via `ax.figure`.

    Saving as vector graphics is also an option. It is not the default because
    file PNG-style image files are more familiar to most folks and more 
    convenient in Jupyter notebooks.

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
        dibase_set1 = set(list(dibase_text1))
        dibase_text2 = ('').join(nt_set.difference(dibase_set1))
        dibase_diff = calc_dibase_diff(dibase_text1,dibase_text2,str(chunk)) # 
        # casting 'Sequence' object to string with `str(chunk)`
        #chunks_diffs[indx] = (GC_vs_AT_diff)
        chunks_diffs.append(dibase_diff)

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
    #df = pd.DataFrame(data, indx, ["GCvsAT"])
    df = pd.DataFrame(chunks_diffs, chunks_midpoints, ["{} vs. {}".format(
        dibase_text1,dibase_text2)])
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
    kwargs['chunk_size'] = chunk_size
    kwargs['overlap_specified'] = overlap_specified
    kwargs['return_fig'] = False #probably don't want plot object returned if 
    # calling script from command line & instead will trigger save of plot image
    #kwargs['return_df'] = False #probably don't want dataframe returned if 
    # calling script from command line
    kwargs['save_vg'] = save_vg
    two_nucleotides_in_proximity_difference_imbalance_plot(
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
    two_nucleotides_in_proximity_difference_imbalance_plot.py',
    description="\
    two_nucleotides_in_proximity_difference_imbalance_plot.py \
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
    def valid_dibase(arg_string):
        if len(arg_string) == 2 and set(
            list(arg_string.upper())).issubset(nt_set):
            return arg_string
        else:
            msg = ("Not a valid set of letters representing two nucleotides: "
                "'{}'. Try something like `GC` or `AT`.".format(arg_string))
            raise argparse.ArgumentTypeError(arg_string)
    parser.add_argument("two_bases", type = valid_dibase, help="Two \
        letters representing one set of bases to check for imbalance vs. the \
        other two remaining out of the possible four from `GATC`. For example, \
        providing `GC` will result in checking for imbalance of nucleotides \
        `G` and `C` in close proximity vs. nucleotides `A` and `T`\
        .", metavar="TWO_BASES")
    parser.add_argument('-bl', '--block_size', action='store', type=int, 
    default= chunk_size, help="OPTIONAL: Use the `--block_size` flag followed \
    by an interger to provide a value to use as the span size (window of \
    basepairs) to analyze instead of the default of '{}'.".format(chunk_size))
    parser.add_argument('-ov', '--overlap_size', action='store', type=int, 
    default= overlap_specified, help="OPTIONAL: Use the `--overlap_size` \
    flag followed by an integer to specify the amount of overlap to use \
    between the \
    analysis windows instead of the default of '{}'.".format(overlap_specified))
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
    dibase_text1 = args.two_bases
    chunk_size=args.block_size
    overlap_specified = args.overlap_size
    save_vg = args.save_vg

    


    main()

#*******************************************************************************
###-***********************END MAIN PORTION OF SCRIPT***********************-###
#*******************************************************************************
