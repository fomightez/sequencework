#!/usr/bin/env python
# -*- coding: utf-8 -*-


# plot_5prime_end_from_bedgraph.py by Wayne Decatur
__author__ = "Wayne Decatur" #fomightez on GitHub
__version__ = "0.0.1"
__license__ = "MIT"



#*******************************************************************************
# Works in both Python 2.7 and Python 3.6. (****Verified****)
#
# PURPOSE: Work in progress script to take output in bedgraph form and plot
# number of ends per position for a region.
# The data from two sources needed to run this script: 1. bedgraph file for 
# 'forward' data., and 2. bedgraph-formatted file for 'reverse' data. 
# This data is then used to plot number of reads with starts at each position 
# to produce a bar chart.
# Need because want more flexible output of that which can be combined with 
# other data.
# It is assumed in this script, the provided bedgraphs only deal with one
# chromosome, specifically a region in that chromosome. And so just the 
# coordinates need designating.
# See http://genome.ucsc.edu/goldenPath/help/bedgraph.html about the bedgraph 
# format.
#
# 
# **"Starts" here means 5'-end of reads of that locate ribonucleotides, see 
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1521130 where it says,
# "Supplementary_files_format_and_content: bedGraph, count of Ribo-seq 5' ends 
# per position".
#
#
#
# This is based on `plot_coverage_and_starts.py` 
#
#
#
# Dependencies beyond the mostly standard libraries/modules:
#
#
#
# VERSION HISTORY:
# v.0.1. basic working version
#
#
#
#
#
# TO RUN:
# Example,
# Enter on the command line of your terminal, the line
#----------------------------------------------------
# python plot_5prime_end_from_bedgraph.py 220-516 forw_data.bedgraph rev_data.bedgraph
#-----------------------------------------------------
# 
# where the region to plot is provided first in the format `start-end`, 
# followed by the forward data, and the reverse data.
# Start and end are to be provided in normal terms where something starts with
# the first position being number 1. (I.e., no zero indexing)
# 
#works in Jupyter cell to call when script in home directory in same session:
# %run plot_5prime_end_from_bedgraph.py 4013-4359 forward.bedgraph reverse.bedgraph
# 

#
#
#*******************************************************************************
#


#*******************************************************************************
##################################
#  USER ADJUSTABLE VALUES        #

##################################
#
suffix_for_saving_result = "_ribonucleotide_plot.png"

title_prefix = "Region" # Change to `None` to suppress title

title_main = " HydEn-Seq-detected 5'-DNA ends"

plot_style = "seaborn" #try also `ggplot`,`default`, `bmh` or `grayscale`; use 
# `print(plt.style.available)` after appropriate imports to see others; 
# illustrated at https://matplotlib.org/examples/style_sheets/style_sheets_reference.html


#
#*******************************************************************************
#**********************END USER ADJUSTABLE VARIABLES****************************





























#*******************************************************************************
#*******************************************************************************
###DO NOT EDIT BELOW HERE - ENTER VALUES ABOVE###

import sys
import os
import argparse
import logging
import matplotlib # in order to use `matplotlib.use('Agg')`, need this first, see source of next line
matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend for running on PythonAnywhere. # from https://stackoverflow.com/questions/2801882/generating-a-png-with-matplotlib-when-display-is-undefined after searched error I was getting after upgrading matplolib in my Pythonanywhere account
import matplotlib.pyplot as plt
import seaborn as sns
#from statsmodels.nonparametric.smoothers_lowess import lowess
# from scipy.stats import mode

#DEBUG CONTROL
#comment line below to turn off debug print statements
#logging.basicConfig(stream=sys.stderr, level=logging.DEBUG)

###---------------------------HELPER FUNCTIONS---------------------------------###

def get_scores_from_bedgraph(bedgraph):
    '''
    Take a bedgraph file and returns a dictionary with positions as keys and
    associated score as the values.
    The first position is provided by the first position in the bedgraph file.
    Same for the last.
    Bedgraph files are to be zero-based and so switch to normal coordinates in 
    returned data.
    '''
    score_per_position_dict = {}
    # According to the format, only the specified positions will be graphed. For 
    # use here then any gap, later get a score value of zero. So no need to be 
    # be concerned here.
    with open(bedgraph, 'r') as input_:
        lines_processed = 0
        for line in input_:
            lines_processed += 1
            # according to the format, some text can be on the header to be skipped
            if line.startswith(
                "browser") or line.startswith("#") or line.startswith("track"):
                continue
            else:
                info = line.split()
                if lines_processed == 1:
                    graph_start = int(info[1])
                ''' # This section was from when I didn't realize easier to handle gaps later
                if lines_processed != 1 and int(info[1]) != end:
                    #print ("current 'start'",info[1]) # FOR DEBUGGING
                    #print ("current 'start' expected as",end) # FOR DEBUGGING
                    #sys.stderr.write("\n\\\\\\\***GAP NOTICED after "
                    #    "position {}****////".format(end)) # FOR DEBUGGING
                    # set gap values to zero
                    for x in range(end,int(info[1])):
                        # there cannot already be an entry for a position. If 
                        # so, raise error.
                        assert x+1 not in score_per_position_dict, ("{} "
                            "already has a value assigned in "
                            "`score_per_position_dict`.\nThere cannot be more "
                            "than one score per base.".format(x)) # the `+ 1` 
                        # is so the bedgraph coordinates are switchd to 
                        # normal
                        score_per_position_dict[x+1] = 0
                '''

                start = int(info[1])
                end = int(info[2])
                score = int(info[3])
                for x in range(start,end):
                    # there cannot already be an entry for a position. If so, 
                    # raise error.
                    assert x+1 not in score_per_position_dict, ("{} already has "
                        "a value assigned in `score_per_position_dict`.\nThere "
                        "cannot be more than one score per base.".format(x)) # 
                    # the `+ 1` is so the bedgraph coordinates are switchd to 
                    # normal
                    score_per_position_dict[x+1] = score
        graph_end = end
        #print ("start extracted from bedgraph file=",graph_start) # FOR DEBUGGING
        #print ("end extracted from bedgraph file=",graph_end) # FOR DEBUGGING
    return score_per_position_dict
                    




def generate_output_file_name(prefix, suffix, sample_id=None):
    '''
    Takes a range as an argument and returns string for the name of the
    output file. The generated name is based on provided info.
    If sample_id is added that will get included in file name too.


    Specific examples
    =================
    Calling function with
        ("200-500", "_ribonucleotide_plot.png")
    returns
        "200-500_ribonucleotide_plot.png

    Calling function with
        ("200-500", "_ribonucleotide_plot.png", sample_id = "WT1a")
    returns
        "200-500_WT1a_ribonucleotide_plot.png
    '''
    if sample_id:
        return prefix + "_" +sample_id + suffix
    return prefix + suffix



###--------------------------END OF HELPER FUNCTIONS---------------------------###
###--------------------------END OF HELPER FUNCTIONS---------------------------###









#*******************************************************************************
###------------------------'main' function of script---------------------------##

def plot_5prime_end_from_bedgraph(
    range_info, forward_file, reverse_file, reverse_orientation=False,
    title_prefix=title_prefix, custom_title_prefix=None, sample_id = None, 
    svg = False, limit=None, single_peak_annotation = False, 
    no_peak_annotation = False, no_dist_annotation = False, large_size = False, 
    export_data = False):
    '''
    Main function of script. 
    Plot 5'-ends for forward & reverse strand data provided in Bedgraph format.
    
    '''


    #Read forward bedgraph file making data for position betweent the first and 
    # last coordinate provided in file. Use zero for any positions with no
    # score in between those two.
    # As http://genome.ucsc.edu/goldenPath/help/bedgraph.html says
    # chromosome positions are zero-based.
    score_per_position_forward_dict = get_scores_from_bedgraph(forward_file)

    #Read forward bedgraph file making data for position betweent the first and 
    # last coordinate provided in file. Use zero for any positions with no
    # score in between those two.
    # As http://genome.ucsc.edu/goldenPath/help/bedgraph.html says
    # chromosome positions are zero-based.
    score_per_position_reverse_dict = get_scores_from_bedgraph(reverse_file)


    # convert the range string provided to start and end variables
    info= range_info.split("-")
    start,end = int(info[0]),int(info[1])
    # make a list of the plot coordinates
    plot_positions = list(range(start,end+1))



    # convert dictionary of starts to lists. This is where gaps get filled to
    # zero
    starts_per_position_forward = [score_per_position_forward_dict[x] 
                                if x in score_per_position_forward_dict else 
                                0 for x in plot_positions]
    starts_per_position_reverse = [score_per_position_reverse_dict[x] 
                                if x in score_per_position_reverse_dict else 
                                0 for x in plot_positions]





    # reverse data if plotting reverse positions orientation
    #FOR DEBUGGING
    logging.debug(plot_positions)
    logging.debug(starts_per_position_forward)
    logging.debug(starts_per_position_reverse)
    if reverse_orientation:
        plot_positions.reverse()
        starts_per_position_forward.reverse()
        starts_per_position_reverse.reverse()
    #FOR DEBUGGING
    logging.debug(plot_positions)
    logging.debug(starts_per_position_forward)
    logging.debug(starts_per_position_reverse)






    ### Determine most frequent 5' ends position for forward and reverse strand##
    # (this will be outer most for each)
    max_starts_forward = max(starts_per_position_forward)
    max_starts_reverse = max(starts_per_position_reverse)

    #now get positions that correspond to each of these (this will be outer most for each)
    for indx, value in enumerate(starts_per_position_forward):
        if value == max_starts_forward:
            most_frequent_end_pos_forward = plot_positions[indx]
            break #so only gets first one 

    for indx, value in enumerate(reversed(starts_per_position_reverse)):
        # based on Eugene's answer at 
        # https://stackoverflow.com/questions/34438901/finding-the-last-occurrence-of-an-item-in-a-list-python
        if value == max_starts_reverse:
            most_frequent_end_pos_reverse = plot_positions[len(plot_positions)-(indx +1)]
            break #so only gets first one from far side

    frequent_pos_list = ([most_frequent_end_pos_forward, 
        most_frequent_end_pos_reverse ]) # for later use
    max_ends_list = [max_starts_forward, max_starts_reverse] # for later use

    distance_between_the_most_frequent_pos = (
        max(frequent_pos_list) - min(frequent_pos_list) + 1)
    #FOR DEBUGGING
    logging.debug(max_starts_forward)
    logging.debug(most_frequent_end_pos_forward)
    logging.debug(max_starts_reverse)
    logging.debug(most_frequent_end_pos_reverse)
    logging.debug(distance_between_the_most_frequent_pos)
    # REPORT FINDINGS
    sys.stderr.write("\nPeak number of 5'-ends in forward is {}, which occurs "
        "at position "
        "{}.".format(str(max_starts_forward),str(
        most_frequent_end_pos_forward)))
    sys.stderr.write("\nPeak number of 5'-ends in reverse is {}, which occurs "
        "at position "
        "{}.".format(str(max_starts_reverse),str(
        most_frequent_end_pos_reverse)))
    sys.stderr.write(
        "\nThe distance between the most abundant sites of 5'-ends is {}.\
        ".format(str(distance_between_the_most_frequent_pos)))










    ## ***********EXPORT EXTRACTED INFORMATION OR PLOT IT*****************
    if export_data:
        ## send data output to stdout so easily redirected separate from stderr that was used to make notes.
        sys.stderr.write("\nThe data that would have been used to make the "
            "plot follows (unless ` > [FILE_NAME]` was\nused to redirect the "
                "data to a file; in that case, examine the contents of that "
                "file):\n")
        print ("position,starts_per_position_forward,starts_per_position_reverse") # for column headings for CSV output
        for indx,pos in enumerate(plot_positions):
            print (str(pos) + "," + str(starts_per_position_forward[indx]) + "," + str(starts_per_position_reverse[indx]))




    else:
        ### MAKE PLOT ####

        plt.close()
        plt.style.use(plot_style)
        f = plt.figure()
        #ax = f.add_axes((0.1, 0.09, 0.88, 0.85))
        ax = plt.axes()


        if title_prefix is not None:
            if custom_title_prefix is not None:
                title_prefix = custom_title_prefix
            else:
                title_prefix = title_prefix + ":" + str(
                    plot_positions[0]) + "-"+  str(plot_positions[-1])
            title = title_prefix + title_main
            if sample_id:
                title += " ({})".format(sample_id)
            plt.title(title, fontsize=18)
        if len(plot_positions) > 800:
            plt.vlines(plot_positions, 
                [0], 
                starts_per_position_forward, 
                colors=['C0'], 
                label = "forward 5′-DNA ends") # `0` to `[0]` based on example at 
                #http://matplotlib.org/1.2.1/examples/pylab_examples/vline_demo.html
            plt.vlines(plot_positions, 
                [0], 
                starts_per_position_reverse, 
                colors=['C1'], 
                label = "reverse 5′-DNA ends") # `0` to `[0]` based on example at 
                #http://matplotlib.org/1.2.1/examples/pylab_examples/vline_demo.html
            sys.stderr.write("\nMaking as a vlines plot because more than 800 "
                "positions to plot, and I had seen artifacts on bar plots of "
                "size much greater than several hundred.")
        else:
            plt.bar(plot_positions, 
                starts_per_position_forward, 
                color=['#BC5F5F'], 
                width=1.4, 
                label = "forward 5'-DNA ends") #width from 
                #https://stackoverflow.com/questions/20454120/how-to-remove-gaps-between-bars-in-matplotlib-bar-chart
            plt.bar(plot_positions, 
                starts_per_position_reverse, 
                color=['#7C71AF'], 
                width=1.4, 
                label = "reverse 5'-DNA ends") #width from 
                #https://stackoverflow.com/questions/20454120/how-to-remove-gaps-between-bars-in-matplotlib-bar-chart

        ax.set_xlabel('Genomic Position')
        ax.set_ylabel('# of Ends')
        plt.legend() # from https://stackoverflow.com/questions/19125722/adding-a-legend-to-pyplot-in-matplotlib-in-the-most-simple-manner-possible
        plt.axis('tight')


        if reverse_orientation:
            plt.gca().invert_xaxis() # from 
            # https://stackoverflow.com/questions/2051744/reverse-y-axis-in-pyplot. 
            # Because even if list provided in descending order, matplotlib will plot 
            # x-axis in ascending order by default, resorting the y-axis values 
            # accordingly to maintain corresponding values. (This can easily be checked 
            # by making small, simple patterns lists for positions and 
            # coverage_per_position just before plot.) But sometimes you actually want 
            # other direction, for example, with reverse strand of genomic region.

        if single_peak_annotation and no_peak_annotation:
            sys.stderr.write(
            "\n\n****Did you mean to set both `single_peak_annotation` and "
            "`no_peak_annotation` to\n`True`? `no_peak_annotation` setting "
            "OVERRIDING.\n")

        if not no_peak_annotation:
            good_y_for_1sttxt = max_starts_forward-(
            0.3*max(starts_per_position_forward+starts_per_position_reverse)) # the
            # `+(3.5*max(coverage_per_position))` is to adjust offser relative scale of
            #  the plot because y-axis less compressed for cases of small maximum y 
            # values and vice versa
            good_y_for_2ndtxt = max_starts_reverse+(
                0.4*max(starts_per_position_forward+starts_per_position_reverse))
            if reverse_orientation:
                good_x_for_1sttxt = (
                    most_frequent_end_pos_forward+(0.11*len(plot_positions)))
                good_x_for_2ndtxt = (
                    most_frequent_end_pos_reverse+(0.01*len(plot_positions)))
            else:
                good_x_for_1sttxt = (
                    most_frequent_end_pos_forward+(0.11*len(plot_positions)))
                good_x_for_2ndtxt = (
                    most_frequent_end_pos_reverse-(0.01*len(plot_positions)))
            #Note: I was finding while `shrink` and `width` work for arrows 
            # implemented as on page 
            # https://matplotlib.org/users/annotations_intro.html#annotations-tutorial 
            #, that those two properties cause errors and script failure when 
            # used `arrowstyle="simple"` or `arrowstyle="fancy"`;  see 
            # `plot_coverage_and_starts.py` for more about this.
            connectionstyle_for = "arc3,rad=0.3"
            if reverse_orientation:
                connectionstyle_for = "arc3,rad=-0.3"
            just_forward_annot = (max(
                max_starts_forward,max_starts_reverse) == max_starts_forward 
                and single_peak_annotation)
            if single_peak_annotation:
                just_forward_annot = (max(
                    max_starts_forward,max_starts_reverse) == max_starts_forward 
                    and single_peak_annotation)
                if just_forward_annot:
                    include_forward_annot = True
                    forward_annot_txt = "pos. "+ str(
                    most_frequent_end_pos_forward)+ ",\npeak 5'-ends"
                    include_reverse_annot = False
                else: 
                    include_forward_annot = False
                    include_reverse_annot = True
                    reverse_annot_txt = "pos. "+ str(
                    most_frequent_end_pos_reverse)+ ",\npeak 5'-ends"
            else:
                include_forward_annot = True
                include_reverse_annot = True
                forward_annot_txt = ("pos. "+ str(
                    most_frequent_end_pos_forward)+ ",\npeak 5'-ends\nfor "
                    "forward")
                reverse_annot_txt = ("pos. "+ str(
                    most_frequent_end_pos_reverse)+ ",\npeak 5'-ends\nfor "
                    "reverse")
            if include_forward_annot:
                plt.annotate(forward_annot_txt, 
                    xy=(most_frequent_end_pos_forward, max_starts_forward),
                    xytext=(good_x_for_1sttxt, good_y_for_1sttxt),
                    size=10,
                    color='#000000',
                    alpha = 0.65,
                    horizontalalignment='center',
                    arrowprops=dict(arrowstyle="fancy, head_width = 0.6", 
                        facecolor='#222222', 
                        edgecolor="#111111", 
                        alpha=0.65, 
                        connectionstyle=connectionstyle_for),
                    )
            if include_reverse_annot:
                plt.annotate(reverse_annot_txt, 
                    xy=(most_frequent_end_pos_reverse, max_starts_reverse),
                    xytext=(good_x_for_2ndtxt, good_y_for_2ndtxt),
                    size=10,
                    color='#000000',
                    alpha = 0.65,
                    horizontalalignment='center',
                    arrowprops=dict(arrowstyle="fancy", 
                        facecolor='#222222', 
                        edgecolor="#111111", 
                        alpha=0.65, 
                        connectionstyle="arc3,rad=-0.3"),
                        )




        if not no_dist_annotation:
            midpoint_between_freq_pos = sum(
            frequent_pos_list)/ float(len(frequent_pos_list))# = average or mean of the two
            good_y_for_text = max_starts_reverse+(0.06*max(
                starts_per_position_forward+starts_per_position_reverse)) # the
            # `+0.06*max(coverage_per_position)` is to adjust offset relative 
            # scale of the plot because y-axis less compressed for cases of 
            # small maximum y values and vice versa
            good_y_for_text = good_y_for_2ndtxt+0.33*good_y_for_2ndtxt
            plt.text(midpoint_between_freq_pos, good_y_for_text, str(
                distance_between_the_most_frequent_pos) + " bp span", 
                color = "#111111", alpha=0.65, 
                horizontalalignment='center', size=10)

            good_y_for_span_line = max_starts_reverse+0.04*max(
                starts_per_position_forward+starts_per_position_reverse)
            good_y_for_span_line = good_y_for_2ndtxt+0.3*good_y_for_2ndtxt
            plt.annotate("",
                xy=(most_frequent_end_pos_forward, good_y_for_span_line),
                xytext=(most_frequent_end_pos_reverse, good_y_for_span_line), 
                arrowprops=dict(arrowstyle="|-|",
                              #linewidth = 5,  #based on 
                              #http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.axvline
                              edgecolor = "#111111",
                              alpha=0.65,
                              connectionstyle="arc3,rad=0."), 
                                )
        # It is important that setting limit go after `plt.axis('tight')` 
        # because seems to wipe out setting if `y.lim` placed earlier; found by 
        # chance but confirmed by testing.
        if len(plot_positions) > 800:
            plt.ylim([0,None]) # important to get base of y-axis to hit 0 and 
            #not have some of negative region displayed for vlines plots; the 
            # `None` part is because needs two items and I was wanting what 
            # Matplotlib would do on its own and not the hack solutions I had 
            # come up with (`max(coverage_per_position)` & 
            # `max(coverage_per_position) + 0.10*max(coverage_per_position)`), 
            # and `None` for upper brings back Matplotlib default handling.

        if limit:
            plt.ylim(0,limit)


        ### END OF PLOT MAKING #####


        # save data and give feedback
        if custom_title_prefix:
            fn_prefix = custom_title_prefix
        else:
            fn_prefix = range_info
        output_file_name = generate_output_file_name(
            fn_prefix, suffix_for_saving_result, sample_id=sample_id)
        sys.stderr.write(
            "\n\nPlot image saved to: {}\n".format(output_file_name))
        if large_size:
            plt.savefig(output_file_name, dpi = (1600))  # IF NEED LARGE. Based 
            # on http://scipy-cookbook.readthedocs.io/items/Matplotlib_AdjustingImageSize.html
        else:
            plt.savefig(output_file_name) 
        # plt.savefig(output_file_name[:-4]+".pdf", orientation='landscape') # 
        # UNFORTUNATELY DOES NOT PRODUCE VECTOR GRAPHICS, unlike ReportLab's 
        # pdf output; USE SVG for that and the make PDF later.
        if svg:
            plt.savefig(output_file_name[:-4]+".svg", orientation='landscape') # FOR 
            # VECTOR GRAPHICS; useful if merging into Adobe Illustrator. Based 
            # on https://neuroscience.telenczuk.pl/?p=331 ; I think ReportLab 
            # also outputs SVG?
        




        #plt.show()






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
    kwargs['range_info'] = range_info
    kwargs['forward_file'] = forward_file
    kwargs['reverse_file'] = reverse_file
    kwargs['reverse_orientation'] = reverse_orientation
    kwargs['custom_title_prefix'] = custom_title_prefix
    kwargs['sample_id'] = sample_id
    kwargs['svg'] = svg
    kwargs['limit'] = limit
    kwargs['single_peak_annotation'] = single_peak_annotation
    kwargs['no_peak_annotation'] = no_peak_annotation
    kwargs['no_dist_annotation'] = no_dist_annotation
    kwargs['large_size'] = large_size
    kwargs['export_data'] = export_data
    plot_5prime_end_from_bedgraph(**kwargs)
    # using https://www.saltycrane.com/blog/2008/01/how-to-use-args-and-kwargs-in-python/#calling-a-function
    # to build keyword arguments to pass to the function above
    # (see https://stackoverflow.com/a/28986876/8508004 and
    # https://stackoverflow.com/a/1496355/8508004 
    # (maybe https://stackoverflow.com/a/7437238/8508004 might help too) for 
    # related help)



if __name__ == "__main__" and '__file__' in globals():
    """ This is executed when run from the command line """
    # Code with just `if __name__ == "__main__":` alone will be run if pasted
    # into a notebook. The addition of ` and '__file__' in globals()` is based
    # on https://stackoverflow.com/a/22923872/8508004
    # See also https://stackoverflow.com/a/22424821/8508004 for an option to 
    # provide arguments when prototyping a full script in the notebook.
    ###-----------------for parsing command line arguments-----------------------###
    import argparse
    parser = argparse.ArgumentParser(prog='plot_5prime_end_from_bedgraph.py',
    description="plot_5prime_end_from_bedgraph.py plots data from \
    bedgraph files for a particular region.        \
    **** Script by Wayne Decatur   \
    (fomightez @ github) ***")

    parser.add_argument("range", help="Region to plot indicated by start postion \
        followed by a dash, and then the end position. Example, `200-500`. \
        REQUIRED.", metavar="start-end")
    parser.add_argument("forward_data_file", help="Name of file containing the \
        forward data in bedgraph format. \
        REQUIRED.", metavar="FORWARD_DATA")
    parser.add_argument("reverse_data_file", help="Name of file containing the \
        reverse data in bedgraph format. \
        REQUIRED.", metavar="REVERSE_DATA")
    parser.add_argument("-revo", "--reverse_orientation", help="Plot left to "
        "right to correspond to 5'-> 3' direction on reverse strand.",
        action="store_true")
    parser.add_argument('-ctp', '--custom_title_prefix', action='store', 
        type=str, help="Provide a string after this flag to use that text as \
        the start of the plot title. Typically, the region info is used to \
        form the tile of the plot image to be saved; however, if this is \
        provided it will be instead used as part of the file name.")
    parser.add_argument('-si', '--sample_id', action='store', 
        type=str, help="Provide a string after this flag to use that text as \
        the sample id/strain name to be in parentheses at end of the plot \
        title. If provided, it will also be incorporated into name of saved \
        plot.")
    parser.add_argument('-svg', '--svg', help="Also save as a vector graphics \
        file (SVG) useful for working in subsequent graphics editing.",
        action="store_true") 
    parser.add_argument('-lim', '--limit', action='store', type=int, metavar="INT",
        help="**FOR ADVANCED USE.*** Allows for controlling the upper end of scale \
        for y-axis. Should only be needed when extremes of values in one plot might\
        not match a plot for a related but different sample. Put the upper limit\
        integer after the flag, such as `--limit 82`.")
    parser.add_argument("-spa", "--single_peak_annotation", help="Only show \
        annotation of the maximum peak for both forward and reverse strand \
        data.",action="store_true")
    parser.add_argument("-npa", "--no_peak_annotation", help="Turn off \
        annotation of the peak points of 5'-ends for both forward and reverse \
        data.",action="store_true")
    parser.add_argument("-nda", "--no_dist_annotation", help="Turn off \
        annotation of the distance between peak points of forward and reverse \
        data.",action="store_true")
    parser.add_argument("-ls", "--large_size", help="Print large image of \
        plot.",action="store_true")
    parser.add_argument("-exp_d", "--export_data", help="Export the data as CSV \
        printed to the console as stdout; use a redirect in a shell command to \
        send the output to a file. If the `--reverse_orientation` option is \
        utilized when the script is called then the data is exported as of the \
        complementary strand 5'-end is the starting point. This is a utility \
        feature added to enable easily passing the data mined by this script onto \
        related scripts. This overrides the plotting action, i.e, in this mode no \
        plot will be made, despite the name of the script.",action="store_true")


    #I would also like trigger help to display if no arguments provided because need at least one input file
    if len(sys.argv)==1:    #from http://stackoverflow.com/questions/4042452/display-help-message-with-python-argparse-when-script-is-called-without-any-argu
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    range_info = args.range
    forward_file = args.forward_data_file
    reverse_file = args.reverse_data_file
    reverse_orientation = args.reverse_orientation
    custom_title_prefix = args.custom_title_prefix
    sample_id = args.sample_id
    svg = args.svg
    limit = args.limit
    single_peak_annotation = args.single_peak_annotation
    no_peak_annotation = args.no_peak_annotation
    no_dist_annotation = args.no_dist_annotation
    large_size = args.large_size
    export_data = args.export_data



    main()

#*******************************************************************************
###-***********************END MAIN PORTION OF SCRIPT***********************-###
#*******************************************************************************
