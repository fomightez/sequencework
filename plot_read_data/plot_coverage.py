#!/usr/bin/env python
# -*- coding: utf-8 -*-


# plot_coverage.py by Wayne Decatur
__author__ = "Wayne Decatur" #fomightez on GitHub
__version__ = "0.0.1"
__license__ = "MIT"


#*******************************************************************************
# USES Python 2.7 but should be convertable via 2to3, see https://docs.python.org/3.0/library/2to3.html
#
# PURPOSE: Work in progress script to take output from `samtools depth` command 
# run with `-a` option, and then plot the coverage for a region, like in IGV's 
# read coverage track (see 
# https://software.broadinstitute.org/software/igv/AlignmentData#coverage). 
# Need because want more flexible output of that which can be combined with 
# other data and want to bypass the MEMORY(?) limits we encounter when opening 
# large BAM files of all our SINGLE experimental data with IGV on our computer, 
# plus don't necessarily(?), always(?) want collapsing, i.e., downsampling, IGV's 
# read coverage track does. Plus, need way to plot 5' to 3' of "reverse" strand
# that IGV appears to not have. 
#
# **IMPORTANT**: `samtools depth -a` command DEFAULT output seems to match IGV's
# depth coverage track (see http://software.broadinstitute.org/software/igv/AlignmentData) 
# while the one listed in `samtools mpileup` doesn't quit match because, by 
# default some filtering/quality control happens.
# Also note that in Samtools documentation for says for depth command it will 
# "Truncate reported depth at a maximum of INT reads. [8000]" and so there is
# a limit what it can output, and so indirectly there can be some 
# "downsampling" occuring upstream of the data the script uses.
#
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
# TO RUN:
# Example,
# Enter on the command line of your terminal, the line
#-----------------------------------
# python plot_coverage.py data.coverage
#-----------------------------------
#
# Other example command:
# python plot_coverage.py coverage_wt1.coverage -rev -exp_d
#
#*******************************************************************************
#


#*******************************************************************************
##################################
#  USER ADJUSTABLE VALUES        #

##################################
#
suffix_for_saving_result = "_cov_plot.png"

title_prefix = "Region " # Change to `None` to suppress title

title_suffix = " coverage"

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

#DEBUG CONTROL
#comment line below to turn off debug print statements
#logging.basicConfig(stream=sys.stderr, level=logging.DEBUG)

###---------------------------HELPER FUNCTIONS---------------------------------###


def generate_output_file_name(file_name, suffix):
    '''
    Takes a file name as an argument and returns string for the name of the
    output file. The generated name is based on the original file
    name.

    It also notes in name specific chromsomes or scaffolds if plotting was 
    limited to those.

    Specific example
    ================
    Calling function with
        ("data1.txt", "_plot.png")
    returns
        "data1_plot.png"
    '''
    main_part_of_name, file_extension = os.path.splitext(
        file_name) #from http://stackoverflow.com/questions/541390/extracting-extension-from-filename-in-python
    return main_part_of_name + suffix




def list2text(a_list):
    '''
    a function that takes a lists and makes a string where each item in list
    is on a new line
    '''
    return "\n".join(a_list)

###--------------------------END OF HELPER FUNCTIONS---------------------------###
###--------------------------END OF HELPER FUNCTIONS---------------------------###














#*******************************************************************************
###-----------------for parsing command line arguments-----------------------###
parser = argparse.ArgumentParser(prog='plot_coverage.py',
    description="plot_coverage.py is a script for plotting data from \
    BAM/SAM files obtained using SAMTools.        \
    **** Script by Wayne Decatur   \
    (fomightez @ github) ***")

parser.add_argument("tsv_file", help="Name of file containing the tab-separated data from the command `samtools depth -a -r <chrom><start>-<end> INPUT.bam > output.coverage`. REQUIRED.", type=argparse.FileType('r'), metavar="COVERAGE")
parser.add_argument("-rev", "--reverse_orientation", help="Plot left to right to correspond to 5'-> 3' direction on reverse strand.",action="store_true")
parser.add_argument('-lim', '--limit', action='store', type=int, metavar="INT",
    help="**FOR ADVANCED USE.*** Allows for controlling the upper end of scale \
    for y-axis. Should only be needed when extremes of values in one plot might\
    not match a plot for a related but different sample. Put the upper limit\
    integer after the flag, such as `--limit 82`.") 
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
tsv_file = args.tsv_file
reverse_orientation = args.reverse_orientation
limit = args.limit
export_data = args.export_data











###-----------------Actual Main portion of script---------------------------###


# make some lists to store the data
chrm = []
positions = []
coverage_per_position = []


lines_processed = 0
new_file_text = ""

# open input file and start reading
sys.stderr.write("\nReading input file and parsing..")
#input_file_stream = open(list_file, "r") # Don't need separate open when use `type=argparse.FileType`. It sets everything up automatically and you will actually cause errors if try to open when already open.

### GUIDE TO INPUT FILE ####
# from https://www.biostars.org/p/104063/
# "file will have 3 columns (Chr#, position and depth at that position) 
# like below." BUT THERE IS NO HEADER IN ACTUAL FILE.!!!

for line in tsv_file:
    lines_processed += 1
    line = line.strip() # don't want line endings so I can easily work with it
    items = line.split()
    chrm.append(items[0])
    positions.append(items[1])
    coverage_per_position.append(items[2])

# convert lists of strings to integers
positions = [int(x) for x in positions]
coverage_per_position = [int(x) for x in coverage_per_position]





# Completed scan of input file and therefore close file, alert user as to any
# issues, and write new file.
sys.stderr.write( "\n"+ str(lines_processed) + " lines read from '" + tsv_file.name + "'.")
#input_file_stream.close() # Don't need because argparse handles when use `type=argparse.FileType`.

sys.stderr.write("\nPositions within data span: {} - {}.".format(str(positions[0]),str(positions[-1]) ))

# using positions list to plot and so need data for each position or will result
# in hidden gaps.  
assert len(positions) == ((positions[-1] - positions[0]) + 1), "*** Depth is not reported for all positions ***. SUGGESTED REMEDY: for `samtools depth` use the `-a` option to generate this data."
# alternatively, I could have counted from first and last and entered 0 for any 
# points without data; however, this is not a good idea for two main reasons: 
# 1. What if the missing data happens first or last (or both) positions, then 
# the covered range would be shorter than user had wanted. 2. In genral, best 
# user know about data fully and use appropriate options for commands. 




# reverse data if plotting reverse positions orientation
#FOR DEBUGGING
logging.debug(positions)
logging.debug(coverage_per_position)
if reverse_orientation:
    positions.reverse()
    coverage_per_position.reverse()
#FOR DEBUGGING
logging.debug(positions)
logging.debug(coverage_per_position)





## ***********EXPORT EXTRACTED INFORMATION OR PLOT IT*****************
if export_data:
    ## send data output to stdout so easily redirected separate from stderr that was used to make notes.
    sys.stderr.write("\nThe data that would have been used to make the plot follows (unless ` > [FILE_NAME]` was\nused to redirect the data to a file; in that case, examine the contents of that file):\n")
    print ("position,coverage_per_position") # for column headings for CSV output
    for indx,pos in enumerate(positions):
        print (str(pos) + "," + str(coverage_per_position[indx]))




else:
    ### MAKE PLOT ####

    plt.close()
    plt.style.use(plot_style)
    f = plt.figure()
    #ax = f.add_axes((0.1, 0.09, 0.88, 0.85))
    ax = plt.axes()


    if title_prefix is not None:
        title = title_prefix + chrm[0] + ":" + str(positions[0]) + "-"+  str(
            positions[-1]) + title_suffix
        plt.title(title, fontsize=18)

    # I was seeing in the plot_coverage_and_starts.py some of the bars for read starts getting obscured because worked out width of bar char with size of origin alignments, which are slightly less than 500 bp, but now I am trying some much larger datasets and seeing artifacts again where things I know are there are obscured.
    # Testing a fix of using the type of plot I used on the `plot_expression.py` when much larger than an origin. Starting with this plot because simpler situation to test code.
    if len(positions) > 800:
        plt.vlines(positions, [0], coverage_per_position, colors=['C0']) # `0` to `[0]` based on example at http://matplotlib.org/1.2.1/examples/pylab_examples/vline_demo.html
        sys.stderr.write("\nMaking as a vlines plot because more than 800 positions to plot, and I had seen artifacts on bar plots of size much greater than several hundred.")
    else:
        plt.bar(positions, coverage_per_position, width=1.0) #width from https://stackoverflow.com/questions/20454120/how-to-remove-gaps-between-bars-in-matplotlib-bar-chart

    plt.axes().set_xlabel('Genomic Position')
    plt.axes().set_ylabel('# of Reads')
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

    # It is important that setting limit go after `plt.axis('tight')` because seems to wipe out setting if `y.lim` placed earlier; found by chance but confirmed by testing.
    if len(positions) > 800:
        plt.ylim([0,None]) # important to get base of y-axis to hit 0 and not have some of negative region displayed for vlines plots; the `None` part is because needs two items and I was wanting what Matplotlib would do on its own and not the hack soltions I had some up with (`max(coverage_per_position)` & `max(coverage_per_position) + 0.10*max(coverage_per_position)`), and `None` for upper brings back Matplotlib default handling.

    if limit:
        plt.ylim(0,limit)


    ### END OF PLOT MAKING #####


    # save data and give feedback
    output_file_name = generate_output_file_name(tsv_file.name, suffix_for_saving_result)
    sys.stderr.write("\n\nPlot image saved to: {}\n".format(output_file_name))
    plt.savefig(output_file_name)
    # plt.savefig(output_file_name[:-4]+"LARGE.png", dpi = (1600))  # IF NEED LARGE. Based on http://scipy-cookbook.readthedocs.io/items/Matplotlib_AdjustingImageSize.html
    # plt.savefig(output_file_name[:-4]+".pdf", orientation='landscape') # UNFORTUNATELY DOES NOT PRODUCE VECTOR GRAPHICS, unlike ReportLab's pdf output; USE SVG for that and the make PDF later.
    # plt.savefig(output_file_name[:-4]+".svg", orientation='landscape') # FOR VECTOR GRAPHICS; useful if merging into Adobe Illustrator. Based on https://neuroscience.telenczuk.pl/?p=331 ; I think ReportLab also outputs SVG?
    
    

    #plt.show()





'''
KEEP BECAUSE MAY WANT TO UNCOMMENT (OR COPY BLOCK) AND EDIT CODE SO CAN DO VLINES FOR SMALL REGIONS AS WELL, PLUS IS A RECORD OF WHAT TRIED AND WORKED OUT FOR VLINES PLOT HANDLING:
# I was seeing in the plot_coverage_and_starts.py some of the bars for read starts getting obscured because worked out width of bar char with size of origin alignments, which are slightly less than 500 bp, but now I am trying some much larger datasets and seeing artifacts again where things I know are there are obscured.
# Testing a fix of using the type of plot I used on the `plot_expression.py` when much larger than an origin. Starting with this plot because simpler situation to test code.
if len(positions) > 800:
    plt.close()

    #colors = (['C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'xkcd:magenta'])
    #from itertools import cycle
    #colors = cycle(colors)
    #color = next(colors)

    #cs.extend([color] * len(positons))
    plt.vlines(positions, [0], coverage_per_position, colors=['C0']) # `0` to `[0]` based on example at http://matplotlib.org/1.2.1/examples/pylab_examples/vline_demo.html
    #plt.bar(positions, coverage_per_position) #<== THIS RESULS IN MASSIVE AMOUNTS OF ARTIFACTS THAT MAKE IT LOOK SPIKY. WAY WORSE THAN ANOYANCE WITH vlines bottom not being just at 0.

    #ax = sns.barplot(positions, coverage_per_position, color = "C0") # COULDN'T USE BECAUSE IT WANTED TO PLOT EVERY POSITION NUMBER ON X-AXIS AND TOOK MANY,MANY, TOO MANY MINUTES FOR SEVERAL THOUSAND
    sys.stderr.write("\nAlso making a vlines plot in addition to a bar plot because more than 800 positions to plot, and I had seen artifacts on bar plots of size much greater than several hundred; comparing each is suggested to make sure bar plot renders all data clearly.")
    if title_prefix is not None:
        plt.title(title, fontsize=18)
    plt.axes().set_xlabel('Genomic Position')
    plt.axes().set_ylabel('# of Reads')
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
    if limit:
        plt.ylim(0,limit)
    else:
        plt.ylim([0,None]) # important to get base of y-axis to hit 0 and not have some of negative region displayed for vlines plots; the `None` part is because needs two items and I was wanting what Matplotlib would do on its own and not the hack soltions I had some up with (`max(coverage_per_position)` & `max(coverage_per_position) + 0.10*max(coverage_per_position)`), and `None` for upper brings back Matplotlib default handling.
    output_file_name = generate_output_file_name(tsv_file.name, "_VLINES" + suffix_for_saving_result)
    #fig = ax.get_figure() # need for using `sns.barplot`
    #fig.savefig(output_file_name) 
    plt.savefig(output_file_name) 
    sys.stderr.write("\n\nPlot image with vlines saved to: {}\n".format(output_file_name))
'''





#*******************************************************************************
###-***********************END MAIN PORTION OF SCRIPT***********************-###
#*******************************************************************************
