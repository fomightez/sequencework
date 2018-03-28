#!/usr/bin/env python
# -*- coding: utf-8 -*-


# plot_coverage_and_starts.py by Wayne Decatur
__author__ = "Wayne Decatur" #fomightez on GitHub
__version__ = "0.0.1"
__license__ = "MIT"


#*******************************************************************************
# USES Python 2.7 but should be convertable via 2to3, see https://docs.python.org/3.0/library/2to3.html
#
# PURPOSE: Work in progress script to take output from two sources and plot
# read coverage per position and number of starts per position for a region.
# The data from two sources needed to run this script: 1. output from 
# `samtools depth` command run with `-a` option, and 2. output from `samtools
# mpileup`. 
# This data is then used to plot the coverage for a region, like in IGV's 
# read coverage track (see 
# https://software.broadinstitute.org/software/igv/AlignmentData#coverage).
# And also the number of reads with starts at each position are plotted as well
# in a manner where it is stacked onto the coverage bar to produce a stacked
# bar chart where the reads starting are seen as a component of the total 
# coverage per position.
# Need because want more flexible output of that which can be combined with 
# other data and want to bypass the MEMORY(?) limits we encounter when opening 
# large BAM files of all our SINGLE experimental data with IGV on our computer, 
# plus don't necessarily(?), always(?) want collapsing, i.e., downsampling, IGV's 
# read coverage track does. Plus, need way to plot 5' to 3' of "reverse" strand
# that IGV appears to not have. 
# 
# **"Starts" here means 5'-end of reads as IGV would show as rectangular end of
# the arrows it used to represents reads. Pertinent to this, actual 
# documentation for the pileup format seems to be at 
# http://samtools.sourceforge.net/pileup.shtml and it seems the use of the `^` 
# or `$` symbols comes from their use in `CALF` format and 
# [there](http://www.phrap.org/phredphrap/calf.pdf) specifically states " 
# 'Start' means the beginning, or left-most end, of the read's alignment, not 
# the 5' end of the read".
#
# **IMPORTANT**: `samtools depth -a` command DEFAULT output seems to match IGV's
# depth coverage track (see http://software.broadinstitute.org/software/igv/AlignmentData) 
# while the one listed in `samtools mpileup` doesn't quit match because, by 
# default some filtering/quality control happens. This is why I need data from
# both sources.
# Also note that in Samtools documentation for says for depth command it will 
# "Truncate reported depth at a maximum of INT reads. [8000]" and so there is
# a limit what it can output, and so indirectly there can be some 
# "downsampling" occuring upstream of the data the script uses.
#
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
# TO RUN:
# Example,
# Enter on the command line of your terminal, the line
#----------------------------------------------------
# python plot_coverage_and_starts.py data.coverage data.pileup
#-----------------------------------------------------
# 
# Other example comamnds:
# python plot_coverage_and_starts.py wt3.coverage wt3_pileup.pileup -rev
# python plot_coverage_and_starts.py coverage_wt1.coverage wt1_pileup.pileup -rev -exp_d

#
#
#*******************************************************************************
#


#*******************************************************************************
##################################
#  USER ADJUSTABLE VALUES        #

##################################
#
suffix_for_saving_result = "_covNstarts_plot.png"

title_prefix = "Region " # Change to `None` to suppress title

title_suffix = " coverage and read starts"

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
from scipy.stats import mode

#DEBUG CONTROL
#comment line below to turn off debug print statements
#logging.basicConfig(stream=sys.stderr, level=logging.DEBUG)

###---------------------------HELPER FUNCTIONS---------------------------------###

def count_5prime_ends(bases_string):
    '''
    Takes a bases string for a position represented in a pileup file and 
    determines the number of 5' ends occuring at that position.

    see `BIG ISSUE` and "THE RULES" below for the pseudocode describing this
    determination.

    returns the number (integer) of the 5' ends at the position.
    '''
    if "^" not in bases_string and "$" not in bases_string:
        return 0
    else:
        # since string has at least one "^" or "$", I need to check each and see
        # if it could be the 5'- end of a read. 
        nunber_of_5prime_ends = 0
        for indx,char in enumerate(bases_string):
            if char == "^":
                determining_character = bases_string[indx+2]
                if determining_character in "AGCTN.":
                    nunber_of_5prime_ends += 1
            elif char == "$":
                determining_character = bases_string[indx-1]
                if determining_character in "agctn,":
                    nunber_of_5prime_ends += 1
        return nunber_of_5prime_ends


def slice_it(li, cols=2):
    '''
    gnerator to make n chunks of a list.

    based on http://www.garyrobinson.net/2008/04/splitting-a-pyt.html
    '''
    start = 0
    for i in range(cols):
        stop = start + len(li[i::cols])
        yield li[start:stop]
        start = stop





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
    if large_size:
        suffix =  suffix .replace(".png", "LARGE.png")
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
parser = argparse.ArgumentParser(prog='plot_coverage_and_starts.py',
    description="plot_coverage_and_starts.py plots data from \
    BAM/SAM files obtained using SAMTools.        \
    **** Script by Wayne Decatur   \
    (fomightez @ github) ***")

parser.add_argument("coverage_file", help="Name of file containing the \
    tab-separated data from the command `samtools depth -a -r \
    <chrom><start>-<end> INPUT.bam > output.coverage`. REQUIRED.", type=argparse.FileType('r'), metavar="COVERAGE")
parser.add_argument("pileup_file", help="Name of file containing the \
    tab-separated data from the command `samtools mpileup -r \
    <chrom><start>-<end> INPUT.bam > output.coverage`. REQUIRED.", type=argparse.FileType('r'), metavar="PILEUP")
parser.add_argument("-rev", "--reverse_orientation", help="Plot left to right to correspond to 5'-> 3' direction on reverse strand.",action="store_true")
parser.add_argument("-npa", "--no_peak_annotation", help="Turn off annotation \
    of the peak points of 5'-ends at both ends of the plot.",action="store_true")
parser.add_argument("-nda", "--no_dist_annotation", help="Turn off annotation \
    of the distance between peak points of 5'-ends at both ends of the plot.",action="store_true")
parser.add_argument("-ls", "--large_size", help="Print large image of plot.",action="store_true")
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
coverage_file = args.coverage_file
pileup_file = args.pileup_file
reverse_orientation = args.reverse_orientation
no_peak_annotation = args.no_peak_annotation
no_dist_annotation = args.no_dist_annotation
large_size = args.large_size
limit = args.limit
export_data = args.export_data









###-----------------Actual Main portion of script---------------------------###


# make some lists to store the data
chrm = []
positions = []
coverage_per_position = []
starts_per_position_dict = {}


lines_processed = 0
new_file_text = ""

# open coverage file and start reading
sys.stderr.write("\nReading coverage file and parsing..")
#input_file_stream = open(list_file, "r") # Don't need separate open when use `type=argparse.FileType`. It sets everything up automatically and you will actually cause errors if try to open when already open.

### GUIDE TO INPUT FILE ####
# from https://www.biostars.org/p/104063/
# "file will have 3 columns (Chr#, position and depth at that position) 
# like below." BUT THERE IS NO HEADER IN ACTUAL FILE.!!!

for line in coverage_file:
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
sys.stderr.write( "\n"+ str(lines_processed) + " lines read from '" + coverage_file.name + "'.")
#input_file_stream.close() # Don't need because argparse handles when use `type=argparse.FileType`.

#sys.stderr.write("\nPositions with data span: {} - {}.".format(str(min(positions)),str(max(positions)) ))
sys.stderr.write("\nPositions within data span: {} - {}.".format(str(positions[0]),str(positions[-1]) ))
if reverse_orientation:
    sys.stderr.write("\n\\\\\\\\***THIS IS A REVERSE STRAND FEATURE****////")

# using positions list to plot and so need data for each position or will result
# in hidden gaps.  
assert len(positions) == ((positions[-1] - positions[0]) + 1), "*** Depth is not reported for all positions ***. SUGGESTED REMEDY: for `samtools depth` use the `-a` option to generate this data."
# alternatively, I could have counted from first and last and entered 0 for any 
# points without data; however, this is not a good idea for two main reasons: 
# 1. What if the missing data happens first or last (or both) positions, then 
# the covered range would be shorter than user had wanted. 2. In genral, best 
# user know about data fully and use appropriate options for commands. 


# open pileup file and start reading
sys.stderr.write("\nReading pileup file and parsing..")
#input_file_stream = open(list_file, "r") # Don't need separate open when use `type=argparse.FileType`. It sets everything up automatically and you will actually cause errors if try to open when already open.

### GUIDE TO INPUT FILE ####
# from https://en.wikipedia.org/wiki/Pileup_format:
# Each line consists of 5 (or optionally 6) tab-separated columns:
# 1. Sequence identifier
# 2. Position in sequence (starting from 1)
# 3. Reference nucleotide at that position
# 4. Number of aligned reads covering that position (depth of coverage)
# 5. Bases at that position from aligned reads
# 6. quality of those bases (OPTIONAL)
# Column 5: The bases string:
# . (dot) means a base that matched the reference on the forward strand
# , (comma) means a base that matched the reference on the reverse strand
# AGTCN denotes a base that did not match the reference on the forward strand
# agtcn denotes a base that did not match the reference on the reverse strand
# A sequence matching the regular expression \+[0-9]+[ACGTNacgtn]+ denotes an insertion of one or more bases starting from the next position
# A sequence matching the regular expression -[0-9]+[ACGTNacgtn]+ denotes a deletion of one or more bases starting from the next position
# ^ (caret) marks the start of a read segment and the ASCII of the character following `^' minus 33 gives the mapping quality 
# $ (dollar) marks the end of a read segment
# * (asterisk) is a placeholder for a deleted base in a multiple basepair deletion that was mentioned in a previous line by the -[0-9]+[ACGTNacgtn]+ notation
# Column 6: The base quality string 
# This is an optional column. If present, the ASCII value of the character minus 33 gives the mapping Phred quality of each of the bases in the previous column 5. This is similar to quality encoding in the FASTQ format.
# And there is no header in this file.
###** BIG ISSUE: `start of a read segment` DOESN'T CORRESPOND EXACTLY TO WHAT IGV SHOWS AS A 'START'. See below about sorting what a 'TRUE START' is from this data.
lines_processed = 0
for line in pileup_file:
    lines_processed += 1
    line = line.strip() # don't want line endings so I can easily work with it
    items = line.split()
    current_position= int(items[1])
    bases_string = items[4]
    assert current_position not in starts_per_position_dict, "There should only be one entry per position in the pileup."
    starts_per_position_dict[current_position] = count_5prime_ends(bases_string) # see `BIG ISSUE` and "THE RULES" below


# BIG ISSUE: Cannot just use `bases_string.count("^")` to get starts because 
# only 'TRUE START' for those on forward strand. 'TRUE START' to me is what IGV displays as rectangular end of its arrow segments in the alignment viewer, i.e. the 5' end. 
# So have to assess what strand it is on.
# This is indeed in the pileup format data:
'''
****---***
THE RULES
****---***

Ends whether 5' or 3' are marked by `^` or `$`. (You can count those for number of total read ends at a position, if you wanted that.)
**--------------------------**

`^` two characters before an uppercase AGCTN or `.`(dot) means a 5'-end (TRUE READ START).

`$` one character after a lowercase agctn or `,`(comma) means a 5'-end (TRUE READ START).

**--------------------------**
I won't need this for a script detecting 5'-ends but for sake of completeness , for 3' end...

`^` two characters before a lowercase agctn or `,`(comma) means a 3'-end (true read END).

`$` one character after an uppercase AGCTN or `.`(dot) means a 3'-end (true read END).



SOURCE OF THE "THE RULES"
Based on looking beyond the (maybe poorly written?) wikipedia documentation of 
the pileup format, the actual documentation seems to be at 
http://samtools.sourceforge.net/pileup.shtml and it seems the use of the `^` or 
`$` symbols comes from their use in `CALF` format and it says 
[there](http://www.phrap.org/phredphrap/calf.pdf)  that the `start marker` is 
before the base it refers to while the `stop marker` PRECEDES what it refers to. 
Moreover, looking at real `mpileup`-generated data, I now see that `~` always 
occurs three characters in from the end of the `Bases String` while `$` occurs 
as the second position of that string!!! So for `~` the quality is placed after 
it and then the base indicator with strand/info. For `$`, it simply means the 
one before it is the `start` or `end` depending on the direction which it also 
encode!!! Plus it says 
[in CALF documentation](http://www.phrap.org/phredphrap/calf.pdf) specifically 
states " 'Start' means the beginning, or left-most end, of the read's alignment, 
not the 5' end of the read".
The actual documentation I found is not overtly clear about this stuff in the 
`bases string` vs. the `base quality string` and I partly worked it out by 
looking into that documentation and comparing it to the `Number of aligned reads
covering that position` column suggest it is the case, plus looking at IGV 
output for the area I was looking at.
'''

# convert dictionary of starts to a list
starts_per_position = [starts_per_position_dict[x] if x in starts_per_position_dict else 0 for x in positions]

# Completed scan of input file and therefore close file, alert user as to any
# issues, and write new file.
sys.stderr.write( "\n"+ str(lines_processed) + " lines read from '" + pileup_file.name + "'.")
#input_file_stream.close() # Don't need because argparse handles when use `type=argparse.FileType`.







# reverse data if plotting reverse positions orientation
#FOR DEBUGGING
logging.debug(positions)
logging.debug(coverage_per_position)
logging.debug(starts_per_position)
if reverse_orientation:
    positions.reverse()
    coverage_per_position.reverse()
    starts_per_position.reverse()
#FOR DEBUGGING
logging.debug(positions)
logging.debug(coverage_per_position)
logging.debug(starts_per_position)






### Determine most frequent 5' ends position for 1st and 2nd half #####
# first split in half
num_of_chunks = 2 # because want split roughly in half
pos_chunks_generator = slice_it(positions,num_of_chunks)
starts_chunks_generator = slice_it(starts_per_position,num_of_chunks)
pos_chunks = []
starts_chunks = []
for x in range(0,num_of_chunks):
    pos_chunks.append(next(pos_chunks_generator))
    starts_chunks.append(next(starts_chunks_generator ))

max_ends_first_half = max(starts_chunks[0])
max_ends_second_half = max(starts_chunks[1])

index_max_ends_first_half = starts_chunks[0].index(max_ends_first_half) #so 
# gets first occurence from start side
most_frequent_end_pos_first_half = pos_chunks[0][index_max_ends_first_half]

for indx, value in enumerate(reversed(starts_chunks[1])):
    # based on Eugene's answer at 
    # https://stackoverflow.com/questions/34438901/finding-the-last-occurrence-of-an-item-in-a-list-python
    if value == max_ends_second_half:
        index_max_ends_second_half = len(starts_chunks[1]) - indx - 1
        break #so only gets first one from far side

most_frequent_end_pos_second_half = pos_chunks[1][index_max_ends_second_half]

frequent_pos_list = ([most_frequent_end_pos_first_half, 
    most_frequent_end_pos_second_half ]) # for later use
max_ends_list = [max_ends_first_half, max_ends_second_half] # for later use

distance_between_the_most_frequent_pos = (
    max(frequent_pos_list) - min(frequent_pos_list) + 1)
#FOR DEBUGGING
logging.debug(max_ends_first_half)
logging.debug(most_frequent_end_pos_first_half)
logging.debug(max_ends_second_half)
logging.debug(most_frequent_end_pos_second_half)
logging.debug(distance_between_the_most_frequent_pos)
# REPORT FINDINGS
sys.stderr.write("\nPeak number of 5'-ends in first half is {}, which occurs \
on left side of plot at position {}.".format(str(max_ends_first_half),str(
    most_frequent_end_pos_first_half)))
sys.stderr.write("\nPeak number of 5'-ends in second half is {}, which occurs \
on right side of plot at position {}.".format(str(max_ends_second_half),str(
    most_frequent_end_pos_second_half)))
sys.stderr.write(
    "\nThe distance between the most abundant sites of 5'-ends is {}.\
    ".format(str(distance_between_the_most_frequent_pos)))










## ***********EXPORT EXTRACTED INFORMATION OR PLOT IT*****************
if export_data:
    ## send data output to stdout so easily redirected separate from stderr that was used to make notes.
    sys.stderr.write("\nThe data that would have been used to make the plot follows (unless ` > [FILE_NAME]` was\nused to redirect the data to a file; in that case, examine the contents of that file):\n")
    print ("position,coverage_per_position,starts_per_position") # for column headings for CSV output
    for indx,pos in enumerate(positions):
        print (str(pos) + "," + str(coverage_per_position[indx]) + "," + str(starts_per_position[indx]))




else:
    ### MAKE PLOT ####

    plt.close()
    plt.style.use(plot_style)
    f = plt.figure()
    #ax = f.add_axes((0.1, 0.09, 0.88, 0.85))
    ax = plt.axes()


    if title_prefix is not None:
        #title = title_prefix + chrm[0] + ":" + str(max(positions)) + "-"+  str(
        #   min(positions)) + title_suffix
        title = title_prefix + chrm[0] + ":" + str(positions[0]) + "-"+  str(
            positions[-1]) + title_suffix
        
        plt.title(title, fontsize=18)
    if len(positions) > 800:
        plt.vlines(positions, [0], coverage_per_position, colors=['C0'], label = "coverage") # `0` to `[0]` based on example at http://matplotlib.org/1.2.1/examples/pylab_examples/vline_demo.html
        plt.vlines(positions, [0], starts_per_position, colors=['C1'], label = "starts (5'ends)") # `0` to `[0]` based on example at http://matplotlib.org/1.2.1/examples/pylab_examples/vline_demo.html
        sys.stderr.write("\nMaking as a vlines plot because more than 800 positions to plot, and I had seen artifacts on bar plots of size much greater than several hundred.")
    else:
        plt.bar(positions, coverage_per_position, width=1.0, label = "coverage") #width from https://stackoverflow.com/questions/20454120/how-to-remove-gaps-between-bars-in-matplotlib-bar-chart
        plt.bar(positions, starts_per_position, width=1.0, label = "starts (5'ends)") #width from https://stackoverflow.com/questions/20454120/how-to-remove-gaps-between-bars-in-matplotlib-bar-chart

    plt.axes().set_xlabel('Genomic Position')
    plt.axes().set_ylabel('# of Reads')
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

    if not no_peak_annotation:
        good_y_for_1sttxt = max_ends_first_half+(0.4*max(coverage_per_position)) # the
        # `+(3.5*max(coverage_per_position))` is to adjust offser relative scale of
        #  the plot because y-axis less compressed for cases of small maximum y 
        # values and vice versa
        good_y_for_2ndtxt = max_ends_second_half+(0.4*max(coverage_per_position))
        if reverse_orientation:
            good_x_for_1sttxt = most_frequent_end_pos_first_half-(0.01*len(positions))
            good_x_for_2ndtxt = most_frequent_end_pos_second_half+(0.01*len(positions))
        else:
            good_x_for_1sttxt = most_frequent_end_pos_first_half+(0.01*len(positions))
            good_x_for_2ndtxt = most_frequent_end_pos_second_half-(0.01*len(positions))
        #Note: I was finding while `shrink` and `width` work for arrows implemented 
        # as on page https://matplotlib.org/users/annotations_intro.html#annotations-tutorial 
        #, that those two properties cause errors and script failure when used when 
        # `arrowstyle="simple"` or `arrowstyle="fancy"` is added. KEPT ORIGINALS 
        # BELOW FOR DEMONSTRATION OF DIFFERENCES.
        plt.annotate("pos. "+ str(most_frequent_end_pos_first_half)+ \
            ",\npeak 5'-ends\nfor 1st half", \
            xy=(most_frequent_end_pos_first_half, max_ends_first_half),
            xytext=(good_x_for_1sttxt, good_y_for_1sttxt),
            size=10,
            color='#000000',
            alpha = 0.65,
            horizontalalignment='center',
            arrowprops=dict(arrowstyle="fancy, head_width = 0.6", facecolor='#222222', edgecolor="#111111", alpha=0.65, connectionstyle="arc3,rad=0.3"),
                )
        plt.annotate("pos. "+ str(most_frequent_end_pos_second_half)+ \
            ",\npeak 5'-ends\nfor 2nd half", \
            xy=(most_frequent_end_pos_second_half, max_ends_second_half),
            xytext=(good_x_for_2ndtxt, good_y_for_2ndtxt),
            size=10,
            color='#000000',
            alpha = 0.65,
            horizontalalignment='center',
            arrowprops=dict(arrowstyle="fancy", facecolor='#222222', edgecolor="#111111", alpha=0.65, connectionstyle="arc3,rad=-0.3"),
                )

        '''
        ***CODE FOR PEAK ANNOTATIONS FROM ORIGINAL WHERE `SHRINGK` AND `WIDTH` WORKED FOR ARROWS:**
        plt.annotate("position "+ str(most_frequent_end_pos_first_half)+ \
            ",\nmost frequent read start\nsite for 1st half", \
            xy=(most_frequent_end_pos_first_half, max_ends_first_half),
            xytext=(most_frequent_end_pos_first_half-0.3, max_ends_first_half+13.7),
            size=10,
            horizontalalignment='center',
            # arrowprops=dict(facecolor='black', shrink=0.05),
            arrowprops=dict(facecolor='#444444', shrink=0.027, width = 4, edgecolor="#333333"),

                )
        plt.annotate("position "+ str(most_frequent_end_pos_second_half)+ \
            ",\nmost frequent read start\nsite for 2nd half", \
            xy=(most_frequent_end_pos_second_half, max_ends_second_half),
            xytext=(most_frequent_end_pos_second_half+0.3, max_ends_second_half+13.7),
            size=10,
            horizontalalignment='center',
            # arrowprops=dict(facecolor='black', shrink=0.05),
            arrowprops=dict(facecolor='#444444', shrink=0.027, width = 3, edgecolor="#333333"),
                )
        '''



    if not no_dist_annotation:
        midpoint_between_freq_pos = sum(
        frequent_pos_list)/ float(len(frequent_pos_list))# = average or mean of the two
        good_y_for_text = max(max_ends_list)+(0.06*max(coverage_per_position)) # the
        # `+0.06*max(coverage_per_position)` is to adjust offser relative scale of
        #  the plot because y-axis less compressed for cases of small maximum y 
        # values and vice versa
        plt.text(midpoint_between_freq_pos, good_y_for_text, str(
            distance_between_the_most_frequent_pos) + " bp span", color = "#111111", alpha=0.65, horizontalalignment='center', size=10)

        good_y_for_span_line = max(max_ends_list)+0.04*max(coverage_per_position)
        plt.annotate("",
            xy=(most_frequent_end_pos_first_half, good_y_for_span_line),
            xytext=(most_frequent_end_pos_second_half, good_y_for_span_line), 
            arrowprops=dict(arrowstyle="|-|",
                          #linewidth = 5,  #based on http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.axvline
                          edgecolor = "#111111",
                          alpha=0.65,
                          connectionstyle="arc3,rad=0."), 
                            )
    # It is important that setting limit go after `plt.axis('tight')` because seems to wipe out setting if `y.lim` placed earlier; found by chance but confirmed by testing.
    if len(positions) > 800:
        plt.ylim([0,None]) # important to get base of y-axis to hit 0 and not have some of negative region displayed for vlines plots; the `None` part is because needs two items and I was wanting what Matplotlib would do on its own and not the hack soltions I had some up with (`max(coverage_per_position)` & `max(coverage_per_position) + 0.10*max(coverage_per_position)`), and `None` for upper brings back Matplotlib default handling.

    if limit:
        plt.ylim(0,limit)


    ### END OF PLOT MAKING #####


    # save data and give feedback
    output_file_name = generate_output_file_name(coverage_file.name, suffix_for_saving_result)
    sys.stderr.write("\n\nPlot image saved to: {}\n".format(output_file_name))
    if large_size:
        plt.savefig(output_file_name, dpi = (1600))  # IF NEED LARGE. Based on http://scipy-cookbook.readthedocs.io/items/Matplotlib_AdjustingImageSize.html
    else:
        plt.savefig(output_file_name) 
    # plt.savefig(output_file_name[:-4]+".pdf", orientation='landscape') # UNFORTUNATELY DOES NOT PRODUCE VECTOR GRAPHICS, unlike ReportLab's pdf output; USE SVG for that and the make PDF later.
    # plt.savefig(output_file_name[:-4]+".svg", orientation='landscape') # FOR VECTOR GRAPHICS; useful if merging into Adobe Illustrator. Based on https://neuroscience.telenczuk.pl/?p=331 ; I think ReportLab also outputs SVG?
    




    #plt.show()








'''
KEEP BECAUSE MAY WANT TO UNCOMMENT (OR COPY THE BLOCK) AND EDIT CODE SO CAN DO VLINES FOR SMALL REGIONS AS WELL, PLUS IS A RECORD OF WHAT TRIED AND WORKED OUT FOR VLINES PLOT HANDLING:
# I was seeing in the plot_coverage_and_starts.py some of the bars for read starts getting obscured because worked out width of bar char with size of origin alignments, which are slightly less than 500 bp, but now I am trying some much larger datasets and seeing artifacts again where things I know are there are obscured.
# Testing a fix of using the type of plot I used on the `plot_expression.py` when much larger than an origin. Starting with this plot because simpler situation to test code.
if len(positions) > 800:
    plt.close()

    #colors = (['C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'xkcd:magenta'])
    #from itertools import cycle
    #colors = cycle(colors)
    #color = next(colors)

    #cs.extend([color] * len(positons))
    plt.vlines(positions, [0], coverage_per_position, colors=['C0'], label = "coverage") # `0` to `[0]` based on example at http://matplotlib.org/1.2.1/examples/pylab_examples/vline_demo.html
    plt.vlines(positions, [0], starts_per_position, colors=['C1'], label = "starts (5'ends)") # `0` to `[0]` based on example at http://matplotlib.org/1.2.1/examples/pylab_examples/vline_demo.html

    #plt.bar(positions, coverage_per_position) #<== THIS RESULS IN MASSIVE AMOUNTS OF ARTIFACTS THAT MAKE IT LOOK SPIKY. WAY WORSE THAN ANOYANCE WITH vlines bottom not being just at 0.

    #ax = sns.barplot(positions, coverage_per_position, color = "C0") # COULDN'T USE BECAUSE IT WANTED TO PLOT EVERY POSITION NUMBER ON X-AXIS AND TOOK MANY,MANY, TOO MANY MINUTES FOR SEVERAL THOUSAND
    sys.stderr.write("\nAlso making a vlines plot in addition to a bar plot because more than 800 positions to plot, and I had seen artifacts on bar plots of size much greater than several hundred; comparing each is suggested to make sure bar plot renders all data clearly.")
    if title_prefix is not None:
        plt.title(title, fontsize=18)
    plt.axes().set_xlabel('Genomic Position')
    plt.axes().set_ylabel('# of Reads')
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
    if limit:
        plt.ylim(0,limit)
    else:
        plt.ylim([0,None]) # important to get base of y-axis to hit 0 and not have some of negative region displayed for vlines plots; the `None` part is because needs two items and I was wanting what Matplotlib would do on its own and not the hack soltions I had some up with (`max(coverage_per_position)` & `max(coverage_per_position) + 0.10*max(coverage_per_position)`), and `None` for upper brings back Matplotlib default handling.
    output_file_name = generate_output_file_name(coverage_file.name, "_VLINES" + suffix_for_saving_result)
    #fig = ax.get_figure() # need for using `sns.barplot`
    #fig.savefig(output_file_name) 
    plt.savefig(output_file_name) 
    sys.stderr.write("\n\nPlot image with vlines saved to: {}\n".format(output_file_name))
'''




'''
# UNCOMMENT THIS BLOCK IF WANT A VLINES PLOT FOR THE SMALL REGIONS THAT NORMALLY JUST GIVE BARPLOTS. THE IDEA BEING JUST TO COMPARE AND SEE WHICH ONE MAY LOOK BEST FOR THOSE CLOSER TO CUTOFF SIZE OF "SMALL" REGIONS. THE ONLY ISSUE IS BY DEFAULT IN VLINES PLOT THE LEGEND ISN'T ADJUSTED TO NOT SHOW UP ON TOP OF DATA, LIKE IT IS IN A BARPLOT.
if len(positions) < 801:
    plt.close()

    #colors = (['C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'xkcd:magenta'])
    #from itertools import cycle
    #colors = cycle(colors)
    #color = next(colors)

    #cs.extend([color] * len(positons))
    plt.vlines(positions, [0], coverage_per_position, colors=['C0'], label = "coverage") # `0` to `[0]` based on example at http://matplotlib.org/1.2.1/examples/pylab_examples/vline_demo.html
    plt.vlines(positions, [0], starts_per_position, colors=['C1'], label = "starts (5'ends)") # `0` to `[0]` based on example at http://matplotlib.org/1.2.1/examples/pylab_examples/vline_demo.html

    #plt.bar(positions, coverage_per_position) #<== THIS RESULS IN MASSIVE AMOUNTS OF ARTIFACTS THAT MAKE IT LOOK SPIKY. WAY WORSE THAN ANOYANCE WITH vlines bottom not being just at 0.

    #ax = sns.barplot(positions, coverage_per_position, color = "C0") # COULDN'T USE BECAUSE IT WANTED TO PLOT EVERY POSITION NUMBER ON X-AXIS AND TOOK MANY,MANY, TOO MANY MINUTES FOR SEVERAL THOUSAND
    sys.stderr.write("\nAlso making a vlines plot in addition to a bar plot because I had seen artifacts on bar plot; comparing each is suggested to make sure bar plot renders all data clearly.")
    if title_prefix is not None:
        plt.title(title, fontsize=18)
    plt.axes().set_xlabel('Genomic Position')
    plt.axes().set_ylabel('# of Reads')
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
    if limit:
        plt.ylim(0,limit)
    else:
        plt.ylim([0,None]) # important to get base of y-axis to hit 0 and not have some of negative region displayed for vlines plots
    output_file_name = generate_output_file_name(coverage_file.name, "_VLINES" + suffix_for_saving_result)
    #fig = ax.get_figure() # need for using `sns.barplot`
    #fig.savefig(output_file_name) 
    plt.savefig(output_file_name) 
    sys.stderr.write("\n\nPlot image with vlines saved to: {}\n".format(output_file_name))
'''








#*******************************************************************************
###-***********************END MAIN PORTION OF SCRIPT***********************-###
#*******************************************************************************
