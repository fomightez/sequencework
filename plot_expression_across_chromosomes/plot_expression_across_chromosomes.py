#!/usr/bin/env python

# plot_expression_across_chromosomes.py by Wayne Decatur
__author__ = "Wayne Decatur" #fomightez on GitHub
__license__ = "MIT"


#*******************************************************************************
# Compatible with both Python 2.7 and Python 3.6 (verified); written initially 
# in Python 2.7 to be compatible with Python 3.
#
# PURPOSE: Plot ratio of expression of experimental condition vs. wild-type (or 
# baseline state) for genes in sequential order across chromosomes in the genome. 
# Should serve to highlight regions of deviation characteristic of aneuploidy 
# or segmental duplication/deletion.
# Requires two files: 1. a file of a genome annotation format in order to parse 
# the locations of genes and (approximate) length of chromosomes; 2. a file of 
# expression data to plot. Currently, the data needs to be in a tab-delimited 
# format, but this could easily be adjusted. At this point, to accomodate 
# replicates data it is assumed you have in your data file for each gene the 
# resulting (combined) level metric for the replicates, such as the mean TPM for
# your "wild-type" samples and your mean TPM for your experimental samples. In 
# the future there will be a related script for "raw" results produced by Salmon 
# or HTSeq. The hope being you just need to point the script at the raw data 
# files and it will automagically handle the combining and produce a plot 
# showing the expression of genes across the chromosomes.
#
# There are several optional flags that can be supplied at the time of calling
# the script to control options. These are shown if you invoke with `-help` 
# flag or simply call the script with no additional arguments. Additionally, 
# inside the script there are several `USER ADJUSTABLE VALUES` that can be 
# edited for easy customization.
# Built to be general enough to be easily modified. For example, should be 
# fairly straightforward to change to plotting ratio of sequencing coverage 
# across chromosomes in a genome.
# The plotting approach and other aspects borrow from Brent Pedersen's awesome
# `manhattan-plot.py` script at the link below:
# https://github.com/brentp/bio-playground/blob/master/plots/manhattan-plot.py
# This script should produce a plot similar to a combination of Brent Pedersen's
# `manhattan-plot`and Figure 5B from Thorburn et al. 2013 (PMID: 23468524), see 
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3639041/
#
#
#
#
#
# Dependencies beyond the mostly standard libraries/modules:
# statsmodels, seaborn (although could easily be altered to make unneeded)
#
#
# VERSION HISTORY:
# v.0.1. basic working version
#
# to do:
# - make jupyter nb version (display mock data in it so don't need to have our real data out there?) which should be easy using the ipython based testing code version to get a hard coded one fast.(sort of did with little change needed) BUT NEED BETTER WAY FOR HANDLING HARDCODING OF FILENAMES STILL.
# - add supplying arguments to basic jupyter nb one where first hardcoded
# - make a static ipython notebook version with image stored. (display mock data so don't need to have our real data out there?)
# - make a jupyter nb version with working binder demo (new beta.mybinder.org)(display mock data so don't need to have our real data out there?)
# - fix so that `chr IX` shows up in sequential order in plot? Probably not important since chromsomes don't relate to arbitrary designations we have them anyway.
# - dream feature for jupyter nb version: interactive plot where gene id shows 
# up when hover (and see 
# https://moderndata.plot.ly/bioinformatics-plots-made-in-python-and-r/ for 
# other? inspiration)
# - add ability to supply chromsome data as BEd format? But then where does gene
# location get supplied???
#
#
#
# TO RUN:
# Example,
# Enter on the command line of your terminal, the line
#-----------------------------------
# python plot_expression_across_chromosomes.py genome_annotation.gtf data.tsv
#-----------------------------------
#
#
#*******************************************************************************
#


#*******************************************************************************
##################################
#  USER ADJUSTABLE VALUES        #

##################################
#

# `genome_annotation_fields` to match your genome annotation source data. You'll 
# most likely only need one of these. 
genome_annotation_fields_for_gtf = ("seqname", "source", "feature type", "start", 
    "end", "score", "strand", "frame", "attribute")
genome_annotation_fields_for_gff = ("seqname", "source", "feature type", "start", 
    "end", "score", "strand", "frame", "group")
genome_annotation_fields_for_bed = ("chrom", "chromStart", "chromEnd", "name", 
    "score", "strand", "thickStart", "thickEnd", "itemRgb", "blockCount", 
    "blockSizes", "blockStarts")

suffix_for_saving_result = "_across_chr.png"

title_prefix = "Expression across " # Change to `None` to suppress title

limit_before_rotate = 3 # upper limit of max length of chromosome or scaffold 
# names before x axis labels are all rotated when plotting all genome. Same 
# number will be used as total number of chromosomes (or scaffold) to be 
# plotting before considering length of name again.

y_cutoff = 4 # A limit was added to avoid extreme values compressing the 
# typically important range when log2 used. Adjust that limit here or 
# run with `--no_limits` flag enables. y_cutoff is not used when `no_log` flag 
# used or when values are within +/- this interval.

deviation_factor = 0.25 #measure of degree of discrepancy to use when suggesting 
# aneuploidy at a chromosome or scaffold level when using `--smooth` flag; 0.5
# seems better but keeping at 0.25 for now to be sure to cast wide net
deviation_fraction = 0.51 #fraction of values for a chromosome or scaffold that 
# need to deviate from baseline by deviation_factor before aneuploidy is 
# suggested when using `--smooth` flag

plot_style = "seaborn" #try also `ggplot`,`default`, `bmh` or `grayscale`; use 
# `print(plt.style.available)` after appropriate imports to see others; 
# illustrated at https://matplotlib.org/examples/style_sheets/style_sheets_reference.html

#colors = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9']
# if using `ggplot` or `seaborn` style as coded in script than it has less than 
# 10 colors in the `CN` group and begins cycling internally back to `C0` when 
# hits `C7`(`C6` for seaborn), etc., & so same colors appear several times in 
# row when cycles again to `C0`. Options to fix current situation:
colors = (['C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'xkcd:magenta', 'xkcd:orange', 
    'tab:gray','tab:pink'])  # for use with `seaborn` style
colors = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6'] # for use with ggplot style
colors = (['C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'xkcd:orange','tab:pink', 
    'tab:cyan'])  # for use with `ggplot` style
colors = (['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 
    'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']) 
colors = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5'] # for use with seaborn style
colors = (['C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'xkcd:magenta'])# w/`seaborn`
# see https://matplotlib.org/users/colors.html for other options


#
#*******************************************************************************
#**********************END USER ADJUSTABLE VARIABLES****************************





















#*******************************************************************************
#*******************************************************************************
###DO NOT EDIT BELOW HERE - ENTER VALUES ABOVE###

import sys
import os
import argparse
import pandas as pd
import numpy as np
from itertools import cycle
import matplotlib # in order to use `matplotlib.use('Agg')`, need this first, see source of next line
matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend for running on PythonAnywhere. # from https://stackoverflow.com/questions/2801882/generating-a-png-with-matplotlib-when-display-is-undefined after searched error I was getting after upgrading matplolib in my Pythonanywhere account
import matplotlib.pyplot as plt
import seaborn as sns
from statsmodels.nonparametric.smoothers_lowess import lowess



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
    
    if `limit_to_chrs = ["II","IV"]`,
    Calling function with
        ("data1.txt", "_across_chr.png")
    returns
        "data1_across_chr_II_IV.png"
    '''
    main_part_of_name, file_extension = os.path.splitext(
        file_name) #from http://stackoverflow.com/questions/541390/extracting-extension-from-filename-in-python
    if limit_to_chrs:
        main_part_of_suffix, suffix_file_extension = os.path.splitext(
        suffix)
        suffix = main_part_of_suffix + "_" + "_".join(
            limit_to_chrs) + suffix_file_extension 
    return main_part_of_name + suffix



def extract_gene_ids(row):
    '''
    parses out the gene_id from the attributes list in a line in annoation file 
    formatted in Ensembl format, for example Ensembl-formatted yeast genome from 
    ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Saccharomyces_cerevisiae/Ensembl/R64-1-1/Saccharomyces_cerevisiae_Ensembl_R64-1-1.tar.gz
    Here the annotation file has been read in as a pandas dataframe and the 
    `attributes` field is in column 9, specified by `row[8]` in the code.
    Using that instead of the column name because not called attribures list
    in the gff format but it has `group` field at same position.

    takes a row of a pandas dataframe derived from an annotation file

    returns a string of the systematic gene id
    '''
    return row[8].split("gene_id")[1].split('"')[1].strip()

def calculate_position(row):
    '''
    takes a row of pandas dataframe derived from an annotaion file & calculates 
    the average position of a gene or feature along the chromosome by taking the
    average of the start and end

    returns an integer
    '''
    return (int(row["start"]) + int(row["end"]))/2



def checkIfRomanNumeral(numeral):
    """
    Controls that the userinput only contains valid roman numerals
    function from praveen's answer at 
    https://stackoverflow.com/questions/20973546/check-if-an-input-is-a-valid-roman-numeral
    """
    numeral = numeral.upper()
    valid_roman_numerals = ["M", "D", "C", "L", "X", "V", "I", "(", ")"]
    valid = True
    for letters in numeral:
        if letters not in valid_roman_numerals:
            #print("Sorry that is not a valid roman numeral")
            valid = False
            break
    return valid

def int_to_roman(input):
    """
    from Paul Winkler's http://code.activestate.com/recipes/81611-roman-numerals/
    (had to reindent; was causing text editor to default to wrong spacing
    otherwise. Plus updated some idioms too.)
    Convert an integer to Roman numerals.

    Examples:
    >>> int_to_roman(0)
    Traceback (most recent call last):
    ValueError: Argument must be between 1 and 3999

    >>> int_to_roman(-1)
    Traceback (most recent call last):
    ValueError: Argument must be between 1 and 3999

    >>> int_to_roman(1.5)
    Traceback (most recent call last):
    TypeError: expected integer, got <type 'float'>

    >>> for i in range(1, 21): print int_to_roman(i)
    ...
    I
    II
    III
    IV
    V
    VI
    VII
    VIII
    IX
    X
    XI
    XII
    XIII
    XIV
    XV
    XVI
    XVII
    XVIII
    XIX
    XX
    >>> print int_to_roman(2000)
    MM
    >>> print int_to_roman(1999)
    MCMXCIX
    """
    if type(input) != type(1):
        raise TypeError ("expected integer, got {0}".format(type(input)))
    if not 0 < input < 4000:
        raise ValueError ("Argument must be between 1 and 3999")
    ints = (1000, 900,  500, 400, 100,  90, 50,  40, 10,  9,   5,  4,   1)
    nums = ('M',  'CM', 'D', 'CD','C', 'XC','L','XL','X','IX','V','IV','I')
    result = ""
    for i,e in enumerate(ints):       #formerly `for i in range(len(ints))`
        count = int(input / e) #formerly `count = int(input / ints[i])`
        result += nums[i] * count
        input -= e * count #formerly `input -= ints[i] * count`
    return result

def roman_to_int_if_possible(input):
    """
    modified from Paul Winkler's roman_to_int at
    http://code.activestate.com/recipes/81611-roman-numerals/
    (had to reindent; was causing text editor to default to wrong spacing
    otherwise. Plus updated some idioms too.)

    Try to convert a roman numeral to an integer. 
    Return original input if not possible.
    """
    #if type(input) != type(""):
    #   raise TypeError, "expected string, got %s" % type(input)
    input = input.upper()
    nums = ['M', 'D', 'C', 'L', 'X', 'V', 'I']
    ints = [1000, 500, 100, 50,  10,  5,   1]
    places = []
    for c in input:
        if not c in nums:
            #raise ValueError, "input is not a valid roman numeral: %s" % input
            return input
    for i,e in enumerate(input):  #formerly `for i in range(len(input)):`
        c = e #formerly `c = input[i]`
        value = ints[nums.index(c)]
        # If the next place holds a larger number, this value is negative.
        try:
            nextvalue = ints[nums.index(input[i +1])]
            if nextvalue > value:
                value *= -1
        except IndexError as e:
            # there is no next place.
            pass
        places.append(value)
    sum = 0
    for n in places: sum += n
    # Easiest test for validity...
    if int_to_roman(sum) == input:
        return sum
    else:
        #raise ValueError, 'input is not a valid roman numeral: %s' % input
        return input
    
def seqname_roman_to_numeric(row):
    '''
    takes a row of pandas dataframe derived from an annotaion file and attempts 
    to convert the seqname value in roman numeral form to a numeric

    returns an integer if able to convert; 
    returns original value if unable to convert
    '''
    return roman_to_int_if_possible(row["seqname"])


def seqname_string_to_numeric(row):
    '''
    takes a row of pandas dataframe derived from an annotaion file and attempts 
    to convert the seqname value in string form to a numeric

    returns an integer if able to convert; 
    returns original value if unable to convert
    '''
    try:
        return int(row["seqname"])
    except ValueError as e:
        return row["seqname"]  

def deviates_from_baseline(values):
    '''
    Determines if values for a chromosome or scaffold suggest aneuploidy.

    Specifically, the function assess if a substantial fraction of lowess ys 
    for a chromosome differ from baseline, i.e., are a good fraction greater 
    than 0.25 or below 0.25 whenlog2 is used (default) or greater than 1.25 and
    less than 0.75 when no_log.
    This degree can be adjusted with deviation_factor setting under 
    "USER ADJUSTABLE VALUES". Default setting provided in script is 0.25.

    Returns a boolean value indicating if chromosome or scaffold shows 
    aneuploidy based on a simplistic assessment. Your mileage may vary.
    '''
    baseline = 1.0 if no_log else 0.0
    return float(len([value for value in values if not (baseline
         - deviation_factor < value < baseline
         + deviation_factor)])) / float(len(values)) \
    > deviation_fraction

###--------------------------END OF HELPER FUNCTIONS---------------------------###
###--------------------------END OF HELPER FUNCTIONS---------------------------###














#*******************************************************************************
###-----------------for parsing command line arguments-----------------------###
parser = argparse.ArgumentParser(prog='plot_expression_across_chromosomes.py',
    description="plot_expression_across_chromosomes.py  plots a ratio of \
    expression values across chromosomes or scaffolds of a genome to highlight \
    regions of deviation. Besides the options listed here, there are several \
    `USER ADJUSTABLE VALUES` inside the script that can be edited for easy \
    customization. A similar plot is called a Manhattan plot and this \
    implementation borrows the plotting approach and some of the features from \
    Brent Pedersen's awesome `manhattan-plot.py` script.       \
    **** Script by Wayne Decatur   \
    (fomightez @ github) ***")

parser.add_argument("annotation", help="Name of file containing the genome \
    annotation. REQUIRED. This is needed to determine the order of individual \
    data points along the chromosome and how to display the data across \
    chromosomes or scaffolds.", 
    type=argparse.FileType('r'), metavar="ANNOTATION_FILE")
parser.add_argument("data", help="Name of file containing the summarized data \
    to plot, such as mean TPM or RPKM, etc. in tab-delimited form. REQUIRED. \
    See my script `plot_expression_across_chromosomes_from_raw.py` if you want \
    supply the individual `raw` data files with the level metric for each \
    sample and/or replicate.", 
    type=argparse.FileType('r'), metavar="DATA_FILE")
parser.add_argument('-cols', '--columns', action='store', type=str, 
    default= '1,2,3', help="columns for gene, wild-type (baseline state) \
    expression value, experimental condition expression value, in that order. \
    This flag is used to specify the data in the summary file to be plotted. \
    Separate the column identifiers by commas, without spaces. \
    Default is `1,2,3`, where `1` indicates the first column, i.e., how \
    you'd refer to the columns in natural language (no zero-indexing). ") 
    # based on
    # https://stackoverflow.com/questions/15753701/argparse-option-for-passing-a-list-as-option
parser.add_argument("-l", "--lines",help=
    "add this flag to plot the expression level ratio value as lines \
    extending from the x-axis rather than points in space. (The resulting \
    aesthetic may resemble a city skyline for which the `manhattan plot` is \
    named.)",
    action="store_true")
parser.add_argument('-chr', '--chrs', action='store', type=str, 
    help="use this flag to limit plotting of the data to particular \
    chromosomes or scaffolds you specify immediately following this flag. \
    Separate the chromosome or scaffold identifiers by commas, without spaces. \
    Example use in a command is `--chrs I,IV,XVI`. \
    Default when this optional flag is not called is to plot that data for all \
    chromosomes or scaffolds. ") # based on
    # https://stackoverflow.com/questions/15753701/argparse-option-for-passing-a-list-as-option
parser.add_argument("-nl", "--no_log",help=
    "add this flag to keep the expression level ratio to be plotted in the \
    common base 10 instead of converting to log2.",
    action="store_true")
parser.add_argument("-nlim", "--no_limits",help=
    "add this flag to not impose a limit of above and below {} in plot window \
    when converting to log2. The cutoff can also be adjusted under \
    `user-adjustable settings` in the script. Issuing this flag has no effect \
     if all values are within +/- the cutoff interval or `--no_log` is used."
    .format(y_cutoff),
    action="store_true")
parser.add_argument("-s", "--smooth",help=
    "add this flag to display a smoothing curve fit to the data points \
    (LOWESS) on a per chromosome basis. This option can enhance visualization \
    of deviations characteristic of aneuploidy and copy number variation across \
    the genome, both within and between chromosomes. Additionally, a \
    simplistically-based assesment will be made for aneuploidy at the \
    chromosome or scaffold level and a notice will be made as the program is \
    running if aneuploidy at the chromosome or scaffold level seems indicated \
    by this simple metric. Further examination is warranted regardless of \
    the result this automated assessment.",
    action="store_true")
parser.add_argument('-ed', '--exp_desig', action='store', type=str, 
    default= 'experimental', help="Allows changing the text used in y-axis \
    label to reference experimental sample. Following `--exp_desig` type what \
    you'd like to read there instead of `experimental`.") 
parser.add_argument('-bd', '--base_desig', action='store', type=str, 
    default= 'wild\mathrm{-}type', help="Allows changing the text used in y-axis \
    label to reference wild-type or baseline sample. Following `--base_desig` \
    type what you'd like to read there instead of `wild-type`.") 
parser.add_argument("-svg", "--save_vg",help=
    "add this flag to save as vector graphics \
    (**RECOMMENDED FOR PUBLICATION***) instead of default png. Not default or \
    saved alongside default because file size can get large due to the large \
    number of points.",
    action="store_true")
parser.add_argument("-ndh", "--no_data_header",help=
    "add this flag if there is no data header or no first line of column names \
    in the data file. Otherwise, it is assumed there is and any item read as \
    the first gene identifier from the first line won't be highlighted as \
    missing from annotation.\
    IMPORTANTLY, this only affects feedback provided as script is run. If the \
    first line resembles data, i.e., numbers in specified columns, it will be \
    automagically parsed as if data. Remove the header or column labels line \
    from your summary data file on the off-chance this causes issues in your \
    resulting plot.",
    action="store_true")
parser.add_argument('-ac', '--advance_color', action='store', type=int, 
    default= '0', help="**FOR ADVANCED USE.*** Allows for advancing the color \
    selection iterator the specified number of times. The idea is it allows \
    the ability to control the color of the chromosome when specifying \
    a chromosome or scaffolds to plot so you could make the color match the \
    one used when all chromsome plotted if needed. Supply the number to \
    advance after the flag on the command line. For example, `-ac 4`.") 

#I would also like trigger help to display if no arguments provided because need at least one input file
if len(sys.argv)==1:    #from http://stackoverflow.com/questions/4042452/display-help-message-with-python-argparse-when-script-is-called-without-any-argu
    parser.print_help()
    sys.exit(1)
args = parser.parse_args()
annotaton_file = args.annotation
data_file = args.data
data_columns_to_grab = [int(item) for item in args.columns.split(',')]
no_log = args.no_log
no_data_header = args.no_data_header
lines = args.lines
if args.chrs:
    if "," in args.chrs:
        limit_to_chrs = args.chrs.split(',')
    else:
        # means only one item
        limit_to_chrs = [args.chrs] #has to be a list for passing to Pandas `isin()` 
else:
    limit_to_chrs = args.chrs # will make `limit_to_chrs` as `None`
advance_color_increments = args.advance_color
display_smooth = args.smooth
no_limits = args.no_limits
exp_designation = args.exp_desig
baseline_designation = args.base_desig
save_vg = args.save_vg





###-----------------Actual Main portion of script---------------------------###

# ANNOTATION FILE ACCESSING AND GENOME DATAFRAME INITIAL PREPARATION
# open annotation file and make it a Pandas dataframe
sys.stderr.write("\nReading annotation file and getting data on genes and "
    "chromosomes...")
#determine if annotation is gff or gtf
if "gff" in annotaton_file.name.lower():
    col_names_to_apply = genome_annotation_fields_for_gff 
if "gtf" in annotaton_file.name.lower():
    col_names_to_apply = genome_annotation_fields_for_gtf 

# read in annotation file
init_genome_df = pd.read_csv(
    annotaton_file, sep='\t', header=None, low_memory=False,
    names=col_names_to_apply, comment='#') # comment handling added because I 
# came across gtfs with a header that had `#`  at start of each line. Others 
# must have encountered same because I saw someone dealing with it at 
# https://github.com/shenlab-sinai/ngsplotdb/pull/2/files. I cannot use that 
# solution since I use Pandas read_csv function. (`read_table`) was 
# deprecated recently.

# parse out gene_ids from attribute or group, i.e., 9th column in the annotation file
init_genome_df["gene_id"] = init_genome_df.apply(extract_gene_ids, axis=1)

# copy each row to a new dataframe, unless gene already present. 
# This will give me unique gene_ids for each and I can make that index.
# Because it takes first occurence of each gene, it only has that as start and
# end. Then to get the full range of data for start and end for each gene, I can
# subset the initial dataframe on each gene_id and get the min and max and use
# those values to replace `start` and `end` for the new dataframe. For those
# on Crick strand, it will turn around the start, end information, but that is
# fine since just want an avg relative position and don't care about direction.
genome_df = pd.DataFrame(columns=init_genome_df.columns)
for i, row in init_genome_df.iterrows():
    if not any(genome_df.gene_id == row.gene_id):
        genome_df = genome_df.append(row)
genome_df = genome_df.set_index('gene_id')
for id in list(genome_df.index.values):
    sub_df = init_genome_df.loc[init_genome_df["gene_id"] == id]
    min_val = min(sub_df[["start","end"]].min()) # `sub_df["start","end"].min()` gives values for the two columns and so min() of that gives single value
    max_val = max(sub_df[["start","end"]].max())
    genome_df.loc[id, "start"] = min_val
    genome_df.loc[id, "end"] = max_val
# provide feedback on number of unique genes identified
sys.stderr.write("Information for {0} genes parsed...".format(len(genome_df)))
# calculate average position
# genome_df["position"] = genome_df[["start","end"]].apply(np.mean, axis=1) # gives float and I'd prefer as integer
genome_df["position"] = genome_df.apply(calculate_position, axis=1)
# make a column of chrosomes as numbers, either converting from the string they
# would be by default (I think?) or converting from roman numerals. This will be
# used for sorting later
# First determine if chromosomes are numbers with X and Y (and others?) or if
# in the form of roman numerals
chromosomes_in_roman_num = False
# check if majority look like integers. If that is the case verify most are 
# valid roman numerals just to be sure.
genome_df['seqname'] = genome_df['seqname'].apply(str) # cast to string so 
# string methods like `.isdigit()` continue to work even when chromosomes are 
# actually numbers; faster way from 
# https://stackoverflow.com/questions/17950374/converting-a-column-within-pandas-dataframe-from-int-to-string/44008334
seqname_set = set(genome_df['seqname'].tolist())
chromosomes_in_roman_num = not bool(len([s for s in seqname_set if s.isdigit()]) > len(seqname_set)/2) #checks if most chromosomes seem to be digits and says they are roman numerals if that is false
if chromosomes_in_roman_num:
    most_valid_rom_numerals = bool(len([n for n in seqname_set if checkIfRomanNumeral(n)]) > len(seqname_set)/2)
    # provide feedback if don't seem to be roman numerals
    if not most_valid_rom_numerals:
        sys.stderr.write("***WARNING***The chromosomes seem to not be in numeric form, but most aren't valid roman numerals either.***WARNING***....")
    else:
        sys.stderr.write("The chromosomes appear to be in roman numeral form....")
else:
    sys.stderr.write("The chromosomes appear to be in numeric form....")
# Second, convert if `chromosomes_in_roman_num` otherwise just try to coerce string to integer, allowing not to coerce and not throw an error in case of the sex chromosomes
if chromosomes_in_roman_num:
    # try and convert each row but allow for non type change without error
    genome_df["chr_as_numeric"] = genome_df.apply(
        seqname_roman_to_numeric, axis=1)
else:
    genome_df["chr_as_numeric"] = genome_df.apply(
        seqname_string_to_numeric, axis=1)
longest_chr_or_scaffold = len(max(seqname_set, key=len))
# Prepare genome dataframe for accessing data by adding a column for relative 
# level and fill with 'NaN'
genome_df["level_val"] = np.nan


# ACCESSING EXPRESSION DATA
# open data file and calculate and add the relative information for each gene 
# to the Pandas dataframe.
# make a list for any genes not found in the genome_df to report to user later
# make a list for any genes triggering `ZeroDivisionError: float division by 
# zero`; note that I expanded this list to include those where experimental 
# level was `0.0`, unless using `no_log` flag, because cases like that where
# result was `0.` when converted to numpy array was causing divide by zero 
# issues when converting to log2 later. I'm trying to allow it to stay untouched
# if using the `no_log` flag in order to impose as little constraints as 
# possible to keep script more broadly useful.
sys.stderr.write("Parsing data file...")
not_found_in_annotation = []
cause_div_by_zero = []
lines_processed = 0
data_entered = 0
for line in data_file:
    lines_processed += 1
    items = line.split("\t")
    gene = items[data_columns_to_grab[0]-1] #first entry in data_columns_to_grab should be column for gene identifier; -1 to account for zero indexing of python
    if gene in genome_df.index:
        if float(items[data_columns_to_grab[1]-1]) != 0.0 and ((not no_log and float(items[data_columns_to_grab[2]-1]) != 0.0) or (no_log)):
            genome_df.loc[gene, "level_val"] = float(items[data_columns_to_grab[2]-1])/float(items[data_columns_to_grab[1]-1])
            data_entered += 1
        else:
            cause_div_by_zero.append(gene)
    else:
        if lines_processed == 1 and not no_data_header: continue #the entry for gene on first line isn't found because there is a header line or column names, and so don't collect it in `not_found_in_annotation`
        not_found_in_annotation.append(gene)
# provide general feedback
sys.stderr.write("{0} lines of data parsed and relative levels for {1} genes recorded...".format(lines_processed, data_entered))
# provide feedback on any genes from the data with no correspondence in the annotation
if not_found_in_annotation:
    sys.stderr.write("\nNote that data for {} genes was read but no corresponding annotation entry was found. Those genes are {!r}.".format(
        len(not_found_in_annotation), not_found_in_annotation))
    if data_entered != 0 and ((len(not_found_in_annotation)/data_entered)*100) > 3:
        sys.stderr.write(" Since that is less than 3% of the genes with expression data and so probably won't make much difference to the plot...")
    else:
        sys.stderr.write("***WARNING***...")
# provide feedback on any genes causing zero division error
if cause_div_by_zero:
    sys.stderr.write("\nNote that the value in the data for {} genes would trigger division by zero errors and so they were left out. Those genes are {!r}...".format(
        len(cause_div_by_zero), cause_div_by_zero))





# PREPARE GENOME DATAFRAME NOW HARBORING EXPRESSION DATA FOR PLOTTING
# if limiting to specific chromosomes or scaffolds, discard any other data
if limit_to_chrs:
    genome_df = genome_df.loc[genome_df["seqname"].isin(limit_to_chrs)]
# sort by position on chromosome
genome_df.sort_values(["position"], inplace=True, ascending=True) 
# sort by chromosome
genome_df.sort_values(["chr_as_numeric","position"], inplace=True, ascending=True) # found `genome_df.sort_values(["chr_as_numeric"], inplace=True, ascending=True)` alone after previous sort didn't work but this seemed to result in what I wanted in end
# make a dictionary of dictionaries with details of each chromsomes, 
# specifically length(approximate based on genes/features) and midpoint
# (to be used for tick marks later).
# not all assigned variables for the chromosomes_specs dictionary used here but
# also useful for generating a summary of assignments used & is small and could 
# be more useful down the road.
chr_specs = {}
xs_by_chr = [] #ported over from https://github.com/brentp/bio-playground/blob/master/plots/manhattan-plot.py
# sorting just above makes it easy to get last entry for each chromosome
grouped = genome_df.groupby('seqname')
previous_chr_last_x = 0
for chr, data_per_chr_df in grouped:
    chr_specs[chr] = {}
    chr_length = data_per_chr_df['end'].tolist()[-1]
    chr_midpoint = chr_length/2
    chr_specs[chr]["length"] = chr_length
    chr_specs[chr]["midpoint"] = chr_midpoint
    # next ones just for plotting adjacent the previous chromosome
    chr_specs[chr]["x_start"] = previous_chr_last_x
    x_midpoint = previous_chr_last_x + chr_midpoint
    chr_specs[chr]["x_midpoint"] = x_midpoint
    xs_by_chr.append((chr,x_midpoint))
    previous_chr_last_x = previous_chr_last_x + chr_length
#I THOUGHT THIS WOULD RESULT IN ORDER CURRENTLY IN DATAFRAME BUT IT DOESN'T. IF WERE ABLE TO HAVE IT WORK, IS IT EVEN NECESARY??#sorted_chr_set = set(genome_df['seqname'].tolist()) #should be order as established in SORTED genome_df; earlier set wasn't from the version sorted on chromosome

# filter out anything where "level_val" is `NaN`; waited to do this until after
# making the dictionary of chr_specs so could consider even genes or features 
# where no expression plotted as contributing to position information to keep
# width scale in plot closer to reality for each chromosome.
genome_df.dropna(subset=["level_val"], inplace = True)
# group by chromosomes and iterate making lists of x values that increase across
# all chromosomes and fills in corresponding y values and arranges colors.
# Note even though I iterated on a similar `grouped` by chromosome iterable 
# earlier, I couldn't yet capture x & y values for plot b/c had not yet removed 
# those with no expression values. Not DRY, but will better represent chromosome
# relative scale this way.
grouped_from_filtered = genome_df.groupby('seqname')
xs = []
ys = []
cs = []

colors = cycle(colors)
# handle advancing colors if adjusting based on user provided info
if advance_color_increments:
    for i in range(advance_color_increments): color = next(colors)
# prepare to collect data on per chromosome basis if lowess smoothing flag on
if display_smooth:
    data_by_chr = {}
previous_chr_last_x = 0
for chr, data_per_chr_df in grouped_from_filtered:
    chr_length = chr_specs[chr]["length"] # this may be larger than the final 
    # "end" position in the filtered group, & so would better reflect scale of 
    # each chromosome
    color = next(colors)
    this_chromosomes_xs = [previous_chr_last_x + position for position in data_per_chr_df['position'].tolist()]
    this_chromosomes_ys = [level_val for level_val in data_per_chr_df['level_val'].tolist()]
    assert min(this_chromosomes_xs) >= chr_specs[chr]["x_start"], (
    "x values for individual genes cannot be lower than the \
    `starting x value` for that chromosome.")
    assert max(this_chromosomes_xs) <= (
        chr_specs[chr]["x_start"] + chr_specs[chr]["length"]), ("x values for \
        individual genes cannot be greater than the sum of the `starting x \
        value` for that chromosome plus the length of that chromosome.")
    xs.extend(this_chromosomes_xs)
    ys.extend(this_chromosomes_ys)
    cs.extend([color] * len(data_per_chr_df))
    # collect pertinent data on per chromosome basis if lowess smoothing flag
    if display_smooth:
        data_by_chr[chr] = {}
        data_by_chr[chr]["xs"] = this_chromosomes_xs
        data_by_chr[chr]["ys"] = this_chromosomes_ys
    previous_chr_last_x = previous_chr_last_x + chr_length





# MAKE THE PLOT
# This follows the plotting approach of Brent Pedersen's `manhattan-plot.py`  
# closely, see https://github.com/brentp/bio-playground/blob/master/plots/manhattan-plot.py .

xs = np.array(xs)
ys = np.array(ys) if no_log else np.log2(ys)

plt.close()
plt.style.use(plot_style)
f = plt.figure()
#ax = f.add_axes((0.1, 0.09, 0.88, 0.85))
ax = plt.axes()

if title_prefix is not None:
    if limit_to_chrs:
        title = title_prefix + " ".join(limit_to_chrs)
    else:
        title = title_prefix + "genome"
    plt.title(title, fontsize=18)

# At first I could not get LATeX to work with the y-axis labels even though what I have in the commented out lines below between `$` works in both Jupyter notebooks markdown cells (put between `$$`) and inserted in place of the LATeX code from the tex_demo.py script at https://matplotlib.org/users/usetex.html that use `plt.plot()`
'''
if no_log:
    ax.set_ylabel('$\frac{\text{experimental\ level}}{\text{wild-type\ level}}$')
else:
    ax.set_ylabel('$\log 2\left(\frac{\text{experimental\ level}}{\text{wild-type\ level}}\right)$',
          fontsize=16)
'''
# Then I found https://stackoverflow.com/questions/23824687/text-does-not-work-in-a-matplotlib-label says that the issue is that matplotlib uses a different rendering engine, MathText. Seems for `\text{}` want obscure `\mathrm{}`. Plus if concatenating variable, need `r` again in front of string after the final plus sign in a concatenation command.
if no_log:
    ax.set_ylabel(r'$\frac{\mathrm{'+ exp_designation +r'\ level}}{\mathrm{'+ baseline_designation +r'\ level}}$',
          fontsize=16)
else:
    ax.set_ylabel(r'$\log 2 \left(\frac{\mathrm{'+ exp_designation +r'\ level}}{\mathrm{'+ baseline_designation +r'\ level}}\right)$',
          fontsize=16)
# Below keeps ylabels plain for if LATeX failing
'''
if no_log:
    ax.set_ylabel('experimental level/wild-type level')
else:
    ax.set_ylabel('log2(experimental level/wild-type level)')
'''
if lines:
    ax.vlines(xs, [0], ys, colors=cs, alpha=0.5) # Despite change noted here still doesn't work when log2 used because of the `divide by zero encountered in log2` issue that makes infinites===> `0` to `[0]` based on example at http://matplotlib.org/1.2.1/examples/pylab_examples/vline_demo.html
else:
    ax.scatter(xs, ys, s=12, c=cs, alpha=0.8, edgecolors='none')

# plot 0.05 line after multiple testing.
#ax.axhline(y=-np.log10(0.05 / len(genome_df)), color='0.5', linewidth=2)
plt.axis('tight')
#plt.xlim(0, xs[-1])
if no_log or no_limits or all(-y_cutoff <= n <= y_cutoff for n in ys):
    #plt.ylim(ymin= None, ymax= None) # necessary?
    pass # see comment on above line
else: 
    plt.ylim(ymin=-y_cutoff, ymax=y_cutoff)
    # catch out of bounds points for plotting in a manner to indicate out of 
    # bounds of limits used to avoid compressing important range. This approach
    # is styled on how DESeq2 plotMA handles out of bounds points.  Found to 
    # easily control rotation of the marker in plot call, best to split.
    # prepare to catch points too high
    xs_out_of_bounds_hi = []
    ys_out_of_bounds_hi = []
    cs_out_of_bounds_hi = []
    # prepare catch points too low
    xs_out_of_bounds_lo = []
    ys_out_of_bounds_lo = []
    cs_out_of_bounds_lo = []
    offset = 0.055 #so not cut off at edge
    for indx, pt in enumerate(ys):
        if pt > y_cutoff:
            ys_out_of_bounds_hi.append(y_cutoff - offset)
            xs_out_of_bounds_hi.append(xs[indx])
            cs_out_of_bounds_hi.append(cs[indx])
        elif pt < -y_cutoff:
            ys_out_of_bounds_lo.append(-y_cutoff + offset)
            xs_out_of_bounds_lo.append(xs[indx])
            cs_out_of_bounds_lo.append(cs[indx])
    oobs_marker_sz = 25 # out of bounds marker size
    oobs_width = 1 # out of bounds marker linewwidth
    rotation_for_hi = 0
    rotation_for_lo = 180
    # "Plot" the out-of-bounds points
    plt.scatter(xs_out_of_bounds_hi, ys_out_of_bounds_hi, marker=(3, 0, rotation_for_hi), s=oobs_marker_sz, linestyle='solid',linewidth=oobs_width, facecolors='None', edgecolors=cs_out_of_bounds_hi) # see https://stackoverflow.com/questions/4143502/how-to-do-a-scatter-plot-with-empty-circles-in-python; use `facecolors` for scatter and markerfacecolor for `plt.plot`
    plt.scatter(xs_out_of_bounds_lo, ys_out_of_bounds_lo, marker=(3, 0, rotation_for_lo), s=oobs_marker_sz, linestyle='solid',linewidth=oobs_width, facecolors='None', edgecolors=cs_out_of_bounds_lo) # see https://stackoverflow.com/questions/4143502/how-to-do-a-scatter-plot-with-empty-circles-in-python; use `facecolors` for scatter and markerfacecolor for `plt.plot`
    sys.stderr.write(
        "\n***Notice:***The {} points beyond the bounds of y-axis are drawn as open triangles at the edge; a limit was imposed to avoid extreme values compressing the typically important range; run with `--no_limits` or `--no_log` to see these accurately."
        .format(len(xs_out_of_bounds_hi) + len(xs_out_of_bounds_lo)))
#if ymax is not None: plt.ylim(ymax=ymax)
size_for_xlabels = 8.5
longest_chr_or_scaffold = len(max(seqname_set, key=len))
if limit_to_chrs is None:
    if longest_chr_or_scaffold > limit_before_rotate:
        plt.xticks(
            [c[1] for c in xs_by_chr], [c[0] for c in xs_by_chr], 
            rotation=-90, size=size_for_xlabels)
    else:
        plt.xticks(
            [c[1] for c in xs_by_chr], [c[0] for c in xs_by_chr], size=size_for_xlabels)
else:
    # If only handling a few chromosomes or scaffolds don't bother checking if name exceeds the length. 
    # Deciding what is "few" there will be the same integer used to decide what length of characters
    # is too great to not rotate. This is done just to make it so another variable isn't needed.
    if (len(limit_to_chrs) > limit_before_rotate) and (longest_chr_or_scaffold > limit_before_rotate):
         plt.xticks(
            [c[1] for c in xs_by_chr], [c[0] for c in xs_by_chr], 
            rotation=-90, size=size_for_xlabels)
    else:
        plt.xticks(
            [c[1] for c in xs_by_chr], [c[0] for c in xs_by_chr], size=size_for_xlabels)



# add lowess curve fit to each chromosome if `smooth` flag
if display_smooth:
    chr_deviating_from_baseline = []
    for chr in data_by_chr:
        lowess_ys = lowess(
            data_by_chr[chr]["ys"], data_by_chr[chr]["xs"], 
            return_sorted=False) if no_log else lowess(
            np.log2(data_by_chr[chr]["ys"]), data_by_chr[chr]["xs"], 
            return_sorted=False)
        plt.plot(
            data_by_chr[chr]["xs"],lowess_ys,'gray',linewidth=12, alpha=0.55)
        # assess if chromosome/scaffold deviates from baseline
        if deviates_from_baseline(lowess_ys):
            chr_deviating_from_baseline.append(chr)
    # provide feedback if deviation from norm noted
    if chr_deviating_from_baseline:
        sys.stderr.write(
            "\nAneuploidy at the chromosome or scaffold level is suggested for {}; examine further.".format(" & ".join(chr_deviating_from_baseline)))


output_file_name = generate_output_file_name(data_file.name, suffix_for_saving_result)
if save_vg:
    plt.savefig(output_file_name[:-4]+".svg", orientation='landscape') # FOR VECTOR GRAPHICS; useful if merging into Adobe Illustrator. Based on https://neuroscience.telenczuk.pl/?p=331 ; I think ReportLab also outputs SVG?
    sys.stderr.write("\n\nPlot image saved to: {}\n".format(output_file_name[:-4]+".svg"))
else:
    plt.savefig(output_file_name)
    sys.stderr.write("\n\nPlot image saved to: {}\n".format(output_file_name))
    # plt.savefig(output_file_name[:-4]+".svg")
    # sys.stderr.write("\n\nPlot image saved to: {}\n".format(output_file_name[:-4]+".svg"))
    # plt.savefig(output_file_name[:-4]+".pdf", orientation='landscape') # UNFORTUNATELY DOES NOT PRODUCE VECTOR GRAPHICS, unlike ReportLab's pdf output; USE SVG for that and the make PDF later.
    #plt.show()


#*******************************************************************************
###-***********************END MAIN PORTION OF SCRIPT***********************-###
#*******************************************************************************
