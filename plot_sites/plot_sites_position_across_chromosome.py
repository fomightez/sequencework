#!/usr/bin/env python
# plot_sites_position_across_chromosome.py
__author__ = "Wayne Decatur" #fomightez on GitHub
__license__ = "MIT"
__version__ = "0.1.0"


# plot_sites_position_across_chromosome.py by Wayne Decatur
# ver 0.1
#
#*******************************************************************************
# Verified compatible with both Python 2.7 and Python 3.6; written initially in 
# Python 3. (See below.)
#
# PURPOSE: Takes a dataframe of 'sites', specified with starts and ends, & plots 
# the position/distribution of these 'sites' along the chromosome. Mainly made
# to highlight where different classes of these sites/elements occur
# relative the other classes across a chromsome.
#
# 'sites' is used as a general term for any sequence element(s) larger than a 
# single base that is specified with 'start' and 'end' data. Use it on your
# favorite class of genetic elements or sites, or many of such items. The 
# 'sites' specification is designed to allow multiple classes or types to be
# used and it expects they are named where the type specificer comes before a
# dash and a unique id number (integer) follows. The unique id number only need
# be unique among the class or type since it is only used as part of keeping 
# them distinguished, for example, when labeling a point on a plot. Valid 
# examples include:  `VY-3` or `MRT-18`, without the quotes/ticks.
#
# 
# Specifics needed in dataframes:
# Columns of starts and ends in provided dataframes must be named `start` and 
# `end`. `sys_gene_id` is to be used for the column name for 
# 'sites' identifiers. 
# This is the information needed. Extra columns in the dataframes don't matter.
# Use `df.rename()` to rename columns to match, see 
# https://gist.github.com/fomightez/ef57387b5d23106fabd4e02dab6819b4 .
#
# A number of options for the plot can be adjusted in 'USER ADJUSTABLE VALUES' 
# below. Some advanced or custom plot options require editing within the script
# itslef, for example to control coloring related to strand, etc. .
#
# Written to run from command line or pasted/loaded inside a Jupyter notebook 
# cell. 
# When using in a notebook, if you don't specify dataframe objects, you must
# instead supply strings of file names for the pickled dataframes in the call
# to the main function. 
#
#
#
#
#
# This script based on `plot_sites_position_relative_genes.py` and musings 
# developed in 
# `[REDACTED] relative genes normalized gene transcripts (and termini).md` 
#
#
#
#
#
# Verified compatible with both Python 2.7 and Python 3.6; written initially in 
# Python 3. (Compatibility check required pickled dataframe from Python 3 to be
# unpickled, saved as TSV/CSV, moved to Python 2 environment, read in and 
# re-pickled in Python 2 environment.)
#
#
# Dependencies beyond the mostly standard libraries/modules:
#
#
#
# VERSION HISTORY:
# v.0.1. basic working version
#
# To do:
# - make a demo notebook of use (including scatter with use of output flag 
    # and post to github
# - Add labels or interactive hovers to the individual points indicating sites?
    # Hovers will only work within Jupyter so make optional(?)
# - make so works with multiple chromosomes (Shouldn't be too hard. It is just 
# that for now I am working with one chromosome and not including the chromosome
# identifier as a part of the input dataframe detailing the genes along it.)
#
#
#
#
# TO RUN:
# Examples,
# Enter on the command line of your terminal, the line
#-----------------------------------
# python plot_sites_position_across_chromosome.py sites.pkl
#-----------------------------------
# Issue `plot_sites_position_across_chromosome.py -h` for details.
# 
#
# When using in a notebook, if you don't specify dataframe objects, , you must
# instead supply strings of file names for the pickled dataframes in the call
# to the main function. 
# To use this after pasting or loading into a cell in a Jupyter notebook, in
# the next cell specify the two dataframes then call the main function similar 
# to below:
# %matplotlib notebook
# sites_df = "df_for_merging.pkl"
# plot_save_as_name = "sites_across.png"
# p,df = plot_sites_position_across_chromosome(sites_df,return_df = True)
# p.figure.savefig(plot_save_as_name, dpi = (130)) #optional file save; best 
# #for getting larger version of plot
# df
#
#
#(`df_save_as_name` can be assigned in a cell before calling the function 
# as well in order to control name of pickled dataframe file. And 
# `sites_df = pd.read_pickle("df_for_merging.pkl") can be used in place of 
# `sites_df = "df_for_merging.pkl"`.)
#
# 
#
'''
CURRENT ACTUAL CODE FOR RUNNING/TESTING IN A NOTEBOOK WHEN LOADED OR PASTED IN 
ANOTHER CELL:
%matplotlib notebook
sites_df = "df_for_merging.pkl"
plot_save_as_name = "relative.png" # optional; only needed if saving plot
p,df = plot_sites_position_across_chromosome(sites_df, return_df = True)
p.figure.savefig(plot_save_as_name, dpi = (130)) #optional file save; best 
# for getting larger version of plot
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
plot_save_as_name = 'sites_across_plot.png' # name for plot to be saved
x_axis_label = "position"
plot_title = "sites across mitochondrial genome"
title_size = 12 # font size for title text

# Special output plot settings (typically don't get used/adjusted)
change_y_label_in_swarm = False # Set to true if you want to change y-axis label
# in default plot style (swarm). Assumed column name probably useable as y-axis 
# label in swarm plot, and so changing label optional/unnecessary.
y_axis_label_swarm = "NA" # y-axis label to apply in default plot style (swarm);
# the setting is only used if `change_y_label_in_swarm = True`.
y_axis_label_scatter = "arbitrary element member # (per type)" # Not typically 
# used; ONLY used if plotting `scatter`/`swarm=False` and not default style.    




#
#*******************************************************************************
#**********************END USER ADJUSTABLE VARIABLES****************************


















#*******************************************************************************
#*******************************************************************************
###DO NOT EDIT BELOW HERE - ENTER VALUES ABOVE###

import sys
import os
import pandas as pd
import seaborn as sns




###---------------------------HELPER FUNCTIONS---------------------------------###

def extract_type(item):
    '''
    Takes a site identifier, such as `M2-41`, and returns just the type/class,
    which would be `M2''` in this example, i.e., everything in front of the dash.
    '''
    return item.split("-")[0]


def midpoint(items):
    '''
    takes a iterable of items and returns the midpoint (integer) of the first 
    and second values
    '''
    return int((int(items[0])+int(items[1]))/2)






###--------------------------END OF HELPER FUNCTIONS---------------------------###
###--------------------------END OF HELPER FUNCTIONS---------------------------###

#*******************************************************************************
###------------------------'main' function of script---------------------------##

def plot_sites_position_across_chromosome(
    sites_df, return_df = False, swarm=True):
    '''
    Main function of script. 
    It will take a string of a file name for the pickled dataframe or dataframe 
    object directly, the latter being intended for use of the function in a 
    Jupyter notebook. The former for use of the script on the command line.

    Default is to make a swarmplot where handles not letting points with that 
    would overlap (same or close x value) occlude each other automatically, but 
    originally worked out where plot individual point y-values based on 
    arbitrary member numbers in group as scatter plot and may be useful for 
    interactively labeling or if want to illustrate how patterns emerge better 
    with swarmplot. Plus the scatter plot makes a nice legend that could be 
    possibly be composited to the side of swarmplot since I have not figured out 
    if/how to make one with seaborn's swarmplot. Plus, I worked out a lot 
    `seaborn.lmplot()` settings that could be useful elsewhere. 

    Returns a plot object like 
    `<matplotlib.axes._subplots.AxesSubplot at 0x7f3aa8763128>` that as of now,
    I don't know what to do with when I am in a Jupyter notebook. Only approach 
    I have found to putting a good-sized version of the plot in a notebook is 
    to use `%matplotlib notebook` at start of cell, otherwise it plot is small.
    Optionally, returns a dataframe of the sites with details used in plotting
    added.
    '''
    # Bring in the necessary data:
    #---------------------------------------------------------------------------

    # if strings provided, assume file names and get the pickled dataframes
    from six import string_types #use based on https://stackoverflow.com/a/11301392/8508004
    if isinstance(sites_df, string_types):
        sites_df = pd.read_pickle(sites_df)
    # Otherwise assume it is a dataframe object provided.

    # feedback
    sys.stderr.write("Provided data read...")

    sys.stderr.write("\nArranging data for plotting...")
    # Additional preparation:
    #---------------------------------------------------------------------------

    # Make sure the dataframe has midpoints calculated. To insure, without 
    # needing to come up with fool-proof method to detect and insure correct,
    # best to calculate & give column name unqie identifier not likely to have 
    # been used by other users.
    updated_sites_df = sites_df.copy()
    updated_sites_df["calcd_mdpt"] = sites_df[['start','end']].apply(
        midpoint, axis=1)

    
    #add type column to updated_sites_df. This will be everything in front of 
    # the dash for the site identifiers.
    updated_sites_df["class"] = updated_sites_df['sys_gene_id'].apply(extract_type)

    # add a column that is the number of the site among its class. 
    # For example the second `Ma3` element, should get 2, etc.
    # I am going to do this a rather brutforce way using the data from the
    # 'class' column because I know that will work, rather
    # than figure out how to do this with `groupby`.
    # (Note that although I my site elements already have numbers that specify
    # each within its group, I am adding this to insure it works for even when
    # it isn't already present and numbers maybe just sequential, irregardless
    # of class. Also note, that my groups where separated already as blocks
    # but this approach will work if interleaved, whereas the way I added
    # unique sub-class numbers early in
    # `comparing_[REDACTED]_to_features_read_signal` relied on them being in
    # blocks as a result of coming from BLAST. This method works if intermixed.)
    updated_sites_df=updated_sites_df.rename(columns={'class':'class_'})# change
    # name of column temporarily because `row.class` throws errow since `class`
    # reserved word
    from collections import defaultdict
    tracking_dict = defaultdict(int)
    type_encounters_hitherto = []
    for row in updated_sites_df.itertuples():
        # print(row) # for debugging
        tracking_dict[row.class_]+= 1
        type_encounters_hitherto.append(tracking_dict[row.class_])
    updated_sites_df["class_element_num"] = type_encounters_hitherto
    updated_sites_df=updated_sites_df.rename(columns={'class_':'class'})#revert
    
    #print(updated_sites_df) # for debugging


    # build the plot object
    if swarm:
        p = sns.swarmplot(x="calcd_mdpt", y="class", data=updated_sites_df)

        # -OR- highlight strand with these, in turn:
        # p = sns.swarmplot(x="calcd_mdpt", y="class", hue="strand", data=updated_sites_df, palette="seismic")

        # -OR- if want overlap of points, use `seaborn.stripplot`
        #p = sns.stripplot(x="calcd_mdpt", y="class", data=updated_sites_df)

        # Alternatively, I found via example at 
        # https://elitedatascience.com/python-seaborn-tutorial 
        # where they had `split` which is now `dodge` that I could 
        # keep something like the class coloring and highlight by strand 
        # stratified with:
        #updated_sites_df['strand'] = updated_sites_df.strand.astype('category')
        #p = sns.swarmplot(x="calcd_mdpt", y="strand", hue="class", dodge=True, data=updated_sites_df)
        #p.legend(bbox_to_anchor=(1, 1), loc=2) # otherwise legend was coming 
        # up in middle. When I later when to look at file made using this, I 
        # could see right part of legend was cut off even thouh the tiny version
        # I saved from the cell by right-clicking and choosing 'Save Image As..' 
        # after running looked fine, with no part cut off.



        # Fix the axes labels; most likely only needed for x-axis
        p.set_xlabel(x_axis_label) # based on something more like 
        # https://stackoverflow.com/a/46235777/8508004 than 
        # https://stackoverflow.com/a/36573970/8508004
        # because swarmplot doesn't seem to be a FacetGrid like `seaborn.lmplot`
        # Because column name probably useable as y-axis label in swarm plot, 
        # make changing label optional.
        if change_y_label_in_swarm:
            p.set_ylabel(y_axis_label_swarm)

        # Possible to turn off the y-axis line?
        sns.despine(left=True) # from https://seaborn.pydata.org/tutorial/aesthetics.html

        # Add a title to the plot
        p.set_title(plot_title, fontsize=title_size) #also works with more basic
        # seaborn approach (like https://stackoverflow.com/a/32724156/8508004 )
        # because swarmplot doesn't seem to be a FacetGrid like `seaborn.lmplot`




        # Add labels or hovers to the individual points indicating sites?
        # hovers will only work within Jupyter so make optional(?)
    else:
        '''#first test visualization:
        p = sns.lmplot(
            x="position", y="calcd_mdpt", hue="class", 
            data=updated_sites_df, fit_reg=False)
        '''
        p = sns.lmplot(
            x="calcd_mdpt", y="class_element_num", hue="class", 
            data=updated_sites_df, fit_reg=False)
        '''
        # one option that splits things out
        #sns.set() # from https://seaborn.pydata.org/tutorial/aesthetics.html
        sns.set_style("darkgrid")
        sns.set_context("talk")
        p = sns.factorplot(x="calcd_mdpt", y="class",
                           #hue="smoker", col="time",
                           data=updated_sites_df, kind="swarm",);
        '''

        # Fix the axes labels
        p.axes[0,0].set_xlabel(x_axis_label) # based on 
        # https://stackoverflow.com/a/36573970/8508004
        # because seaborn.lmplot is actually a FacetGrid
        p.axes[0,0].set_ylabel(y_axis_label_scatter)

        # Because the 'numbers' involved are arbitrary and just added so those 
        # that overlap on x-axis don't occlude each other...
        # Turn off the y-axis ticks
        p.set(yticks=[]) # based on https://stackoverflow.com/a/24497489/8508004
        # I also found that `p.axes[0,0].set(yticks=[])` based on editing 
        # axes labels also works too.

        # Add a title to the plot
        # Based on https://stackoverflow.com/a/49389150/8508004
        # Access the figure
        fig = p.fig 
        # Add a title to the Figure
        fig.suptitle(plot_title, fontsize=title_size)

        # Didn't seem to need the solution on the next line (or the solution in 
        # https://stackoverflow.com/a/34579525/8508004) so that the legend wasn't 
        # being placed within the plot over points after I fixed script a bit 
        # and tested next day. So maybe that was an artefact, but keeping here in
        # case the issue arises later because this did cause changing of position
        # although it also produced artifacts at the time.
        # p.ax.legend(loc='center right', bbox_to_anchor=(1.25, 0.5), ncol=1) # solution 
        # from https://stackoverflow.com/a/42879824/8508004 to get legend outside
        # without needing `import matplotlib.pyplot as plt`, needed by https://stackoverflow.com/a/34579525/8508004
        # BUT BOTH SOLUTIONS CAUSING ME TO SEE TWO LEGENDS??! Maybe because I call
        # `p` twice?? I made first a different name and still two?!? Or maybe 
        # because mine is so large and I am using a wrong parameter.
        # See block of code above built for `y="strand", hue="class", dodge=True,` 
        # where did need to adjust legend.


        # Add labels or hovers to the individual points indicating sites?
        # hovers will only work within Jupyter so make optional(?)

    


    # Provide feedback
    sys.stderr.write( "\nPlot generated.")

    #return plot and optionally df
    if return_df:
        return p, updated_sites_df
    else:
        return p
###--------------------------END OF MAIN FUNCTION----------------------------###
###--------------------------END OF MAIN FUNCTION----------------------------###















#*******************************************************************************
###------------------------'main' section of script---------------------------##

def main():
    """ Main entry point of the script """
    # placing actual main action in a 'helper'script so can call that easily 
    # with a distinguishing name in Jupyter notebooks, where `main()` may get
    # assigned multiple times depending how many scripts imported/pasted in.
    if scatter:
        p = plot_sites_position_across_chromosome(sites_df,swarm=False)
    else:
        p = plot_sites_position_across_chromosome(sites_df)

    
    #---------------------------------------------------------------------------
    # Saving and feedback...
    # Store plot image (and svg?) because this will most likely be the access
    # point when running on command line or equivalent (i.e., using `%run` in
    # a Jupyter notebook).
    #---------------------------------------------------------------------------

    
    if not scatter:
        # For swarmplot, use:
        #p.figure.savefig(plot_save_as_name) # from 
        # https://stackoverflow.com/a/41452422/8508004

        p.figure.savefig(plot_save_as_name, dpi = (130)) # see 
        # https://stackoverflow.com/a/41452422/8508004 for command; addition of 
        # dpisetting based on solution worked out with matplotlib earlier from 
        # http://scipy-cookbook.readthedocs.io/items/Matplotlib_AdjustingImageSize.html
        # , see `plot_panel_bar_plots_with_fit.py` for another example use.

        p.figure.savefig(plot_save_as_name[:-4]+".svg", orientation='landscape') # FOR 
        # VECTOR GRAPHICS; useful if merging into Adobe Illustrator. Based on 
        #https://neuroscience.telenczuk.pl/?p=331;I think ReportLab outputs SVG too?


    else:
        # For seaborn.lmplot, use:
        # However, also noting here that although I am only using the `p.savefig`
        # command below when using `%run` to run the the script from the 
        # command line in a notebook cell. I am also seeing the plot in the notebook
        # output below the cell in a freshly-started instance and I thought I was 
        # not seeing it yesterday.
        
        # Note that the, first time in fresh Jupyter instance on Azure with 
        # default, i.e., without `, dpi = (1600)` in `.savefig` command, default 
        # plot size appeared lrger.)
        #p.savefig(plot_save_as_name) # see https://stackoverflow.com/a/39482402/8508004
        # Note that the, first time in freshly started Jupyter instance on Azure 
        # with  default, i.e., without `, dpi = (1600)` in `.savefig` command,  
        # default plot size appeared larger when I ran script with `%run`,.i.e., 
        # calling from command line. But then shrunk when I re-ran script the 
        # same way the second time.) 

        p.savefig(plot_save_as_name, dpi = (130)) # see 
        # https://stackoverflow.com/a/39482402/8508004 for command; addition of dpi
        # setting based on solution worked out with matplotlib earlier from 
        # http://scipy-cookbook.readthedocs.io/items/Matplotlib_AdjustingImageSize.html
        # , see `plot_panel_bar_plots_with_fit.py` for another example use.

        # Also noting here that although I am only using the `p.savefig`
        # command below when using `%run` to run the the script from the 
        # command line in a notebook cell. I am also seeing the plot in the notebook
        # output below the cell in a freshly-started instance and I thought I was 
        # not seeing it yesterday.

        p.savefig(plot_save_as_name[:-4]+".svg", orientation='landscape') # FOR 
        # VECTOR GRAPHICS; useful if merging into Adobe Illustrator. Based on 
        #https://neuroscience.telenczuk.pl/?p=331;I think ReportLab outputs SVG too?


    #sns.plt.show()
    # `sns.plt.show(g)` from https://stackoverflow.com/a/42879824/8508004 might 
    # be useful when using main function in a cell in a notebook. Noting that
    # now here.

    
    sys.stderr.write("\nPlot image saved to: {}\n".format(plot_save_as_name))
    sys.stderr.write(
            "Plot image saved to: {}\n".format(plot_save_as_name[:-4]+".svg"))


    #------------------END of 'Saving and feedback' section---------------------
    #---------------------------------------------------------------------------





        







if __name__ == "__main__" and '__file__' in globals():
    """ This is executed when run from the command line """
    # Code with just `if __name__ == "__main__":` alone will be run if pasted
    # into a notebook. The addition of ` and '__file__' in globals()` is based
    # on https://stackoverflow.com/a/22923872/8508004
    # See also https://stackoverflow.com/a/22424821/8508004 for an option to 
    # provide arguments when prototyping a full script in the notebook.
    ###-----------------for parsing command line arguments-----------------------###
    import argparse
    parser = argparse.ArgumentParser(prog='plot_sites_position_across_chromosome.py',
        description="plot_sites_position_across_chromosome.py Takes a \
        dataframe of 'sites', specified with starts and ends, and plots the \
        position/distribution of \
        these 'sites' along the chromosome. Mainly made to highlight where \
        different classes of these sites/elements occur relative the other \
        classes other across a chromsome. 'sites' \
        is meant as a generalized term and they way they are specified allows \
        using multiple classes or types in the same plot. The different types \
        or classes are plotted in different color.\
        **** Script by Wayne Decatur   \
        (fomightez @ github) ***")

    parser.add_argument("sites_data", help="sites dataframe file (pickled).\
        Start and end of the sites need to specified. The identifiers for the \
        sites should begin with a alphanumeric string of characters preceding \
        a dash. After the dash should be an integer with a unique identifing \
        whole number (integer). Valid examples: `V-3` or `M1-18`, without the.\
        quotes/ticks.\
        ", metavar="SITES_DATA")
    parser.add_argument('-o', '--output', action='store', type=str, 
    default= plot_save_as_name, help="OPTIONAL: Set file name for saving plot\
    . File extension optional/moot. \
    If none provided, '{}' will be used.".format(plot_save_as_name)) 
    parser.add_argument("-sc", "--scatter",help=
    "add this flag to plot sites as points where y-value is based on arbitrary\
    membership value in site class, instead of a swarm plot (Not Recommended).",
    action="store_true")




    #I would also like trigger help to display if no arguments provided because 
    # need at least one for url
    if len(sys.argv)==1:    #from http://stackoverflow.com/questions/4042452/display-help-message-with-python-argparse-when-script-is-called-without-any-argu
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()
    sites_df= args.sites_data
    scatter = args.scatter
    if args.output != plot_save_as_name:
        #print ("output argument detected") # For Debugging
        # because I want to make sure file extension matches what will be 
        # needed later, delete anything user put and put what I need.
        if '.' in args.output:
            main_part_of_name, file_extension = os.path.splitext(args.output) #from 
            # http://stackoverflow.com/questions/541390/extracting-extension-from-filename-in-python
            plot_save_as_name = main_part_of_name + ".png"
        else:
            plot_save_as_name = args.output + ".png"

    main()

#*******************************************************************************
###-***********************END MAIN PORTION OF SCRIPT***********************-###
#*******************************************************************************
