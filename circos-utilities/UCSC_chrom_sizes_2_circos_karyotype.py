#!/usr/bin/env python
# UCSC_chrom_sizes_2_circos_karyotype.py
__author__ = "Wayne Decatur" #fomightez on GitHub
__license__ = "MIT"
__version__ = "0.1.0"


# UCSC_chrom_sizes_2_circos_karyotype.py by Wayne Decatur
# ver 0.1
#
#*******************************************************************************
# Verified compatible with both Python 2.7 and Python 3.6; written initially in 
# Python 3.
#
# PURPOSE: Takes a URL for a UCSC `chrom.sizes` file and makes a `karyotype.tab` 
# file from it for use with Circos.
# Note: to determine the URL, google `YOUR_ORGANISM genome UCSC chrom.sizes`, 
# where you replace `YOUR_ORGANISM` with your organism name and then
# adapt the path you see in the best match to be something similar to 
# "http://hgdownload-test.cse.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.chrom.sizes"
# -or-
# "http://hgdownload-test.cse.ucsc.edu/goldenPath/canFam2/bigZips/canFam2.chrom.sizes"
# 
# IMPORTANTLY, this script is intended for organisms without cytogenetic bands, 
# such as dog, cow, yeast, etc..
# Acquiring the cytogenetic bands information is described at 
# http://circos.ca/tutorials/lessons/ideograms/karyotypes/ , about halfway down 
# the page where it says, "obtain the karyotype structure from...". 
# Unfortunately, it seems the output directed to by those instructions is not
# directly useful in Circos(?). Fortunately, though as described at 
# http://circos.ca/documentation/tutorials/quick_start/hello_world/ 
# ,"Circos ships with several predefined karyotype files for common sequence 
# assemblies: human, mouse, rat, and drosophila. These files are located in 
# data/karyotype within the Circos distribution."
#
# Written to run from command line or pasted/loaded inside a Jupyter notebook 
# cell. 
#
#
#
# This script based on work and musings developed in 
# `Trying to convert k75.Umap.bedGraph to bigwig file that works at SGD jbrowse.md` 
# (specifically use of chrom.sizes) and 
# `Resources in regards to plotting information on presence or absence of signal on circular chromosome circos.md` 
# (where was describing issues with getting karyotype) and 
# http://circos.ca/tutorials/course/handouts/session-4.pdf (that shows first 
# part of Saccharomyces cerevisiae karyptype on page 6).
#
# Example input from 
# http://hgdownload-test.cse.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.chrom.sizes:
'''
chrIV   1531933
chrXV   1091291
chrVII  1090940
chrXII  1078177
chrXVI  948066
chrXIII 924431
chrII   813184
chrXIV  784333
chrX    745751
chrXI   666816
chrV    576874
chrVIII 562643
chrIX   439888
chrIII  316620
chrVI   270161
chrI    230218
chrM    85779
'''

#
#Example output (tab-separated):
'''
chr -   Sc-chrIV    chrIV   0   1531933 black
chr -   Sc-chrXV    chrXV   0   1091291 black
chr -   Sc-chrVII   chrVII  0   1090940 black
chr -   Sc-chrXII   chrXII  0   1078177 black
chr -   Sc-chrXVI   chrXVI  0   948066  black
chr -   Sc-chrXIII  chrXIII 0   924431  black
chr -   Sc-chrII    chrII   0   813184  black
chr -   Sc-chrXIV   chrXIV  0   784333  black
chr -   Sc-chrX chrX    0   745751  black
chr -   Sc-chrXI    chrXI   0   666816  black
chr -   Sc-chrV chrV    0   576874  black
chr -   Sc-chrVIII  chrVIII 0   562643  black
chr -   Sc-chrIX    chrIX   0   439888  black
chr -   Sc-chrIII   chrIII  0   316620  black
chr -   Sc-chrVI    chrVI   0   270161  black
chr -   Sc-chrI chrI    0   230218  black
chr -   Sc-chrM chrM    0   85779   black
'''

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
# - probably would be nice to add automated handling of ordering by increasing 
# chromosome number. (I've used detection of roman numerals before, see 
# `plot_expression_across_chromosomes.py) Because would need to be able to 
# store and sort, probably putting the chromosomes and lengths in a dataframe 
# instead would be a good route. Then could write a function to iterrows and 
# write the output lines.
# - possible to do: automate making ones for ones with cytogenetic bands, or is
# there not enough aside from the ones included?
#
#
#
#
# TO RUN:
# Examples,
# Enter on the command line of your terminal, the line
#-----------------------------------
# python UCSC_chrom_sizes_2_circos_karyotype.py http://hgdownload-test.cse.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.chrom.sizes

#-OR-
# python UCSC_chrom_sizes_2_circos_karyotype.py http://hgdownload-test.cse.ucsc.edu/goldenPath/canFam2/bigZips/canFam2.chrom.sizes dog_karyotype.tab --species_code dog
#-----------------------------------
# Issue `python UCSC_chrom_sizes_2_circos_karyotype.py -h` for details.
# 
#
# To use this after pasting or loading into a cell in a Jupyter notebook, in
# the next cell define the URL and then call the main function similar to below:
# url = "http://hgdownload-test.cse.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.chrom.sizes"
# UCSC_chrom_sizes_2_circos_karyotype(species_code)
#
#(`species_code_hardcoded` and `output_file_name `can be assigned in a cell 
# before calling the function as well.)
#
# Note that `url` is actually not needed if you are using the yeast one because 
# that specific one is hardcoded in script as default.
# In fact due to fact I hardcoded in defaults, just `main()` will indeed work 
# for yeast.
#
# 
#
'''
CURRENT ACTUAL CODE FOR RUNNING/TESTING IN A NOTEBOOK WHEN LOADED OR PASTED IN 
ANOTHER CELL:
UCSC_chrom_sizes_2_circos_karyotype()

-OR, just-

main()

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
## default URL
url = "http://hgdownload-test.cse.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.chrom.sizes" 
output_file_name = "karyotype.tab"

species_code_hardcoded = None # replace `None` with what you want to use,
# with flanking quotes if something appropriate is not being extracted from the
# provided URL to be used as the species code.





#
#*******************************************************************************
#**********************END USER ADJUSTABLE VARIABLES****************************


















#*******************************************************************************
#*******************************************************************************
###DO NOT EDIT BELOW HERE - ENTER VALUES ABOVE###

import sys
import os




###---------------------------HELPER FUNCTIONS---------------------------------###




def make_and_save_karyotype(chromosomes_and_length, species_code):
    '''
    Takes a dictionary of chromosome identifiers and length and makes a karyotype
    file with that information.

    Result will look like this at start of output file:
    chr - Sc-chrIV chrIV 0 1531933 black
    chr - Sc-chrXV chrXV 0 1091291 black
    ...

    Function returns None.
    '''
    # prepare output file for saving so it will be open and ready
    with open(output_file_name, 'w') as output_file:
        for indx,(chrom,length) in enumerate(chromosomes_and_length.items()):
            next_line = ("chr\t-\t{species_code}-{chrom}\t{chrom}\t0"
                "\t{length}\tblack".format(
                species_code=species_code,chrom=chrom, length=length))
            if indx < (len(chromosomes_and_length)-1):
                next_line += "\n" # don't add new line character to last line
            # Send the built line to output
            output_file.write(next_line)
    sys.stderr.write( "\n\nThe karyotype file for {} chromosomes has been saved "
            "as a file named"
            " '{}'.".format(len(chromosomes_and_length),output_file_name))


def extract_species_code_fromUCSC_URL(url):
    '''
    Take something like:
    http://hgdownload-test.cse.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.chrom.sizes

    And return:
    sacCer

    Note:
    I decided to use `''.join([i for i in s if not i.isdigit()])`, where s is a
    aprovided string, to toss digits.
    '''
    species_code = url.split("goldenPath")[1].split("/")[1]
    return ''.join([i for i in species_code if not i.isdigit()]) # remove digits


def UCSC_chrom_sizes_2_circos_karyotype():
    '''
    Main function of script. Will use url to get `chrom.sizes` file from UCSC 
    and use that to make a karyotype file for use in Circos.
    Saves the file as tab-separated values with the extension `.tab`, by
    default, to be consistent with what Circos ecosystem seems to use.

    Returns None
    '''
    # Get data from URL.
    chromosomes_and_length = {}
    # Getting html originally for just Python 3, adapted from 
    # https://stackoverflow.com/a/17510727/8508004 and then updated from to 
    # handle Python 2 and 3 according to same link.
    try:
        # For Python 3.0 and later
        from urllib.request import urlopen
    except ImportError:
        # Fall back to Python 2's urllib2
        from urllib2 import urlopen
    html = urlopen(url)
    for line in html.read().splitlines():
        #chromosome, chr_len, *_ = line.strip().split()
        # that elegant unpack above is based on 
        # https://stackoverflow.com/questions/11371204/unpack-the-first-two-elements-in-list-tuple
        # , but it won't work in Python 2. From same place, one that works in 2:
        chromosome, chr_len = line.strip().split()[:2]
        chromosomes_and_length[chromosome.decode(
            encoding='UTF-8')] = chr_len.decode(encoding='UTF-8')



    # Parse the URL for a genus/species -type identifier. (If one not provided.)
    # Note part of keeping URL separate is so that I parse it to parse out from URL 
    # first part of genus-species identifier. Here in development version that is 
    # `sacCer3`, for yeast Saccharmyces cerevisiae. Parsing
    # because of advice [here](http://circos.ca/documentation/tutorials/ideograms/karyotypes/), 
    # "Even when working with only one species, prefixing the chromosome with a 
    # species code is highly recommended - this will greatly help in creating 
    # more transparent configuration and data files."
    if species_code_hardcoded:
        species_code = species_code_hardcoded
        sys.stderr.write( "\nThe following "
                "species code will be used in the ID column "
                "in the\nproduced karyotype file: '{}'.".format(species_code))
    else:
        species_code = extract_species_code_fromUCSC_URL(url)
        if species_code == "sacCer":
            species_code = "Sc" # CUSTOMIZING; I'd prefer to use this for yeast.
        sys.stderr.write( "\nBased on the provided URL, the following "
                "species code will be used in the\nID column "
                "in the karyotype file: '{}'.\n"
                "If that is not suitable, you can re-run the script and "
                "provide one when calling\nthe script using the "
                "`--species_code` flag. Alternatively, edit "
                "the produced file with find/replace.".format(species_code))
    # With the approach in that above block, I can expose `species_code` to 
    # setting for advanced use without it being required and without need to be 
    # passed into the function.




    # Now use the data to make a karyotype file as described at 
    # http://circos.ca/documentation/tutorials/ideograms/karyotypes/ and like 
    # on page 6 of http://circos.ca/tutorials/course/handouts/session-4.pdf
    make_and_save_karyotype(chromosomes_and_length, species_code)










###--------------------------END OF HELPER FUNCTIONS---------------------------###
###--------------------------END OF HELPER FUNCTIONS---------------------------###















#*******************************************************************************
###------------------------'main' secion of script---------------------------###

def main():
    """ Main entry point of the script """
    # placing actual main action in a 'helper'script so can call that easily 
    # with a distinguishing name in Jupyter notebooks, where `main()` may get
    # assigned multiple times depending how many scripts imported/pasted in.
    UCSC_chrom_sizes_2_circos_karyotype()
        







if __name__ == "__main__" and '__file__' in globals():
    """ This is executed when run from the command line """
    # Code with just `if __name__ == "__main__":` alone will be run if pasted
    # into a notebook. The addition of ` and '__file__' in globals()` is based
    # on https://stackoverflow.com/a/22923872/8508004
    # See also https://stackoverflow.com/a/22424821/8508004 for an option to 
    # provide arguments when prototyping a full script in the notebook.
    ###-----------------for parsing command line arguments-----------------------###
    import argparse
    parser = argparse.ArgumentParser(prog='UCSC_chrom_sizes_2_circos_karyotype.py',
        description="UCSC_chrom_sizes_2_circos_karyotype.py takes a URL for a \
        UCSC chrom.sizes file and makes a karyotype.tab file.      \
        **** Script by Wayne Decatur   \
        (fomightez @ github) ***")

    parser.add_argument("URL", help="URL of chrom.sizes file at UCSC. \
        ", metavar="URL")
    parser.add_argument('-sc', '--species_code', action='store', type=str, 
        default= species_code_hardcoded, help="**OPTIONAL**Identifier \
        to use in front of chromosome names. An attempt will be made to extract \
        one and that is why this is optional.")
    parser.add_argument("output", nargs='?', help="**OPTIONAL**Name of file \
        for storing the karyotype. If none is provided, the karyotype will be \
        stored as '"+output_file_name+"'.", 
        default=output_file_name , metavar="OUTPUT_FILE")

    # See
    # https://stackoverflow.com/questions/4480075/argparse-optional-positional-arguments 
    # and 
    # https://docs.python.org/2/library/argparse.html#nargs for use of `nargs='?'` 
    # to make output file name optional. Note that the square brackets
    # shown in the usage out signify optional according to 
    # https://stackoverflow.com/questions/4480075/argparse-optional-positional-arguments#comment40460395_4480202
    # , but because placed under positional I added clarifying text to help 
    # description.
    # IF MODIFYING THIS SCRIPT FOR USE ELSEWHERE AND DON'T NEED/WANT THE OUTPUT 
    # FILE TO BE OPTIONAL, remove `nargs` (& default?) BUT KEEP WHERE NOT
    # USING `argparse.FileType` AND USING `with open` AS CONISDERED MORE PYTHONIC.



    #I would also like trigger help to display if no arguments provided because 
    # need at least one for url
    if len(sys.argv)==1:    #from http://stackoverflow.com/questions/4042452/display-help-message-with-python-argparse-when-script-is-called-without-any-argu
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    url= args.URL
    output_file_name = args.output
    species_code = args.species_code


    main()

#*******************************************************************************
###-***********************END MAIN PORTION OF SCRIPT***********************-###
#*******************************************************************************
