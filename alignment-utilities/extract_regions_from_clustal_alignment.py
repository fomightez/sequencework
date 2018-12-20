#!/usr/bin/env python
# extract_regions_from_clustal_alignment.py
__author__ = "Wayne Decatur" #fomightez on GitHub
__license__ = "MIT"
__version__ = "0.1.0"


# extract_regions_from_clustal_alignment.py by 
# Wayne Decatur
# ver 0.1
#
#*******************************************************************************
# Verified compatible with both Python 2.7 and Python 3.6; written initially in 
# Python 3. 
#
#
# PURPOSE: Takes a text document of an alignment in CLUSTAL format, presumably
# a large region, such as an entire chromosome and extracts the sequences that
# corresponds to a specific region according to the ungapped numbering of the 
# positions (coordinates) for the sequence on the top line, presumably the 
# reference. (The alignment can also be provided as a Python string to the 
# script or core function when running in a computational environment, such as 
# Jupyter notebook. MAYBE THIS DOESN't WORK YET!!!!)
# 
# This is useful after you have performed a large alignment, say of an entire 
# chromosome, in order to have individual occurences of related segments fall 
# linearly with where they match up along the span of the sequence aligned.
# Also collects the raw (ungapped) sequences  and notes what range they cover
# in each FASTA entry. Useful if you need to reprocess that block of the 
# sequences for some reason later. 
#
# The actual coordinates of the extracted alignment segment are in the 
# description line of the extracted sequence FASTA entries, as well as a text
# (tab-delimited) file as a table of the start and end of actual segments 
# represented for each sequence.
#
# It is designed to handle gaps indicated with either periods or dashes.
# For those on the reverse strand, there is a related script that takes the
# alignment block and makes the reverse complement of that, called 
# `reverse_complement_of_clustal_alignment.py`.
# 
# I'll mention it here in case there are concerns from folks that know that 
# python/bopython are zero-indexed and wondering what numbering system might be
# used in the resulting data. Every effort has been made to insure that 
# numberings returned should match those provided; nothing should be shifted by 
# one relative the input. In other words, if the user provides data using the
# typcial, conventional numbering where the first residue is numbered one, the 
# generated coordinate data will also be in the 'common language' form where the
# residue #1 corresponds to the first in the sequence.
#
#
# 
#
# Written to run from command line or imported into/pasted/loaded inside a 
# Jupyter notebook cell. When doing in Jupyter (IPython, I believe) you can skip
# the file save intermediate, see https://git.io/vh8Mi for similar advanced 
# examples.
#
#
#
#
#
#
#
#
#
#
# If you are Wayne, see `Collecting yeast XXXXX XXXXXX (XXXXXX) May 2018.md` for 
# impetus behind this script.
#
#
#
# Dependencies beyond the mostly standard libraries/modules:
# Biopython
#
#
# VERSION HISTORY:
# v.0.1. basic working version. But VERY SLOW IF THE SEQUENCES ALIGNED OVER 10K.
# Kept file of that version though as parts might be useful when looking to 
# avoid biopython in order to limit dependencies. (Cannot imagine when I might 
# need that though given convenience of MyBinder.org.)
# v.0.2. converted handling alignments and extraction part to Biopython. Last 
# part already made use of biopython in version 0.1.

#
# To do:
# - why do I have `name_basis` in main function and even a special example call
# to function. I don't see where it gets used here?? Just not implemented yet, 
# but needed?
#
#
#
#
# TO RUN:
# Examples,
# Enter on the command line of your terminal, the line
#-----------------------------------
# python extract_regions_from_clustal_alignment.py ALIGNMENT_TEXT_FILE 101-200
#-----------------------------------
#
# Issue `extract_regions_from_clustal_alignment.py -h` for details.
# 
#
#
# To use this after importing/pasting or loading into a cell in a Jupyter 
# notebook, specify the alignment file (or results as a string) and the 
# range in the call to the main function similar to below:
# extract_regions_from_clustal_alignment("test.clustal",region_str="101-200")
# 
#-or-
# To specify OPTIONAL 'file name'-like string to base names of generated files 
# on, specify `name_basis` like on next line when calling the function:
# extract_regions_from_clustal_alignment("test.out",region_str="101-200",name_basis="rt_arm_chrV_sensu.clustal")
# **`name_basis` ignored when a real file is supplied; presently, 
# "alignment.clustal" is default.** 
#
# 
#
'''
CURRENT ACTUAL CODE FOR RUNNING/TESTING IN A NOTEBOOK WHEN IMPORTED/LOADED OR 
PASTED IN ANOTHER CELL:
extract_regions_from_clustal_alignment("test.clustal",region_str="101-200")
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


suffix_for_saving = "_extracted" #to be used for naming the output automatically
# when running script from command line to act on an input file

coordinates_delimiter_default = "-" #change to ":" to use a colon to specify the 
# positions range to span. Mainly meant for advanced/power users because for
# the command line you can just use the `--use_colon` (or `-uc`) flag. And if
# using Jupyter cell you can specify `use_colon = True` when calling the main 
# function.



#
#*******************************************************************************
#**********************END USER ADJUSTABLE VARIABLES****************************













#*******************************************************************************
#*******************************************************************************
###DO NOT EDIT BELOW HERE - ENTER VALUES ABOVE###

import sys
import os
from collections import Counter
from Bio import SeqIO
from Bio import AlignIO
from Bio.Alphabet import generic_dna 
from Bio.Seq import Seq 
from Bio.SeqRecord import SeqRecord 
from Bio.Align import MultipleSeqAlignment 



###---------------------------HELPER FUNCTIONS---------------------------------###


def column_from_residue_number(aln, id, res_no):
    '''
    Takes an alignment an identifier and a residue number and returns the column
    that the specified residue number occurs in.

    from https://www.biostars.org/p/93805/#93811
    '''
    rec = next((r for r in aln if r.id == id), None)
    j = 0
    for i, res in enumerate(rec.seq):
        if res!='-':
            if j==res_no:
                return i
            j+=1

def get_aln_index_and_real_pos(sequence):
    '''
    a generator that takes a sequence and then returns next position in 
    alignment (equivalent to index, i.e., zero-indexed) and the actual position 
    in the contiguous sequence to which that corresponds. So the second-value
    is the ungapped position in common terms.
    For example, if sequence is `a-b`,the first values returned are `0,1`.
    The second values returned are `1,1`. And third returned values are `2,2`.
    '''
    indx = 0
    while True:
        yield indx, len(sequence[:indx+1].replace("-","")) #because second 
        # value after colon means up to but not including, the `+1` there allows
        # getting first character when index is set to zero
        indx+=1

def generate_output_file_name(file_name,suffix_for_saving, fa=False, tb = False):
    '''
    Takes a file name as an argument and returns string for the name of the
    output file. The generated name is based on the original file
    name.

    Added `fa` parameter so can say if FASTA and then make file extension `.fa`.
    Added `tb` parameter so can say if tab-delimited table and then make file 
    extension `.tsv`.

    Specific example
    =================
    Calling function with
        ("alignment.clustal")
    returns
        "alignment_ADJUSTED.clustal"
    '''
    main_part_of_name, file_extension = os.path.splitext(
        file_name) #from 
    #http://stackoverflow.com/questions/541390/extracting-extension-from-filename-in-python
    if '.' in file_name:  #I don't know if this is needed with the os.path.splitext method but I had it before so left it
        if fa:
            return main_part_of_name + suffix_for_saving  + ".fa"
        elif tb:
            return main_part_of_name + suffix_for_saving  + ".tsv"
        else:
            return main_part_of_name + suffix_for_saving  + file_extension
    else:
        if fa:
            return file_name + suffix_for_saving + ".fa"
        elif tb:
            return main_part_of_name + suffix_for_saving  + ".tsv"
        else:
            return file_name + suffix_for_saving + ".clustal"


def any_first_words_occur_four_times(repeated_words_list):
    '''
    Return True if any words in the list occur four times
    '''
    most_common,num_most_common = Counter(
        repeated_words_list).most_common(1)[0] # based on
        # https://stackoverflow.com/a/6987358/8508004
    return num_most_common >= 4

###--------------------------END OF HELPER FUNCTIONS---------------------------###
###--------------------------END OF HELPER FUNCTIONS---------------------------###



#*******************************************************************************
###------------------------'main' function of script---------------------------##

def extract_regions_from_clustal_alignment(
    alignment, region_str, use_colon = False,
    suffix_for_saving = suffix_for_saving, name_basis="alignment.clustal"):
    '''
    Main function of script. 
    It will take an alignment text file in CLUSTAL format, presumably a large 
    region, such as an entire chromosome and extracts the sequences that 
    corresponds to a specific region according to the ungapped numbering of the 
    positions (coordinates) for the sequence on the top line, presumably the 
    reference. It saves these in aligned and unaligned format. Also reports the
    coordinates from where they were extracted.

    The option to provide the results as a string is to handle where sending the
    data directly from shell script to Python without a typical file 
    intermediate, see the advanced notebook at https://git.io/vpr7i for 
    examples. The obvious use case for that is when working in the Jupyter 
    # notebook environment.
    '''
    # Bring in the necessary alignment data:
    #---------------------------------------------------------------------------

    try:
        with open(alignment, 'r') as the_file:
            file_name = alignment
            #alignment = the_file.read() #<--from when not using Bopiython to 
            # handle alignment file and extraction
            alignment = AlignIO.read(the_file, "clustal")
    except (TypeError,OSError,IOError) as e:
        file_name = name_basis
        pass # pass because instead just want to use results as a string because
        # probably provided as a string directly piped from PatMatch to Python
        # without file intermediate. The results variable will already be a
        # string in that case and so ready to go and don't need to read file.
        # The `except` above the pass has three errors it excepts because the
        # first was a TypeError seen when I tried to pass a file-like string to
        # `with open` because it seems `with open` method is incompatible with
        # use of StringIO, I think, which I usually use to try to pass things 
        # associated with file methods string. (I qualifed it with 'I think' b/c 
        # questions on stackoverlflow seemed to agree but I didn't try every 
        # possibility because realized this would probably be a better way to 
        # handle anyway.) That TypeError except got me to the next issue which
        # was trying the string as a file name and getting it was too long, and 
        # so I added the `OSError` catch and that seemed to make passing a 
        # string into the function work. `IOError` seemed to handle that same
        # thing in Python 2.7.
        # Note "FileNotFoundError is a subclass of OSError"
        # (https://stackoverflow.com/a/28633573/8508004)

    # feedback
    sys.stderr.write("Alignment file read...")



    '''
    # Go through multiple sequence alignment file parsing it to collect the 
    # identifiers. Also, identify the first identifier in each alignment block.
    #---------------------------------------------------------------------------
    first_words = []
    first_id_needed = True
    for line in alignment.split("\n"):
        if line and not line.isspace():
            # originally had just `if line:` on above line but I was seeing
            # lines that were all spaces for the length of the line made by 
            # Clustal Omega at the bottom of a sequence blocks. And these 
            #  space-filled' were causing out of index errors with the next 
            # line, `first_word = line.split()[0]` because nothing in first 
            # position if default `split` run on such a line. These 
            # 'space-filled' lines were not encountered when I developed with 
            # Kalign-generated MSAs.
            first_word = line.split()[0]
            if first_word in first_words:
                if first_id_needed:
                    first_identifier = first_word
                    first_id_needed = False
                first_words.append(first_word)
                if any_first_words_occur_four_times(first_words):
                    # want to collect those that occur at leas three times 
                    # because if an identifier has occured four times, than 
                    # others should number three and don't want anything 
                    # occurring less since usually there is a header line in
                    # Clustal alignments describing source (and/or version)
                    the_count = Counter(first_words)
                    aln_ids = [k for k, v in the_count.items() if v > 2] # based
                    # on https://stackoverflow.com/a/26773120/8508004 and
                    # https://stackoverflow.com/a/30418498/8508004, to work in 
                    # 2.7 and 3
                    break # because have gone far enough to collect the identifiers
            else:
                first_words.append(first_word)
    # feedback
    sys.stderr.write(
        "top line identifier determined as '{}'...".format(first_identifier))
    # NOTE WILL GET `UnboundLocalError: local variable 'first_identifier' 
    # referenced before assignment` if no alignment/file present!!!!




    # Go through multiple sequence alignment file parsing it as needed to 
    # accumulate full alignment for each sequence identifier code.
    # Since collected identifiers above, can use them.
    #---------------------------------------------------------------------------
    alignment_dict = {}
    for line in alignment.split("\n"):
        # determing if is one of the alignment lines
        for id_ in aln_ids:
            if line.strip().startswith(id_):
                if id_ in alignment_dict:
                    alignment_dict[id_] += line[len(
                        id_):].strip().split()[0].strip()
                else:
                    alignment_dict[id_] =  line[len(
                        id_):].strip().split()[0].strip()

    # Because: "Note the website should have an option about showing gaps as 
    # periods (dots) or dashes, we’ve shown dashes above." - SOURCE: 
    # http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc75
    # and I am adjusting sequence to allow for that by changing from '.' to
    # dashes.
    # First I had check if any have periods, but maybe just faster/easier to 
    # run `.replace()` on all anyway? Rather then check and then process.
    '''
    #if any([True if '.' in v else False for k,v in alignment_dict.items()]):
    #    alignment_dict= {k: v.replace(".","-") for k,v in alignment_dict.items()}
    '''
    alignment_dict= {k: v.replace(".","-") for k,v in alignment_dict.items()}

    # feedback
    sys.stderr.write("\nindividual lines for each sequence identifier parsed...")
    # sanity check
    # Verify that alignments are equal in length
    lengths_of_alignments = [len(v) for v in alignment_dict.values()]
    assert lengths_of_alignments.count(
        lengths_of_alignments[0]) == len(lengths_of_alignments), ("The length "
        "of all parsed alignments should be the same.") # see where 
        # `conservation_line` being made in 
        # `pretty_msa_maker_from_clustal_nucleic.py ` for mechanism used in 
        # assertion test involving Ivo van der Wijk's solution from 
        # https://stackoverflow.com/a/3844948/8508004
    '''





    # Parse the region_str to get the start and end positions of the reference 
    # sequence to specify what corresponding segment to extract from each of 
    # the aligned sequences.
    #---------------------------------------------------------------------------
    if use_colon:
        coordinates_delimiter= ":"
    else:
        coordinates_delimiter = coordinates_delimiter_default
    region_str_parts = region_str.split(coordinates_delimiter)
    start, end = int(region_str_parts[0]), int(region_str_parts[1])
    # sanity checks
    assert start < end, (
    "The user-supplied 'start' ({}) must be less than "
    "'end' ({}).".format(start,end))
    ref_seq_len = len(alignment[0].seq.ungap("-"))
    assert end < ref_seq_len, (
    "The user-supplied 'end' ({}) cannot be larger than the length of the "
    "sequence on the top line, which is '{}'.".format(end,ref_seq_len))



    # Collect the aligned sequences corresponding to the specified region in the 
    # reference sequence where numbering given in terms of ungapped, contiguous
    # sequence of reference sequence.
    #---------------------------------------------------------------------------
    # first need to determine location corresponding to start and end of 
    # contiguous sequence in the alignment
    aln_start_col = column_from_residue_number(
        alignment, alignment[0].id, start-1)  #-1 because number being supplied 
    # is in terms of 'natural' sequence numbering where first residue is 
    # NUMBER 1 and need to convert it to zero indexed bioppython.
    aln_end_col = column_from_residue_number(
        alignment, alignment[0].id, end-1)  #-1 because number being supplied 
    # is in terms of 'natural' sequence numbering where first residue is 
    # NUMBER 1 and need to convert it to zero indexed bioppython.
    
    # Now use the corresponding alignment start and end points to slice the 
    # aligned sequences.
    # Also collect the actual coordinates for each sequence in the aligned block
    # extracted because should be useful (Intended to be used with 
    # `pretty_msa_maker_from_clustal_nucleic` script)
    aligned_segmnts_dict = {}
    actl_coords_in_segments = {}
    for record in alignment:
        aligned_segmnts_dict[record.id] = record.seq[aln_start_col:aln_end_col+1]
        actual_start_coordnt = len(record.seq[:aln_start_col+1].ungap("-"))
        actual_end_coordnt = len(record.seq[:aln_end_col+1].ungap("-")) 
        actl_coords_in_segments[record.id] = (
            actual_start_coordnt,actual_end_coordnt)
    assert actl_coords_in_segments[alignment[0].id] == (start,end) , (
    "actl_coords_in_segments may be getting calculated wrong because for '{}'\n"
    "they should match start and end that they user-provided, and they don't".format(
        alignment[0].id))






    # Save aligned sequence region in clustal format
    # based on http://biopython.org/DIST/docs/api/Bio.Align-pysrc.html
    # where it says "You would normally load a MSA from a file using 
    # Bio.AlignIO, but you can do this from a list of SeqRecord objects too:"
    records = []
    for id_ in aligned_segmnts_dict:
        #records.append(
        #    SeqRecord(Seq(aligned_segmnts_dict[id_], generic_dna), id_)) # based
        ## on https://www.biostars.org/p/48797/; from when was using a string 
        # in v.0.1
        records.append(
            SeqRecord(aligned_segmnts_dict[id_], id_)) # based
        # on https://www.biostars.org/p/48797/
    align = MultipleSeqAlignment(records,
                             annotations={},
                             column_annotations={})
    out_alignment_name = generate_output_file_name(file_name,suffix_for_saving)
    AlignIO.write(align, out_alignment_name, "clustal") # based on 
    # https://biopython.org/wiki/AlignIO

    # Feedback
    sys.stderr.write("\n\nExtracted alignment segment saved as "
        "'{}'.".format(out_alignment_name))

    # Make raw, ungapped fasta of each aligned sequence.
    # Seems  need `.ungap()` method seq.io includes
    records = []
    for id_ in aligned_segmnts_dict:
        id_descript = "actual_extracted_coordinates: {}-{}".format(
            *actl_coords_in_segments[id_])
        #records.append(
        #    SeqRecord(
        #    Seq(aligned_segmnts_dict[id_], generic_dna).ungap(
        #    "-"), id_, description=id_descript))# based
        ## on https://www.biostars.org/p/48797/ and `.ungap()` method, see
        ## https://github.com/biopython/biopython/issues/1511 , and `description`
        ## from what I've seen for `id` plus https://biopython.org/wiki/SeqIO; 
        # from when was using strings in v.0.1
        records.append(
            SeqRecord(
            aligned_segmnts_dict[id_].ungap(
            "-"), id_, description=id_descript))# based
        # on https://www.biostars.org/p/48797/ and `.ungap()` method, see
        # https://github.com/biopython/biopython/issues/1511 , and `description`
        # from from what I seen for `id` plus https://biopython.org/wiki/SeqIO

    # save records as one multi-fasta file
    ungapped_file_name = generate_output_file_name(
        file_name,suffix_for_saving+"_ungapped",fa=True)
    SeqIO.write(records,ungapped_file_name, "fasta");
    # Feedback
    sys.stderr.write("\n\nExtracted sequences saved in ungapped form as\n"
        "FASTA-formatted in '{}'.".format(ungapped_file_name))

    # Save a table of the actual coordinates
    import pandas as pd
    starts, ends = zip(
        *[actl_coords_in_segments[k] for k in actl_coords_in_segments]) # based 
    # on https://stackoverflow.com/a/7558990/8508004 but I couldn't figure easy 
    # way to include the key in there so separated out
    coords = {
        'id': list(actl_coords_in_segments.keys()),
        'start': starts,
        'end': ends,
        }
    actual_coords_df = pd.DataFrame(
        coords,columns = ['id','start','end']) 
    coords_file_name = generate_output_file_name(
        file_name,suffix_for_saving+"_coords",tb=True)
    actual_coords_df.to_csv(coords_file_name, sep='\t',index = False)
    # Feedback
    sys.stderr.write("\n\nExtracted sequences coordinates in matching ungapped "
        "form saved as\n"
        "a table '{}'.".format(coords_file_name))



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
    kwargs['region_str'] = region_str
    kwargs['use_colon'] = use_colon
    kwargs['suffix_for_saving'] = suffix_for_saving
    extract_regions_from_clustal_alignment(alignment,**kwargs)
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
    parser = argparse.ArgumentParser(prog=
        'extract_regions_from_clustal_alignment.py',
        description="extract_regions_from_clustal_alignment.py \
        takes a text document of an alignment in CLUSTAL format, presumably a \
        large region, such as an entire chromosome and extracts the sequences \
        that corresponds to a specific region according to the ungapped \
        numbering of the positions (coordinates) for the sequence on the top \
        line, presumably the reference.\
        **** Script by Wayne Decatur   \
        (fomightez @ github) ***")

    parser.add_argument("align_file", help="Name of file of alignmnet text \
        file.\
        ", metavar="ALIGNMENT_FILE")

    parser.add_argument("region", help="Region to cover in extracting. Provide \
        start position and end position coordinates separated by the region \
        delimiter which is '{}' by default. You can use `--use_colon` flag to \
        change to a colon, for using something like, `201:405` instead of \
        `201{}405`. (Coordinates are meant to refer to 'common' numbering \
        scheme where first residue is numbered one, etc.)\
        ".format(coordinates_delimiter_default,coordinates_delimiter_default), 
        metavar="REGION_START-REGION_END")

    parser.add_argument("-uc", "--use_colon",help=
    "Add this flag to be able to specify that you want to use a colon in for \
    specifying the region to extact the corresponding aligned sequences.",
    action="store_true")

    parser.add_argument('-os', '--output_suffix', action='store', type=str, 
    default= suffix_for_saving, help="OPTIONAL: Set a suffix for including in \
    file name of output. \
    If none provided, '{}' will be used.".format(suffix_for_saving))



    #I would also like trigger help to display if no arguments provided because 
    # need at least one for url
    if len(sys.argv)==1:    #from http://stackoverflow.com/questions/4042452/display-help-message-with-python-argparse-when-script-is-called-without-any-argu
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    alignment = args.align_file
    region_str = args.region
    use_colon = args.use_colon
    suffix_for_saving = args.output_suffix


    main()

#*******************************************************************************
###-***********************END MAIN PORTION OF SCRIPT***********************-###
#*******************************************************************************
