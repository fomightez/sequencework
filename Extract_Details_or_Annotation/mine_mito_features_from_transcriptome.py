#!/usr/bin/env python
# mine_mito_features_from_transcriptome.py
__author__ = "Wayne Decatur" #fomightez on GitHub
__license__ = "MIT"
__version__ = "0.1.0"


# mine_mito_features_from_transcriptome.py by Wayne Decatur
# ver 0.1
#
#*******************************************************************************
# USES Python 2.7 but should be convertable via 2to3, see https://docs.python.org/3.0/library/2to3.html
# -or-after tested, below lines
# Compatible with both Python 2.7 and Python 3.6 (verified); written initially 
# in Python 2.7 to hopefully be convertible to Python 3.
#
# PURPOSE: Takes a transcriptome file of entries in FASTA format and mines the 
# mitochondrial details from definition lines like below:
#--------------------------------------------------------
#>Q0065 cdna chromosome:R64-1-1:Mito:13818:21935:1 gene:Q0065 gene_biotype:protein_coding transcript_biotype:protein_coding gene_symbol:AI4 description:Endonuclease I-SceII; encoded by a mobile group I intron within the mitochondrial COX1 gene; intron is normally spliced by the BI4p maturase but AI4p can mutate to acquire the same maturase activity [Source:SGD;Acc:S000007264]
#>Q0143 cdna chromosome:R64-1-1:Mito:51277:51429:1 gene:Q0143 gene_biotype:protein_coding transcript_biotype:protein_coding description:Dubious open reading frame; unlikely to encode a functional protein, based on available experimental and comparative sequence data [Source:SGD;Acc:S000007277]
#--------------------------------------------------------
# Also needs the exact length of the mitochondrial chromosome.
#
# Returns a dataframe with the details for the mito features.
# Specifically, in this case systematic_id, start, end, midpoint, gene_symbol
# (standard name) if there is one noted
# midpoint to be used soon(?) for sorting and, much later, labeling.
#
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
# To do:
# -
#
# TO RUN:
# Example,
# Enter on the command line of your terminal, the line
#-----------------------------------
# pyhton mine_mito_features_from_transcriptome.py transcriptome_file.fa mito_df.pkl
#-----------------------------------
# Where you replace `transcriptome_file.fa` with the specific file name of your
# transcriptome. And replace `mito_df.pkl` with name of 
# file for storing results.
# Issue `mine_mito_features_from_transcriptome.py -h` for details.
#
#
#*******************************************************************************
#


#*******************************************************************************
##################################
#  USER ADJUSTABLE VALUES        #

##################################
#
output_file_name = "mito_transcripts_dataframe.pkl"

sort_on_midpoint = True # sort dataframe on gene midpoint location






#
#*******************************************************************************
#**********************END USER ADJUSTABLE VARIABLES****************************


















#*******************************************************************************
#*******************************************************************************
###DO NOT EDIT BELOW HERE - ENTER VALUES ABOVE###

import sys
import os
import pandas as pd



###---------------------------HELPER FUNCTIONS---------------------------------###
def with_string_mine_fasta_def_line(def_line, string, list_of_indices2return=[1], return_1stWord = True):
    '''
    Takes a definition line of a fasta-formatted sequence entry and splits it
    on the string provided.
    It also takes a list of integers to specify the indices of the items 
    produced by the split to return in the order of the indices provided.
    You can also specify to return just the first word in each resulting split 
    (default option) or the entire string.
    
    

    Returns a list of the items produced by splitting that correspond to the
    indices in the list provided.
    The returned list of items will match the order of the provided indices.

    For indices that don't have associated data the Python `None` type (not same
    as the string "None") will be returned in the list.

    By default it just returns the 1st word of the split. That can be changed 
    by setting `return_1stWord` to `False` when calling the script.

    Examples: 
    provided `line` is the following line:
    >Q0065 cdna chromosome:R64-1-1:Mito:13818:21935:1 gene:Q0065 gene_biotype:protein_coding transcript_biotype:protein_coding gene_symbol:AI4 description:Endonuclease I-SceII; encoded by a mobile group I intron within the mitochondrial COX1 gene; intron is normally spliced by the BI4p maturase but AI4p can mutate to acquire the same maturase activity [Source:SGD;Acc:S000007264]
    Also needs the exact length of the mitochondrial chromosome.

    Calling with:
    with_string_mine_fasta_def_line(line, ":")

    Returns 
    ['R64-1-1']

    Calling with:
    with_string_mine_fasta_def_line(line, ":",[1,5,3,4])

    Returns 
    ['R64-1-1', '1', '13818', '21935']
    '''
    # first make sure string is even in the definition line and if isn't
    # return a list of the appropriate size containing `None` type.
    if string not in def_line:
        return [None]*len(list_of_indices2return)
    split_line = def_line.split(string)
    # Want to then return the items in the split that correspond to the indices 
    # the user-provided in list_of_indices2return, and for provided indices that 
    # don't have associated data the Python `None` type will be returned in the 
    # list.
    if return_1stWord:
        return [split_line[i].split()[0] if 0 <= i < len(
            split_line) else None for i in list_of_indices2return]
    else:
        return [split_line[i].split() if 0 <= i < len(
            split_line) else None for i in list_of_indices2return]



def with_descriptor_mine_fasta_def_line(def_line, list_of_descriptors):
    '''
    Takes a definition line of a fasta-formatted sequence entry and mines it for
    the data associated with provided descriptors.

    Returns data associated with provided descriptors in a dictionary with the 
    descriptors as keys and the associated data as keys.
    For descriptors that don't occur in a line, the Python `None` type (not same
    as the string "None") will be assigned as the value.
    The descriptors in the returned dictionary will be stripped of flanking 
    spaces if they happened to have been supplied as part of the descriptor;
    but the spaces will have been used as part of the descriptor already to 
    identify the right descriptor in the course evaluating matches.

    Because this is a simple helper function, fow now it will just return the 
    word right after the descriptor indicator.


    Example: 
    provided `line` is the following line:
    >Q0065 cdna chromosome:R64-1-1:Mito:13818:21935:1 gene:Q0065 gene_biotype:protein_coding transcript_biotype:protein_coding gene_symbol:AI4 description:Endonuclease I-SceII; encoded by a mobile group I intron within the mitochondrial COX1 gene; intron is normally spliced by the BI4p maturase but AI4p can mutate to acquire the same maturase activity [Source:SGD;Acc:S000007264]
    Also needs the exact length of the mitochondrial chromosome.

    Calling with:
    mine_fasta_definition_line(line, [" gene:", "gene_symbol:"])

    Return 
    {"gene":"Q0065","gene_symbol":"AI4"}
    '''
    result = {}
    for descriptor in list_of_descriptors:
        if descriptor in def_line:
            result[descriptor.strip()] = with_string_mine_fasta_def_line(
            def_line, descriptor)[0] #it should only return one item when 
            # called this way, but will be in a list, and so want first index.
        else:
            result[descriptor.strip()] = None
    return result



def mine_mito_features(tx_file_name, pickle_df = True):
    '''
    Takes a transcriptome file of entries in FASTA format and mines the 
    mitochondrial details from definition lines like below:
    --------------------------------------------------------
    >Q0065 cdna chromosome:R64-1-1:Mito:13818:21935:1 gene:Q0065 gene_biotype:protein_coding transcript_biotype:protein_coding gene_symbol:AI4 description:Endonuclease I-SceII; encoded by a mobile group I intron within the mitochondrial COX1 gene; intron is normally spliced by the BI4p maturase but AI4p can mutate to acquire the same maturase activity [Source:SGD;Acc:S000007264]
    >Q0143 cdna chromosome:R64-1-1:Mito:51277:51429:1 gene:Q0143 gene_biotype:protein_coding transcript_biotype:protein_coding description:Dubious open reading frame; unlikely to encode a functional protein, based on available experimental and comparative sequence data [Source:SGD;Acc:S000007277]
    --------------------------------------------------------
    Optionally you can tell it not to pickle the dataframe. Set 
    `pickle_df = False` when calling the script to not save the dataframe object.

    Returns a dataframe with the details for the mito features.
    Specifically, in this case systematic_id, start, end, strand, midpoint, 
    gene_symbol (standard name), if there  a gene_symbol noted. Otherwise put 
    `None` there

    midpoint will simply be integer close to midpoint for now.
    midpoint to be used soon(?) for sorting and, much later, labeling.

    By default, a picked version of the dataframe will be saved for use in
    other places.

    
    was considering breaking into 50 bp pieces(?) at this point, which meant
    that there'd be a row for each 50 bp and I'd note the feature that 
    overlaps but there may be more than one feature, plus probably
    best to just mine data first and do that at a later step??

    Curious things noted:
    Cannot use `gene:` as a signal to mine systematic ID because I see although there are 7135 entries when I search `^>` with Regex on in sublime text (which agrees with Salmon `quant.sf` results file), in the current transcriptome in use, there are 7136 instances returned if search `gene:` and so it seems one has `gene:` descriptor twice?? And only 7134 if I use ` gene:` where there is a space first? So one lacks? Yes, Narna virus lacks ` :gene`. So can use ` :gene` .
    (Note that `gene:` won't work because of lines like,">YCL074W cdna chromosome:R64-1-1:III:2824:3750:1 gene:YCL074W gene_biotype:pseudogene transcript_biotype:pseudogene description:Pseudogene: encodes fragment of Ty Pol protein [Source:SGD;Acc:S000000579]". where `pseudogene:` occurs. Found by using regex with `gene:.*gene:`.)
    (Also note that `>` alone doesn't work to precisely count fasta entries because of occurences of `>` to indicate italics in lines like, ">YIL156W-B cdna chromosome:R64-1-1:IX:47690:47973:1 gene:YIL156W-B gene_biotype:protein_coding transcript_biotype:protein_coding description:Putative protein of unknown function; originally identified based on homology to <i>Ashbya gossypii</i> and other related yeasts [Source:SGD;Acc:S000028511]".)
    '''
    # prepare output file for saving so it will be open and ready
    # read in the text file but skip header (see http://stackoverflow.com/questions/4796764/read-file-from-line-2-or-skip-header-row)
    with open(tx_file_name, 'r') as input_handler:
        # prepare to give feeback later
        sequence_entries = 0
        #initialize places for storing collected information
        mito_transcripts = {
        'sys_gene_id':[],
        'gene_symbol':[],
        'start':[],
        'end':[],
        'strand':[],
        'midpoint':[]
        }
        for line in input_handler:
            line = line.strip() # to remove anything odd
            # when a definition line is found, check if for a mitochondrial
            # chromosome gene / feature
            # But first for sake of counting check if this is a description line
            if line.startswith('>'):
                sequence_entries += 1
                #if line.startswith('>') and ':Mito:' in line:
                if ':Mito:' in line:
                    mined = with_descriptor_mine_fasta_def_line(
                        line, [">", " gene:", "gene_symbol:"])
                    # I am curious if the identifier after `>` always agrees 
                    # with what is after ` gene:` and so going to check, for now
                    if mined['>'] != mined['gene:']:
                        warning = (
                            "\nWARNING: the indentifier '{}' that follows the `>` "
                            "does not match the indentifier '{}' that follows "
                            "` gene:`  in the description line.\n")
                        sys.stderr.write(
                            warning.format(mined['>'],mined['gene:']))
                    mito_transcripts['sys_gene_id'].append(mined['gene:'])
                    mito_transcripts['gene_symbol'].append(mined['gene_symbol:'])
                    mined_items = with_string_mine_fasta_def_line(line, ":",[3,4,5])

                    # append mined items to the target lists, based on 
                    # https://stackoverflow.com/a/23400200/8508004 because it seemed 
                    # there must be a better way to append items from a list to 
                    # multiple lists than one at a time.
                    target_lists = [
                        mito_transcripts['start'], mito_transcripts['end'], 
                        mito_transcripts['strand']
                        ]
                    for item, lizt in zip(mined_items, target_lists):
                        lizt.append(item)
                    mito_transcripts['midpoint'].append(
                        int((int(mined_items[0])+int(mined_items[1]))/2))



        # Use the mined data in mito_transcripts dict to make a dataframe
        columns = ['sys_gene_id','gene_symbol',
        'start','end','midpoint','strand']
        df = pd.DataFrame(mito_transcripts, columns = columns)
        # Completed scan of input file & made dataframe and so give feedback.
        sys.stderr.write( "The "+ str(
            sequence_entries) + " transcripts listed in '" + tx_file_name + 
            "' have been processed to make a dataframe of "+ str(len(df)) +
            " mitochondrial genes.")


        # Fix some of the gene names or symbols (standard names)
        # For additional details of mitochondrial gene names and sources 
        # (beyond the transcriptome used here which is well-described in 
        # `Compiling RNA-Seq signal info for mitochondria.ipynb`) used in 
        # various places, such as YeastMine and the Turk et al 2013 paper:
        # See section 'Extracting mito data from quantified mito data' of 
        # `analyzing mito landscape paper data as well as normal BY4741 RNA-seq with Salmon.md`
        # Also see bottom section of 
        # `Generating TPM plots with WT and mutant mito data salmon mapped total transcriptome PLUS origins.md`
        # plus `mitochondrial-encoded gene names.txt`
        # Note that I had several of these edits done in my scripts to plot mito TPM too.
        #
        # from ttp://www.yeastgenome.org/locus/S000007335/overview
        # > "2012-10-26 Nomenclature: The systematic name of this gene was formerly tT(XXX)Q2, based on its identification in the Genomic tRNA Database as a possible pseudogene with undetermined anticodon. The feature name was changed to tT(UAG)Q2 based on literature describing it as a mitochondrial threonine tRNA isoform with UAG anticodon."
        # `tT(XXX)Q2` in my transcriptome is `tT(UAG)Q2` in most modern designations 
        # Because I want to leave `tT(XXX)Q2` to match the transcriptome, I am
        # going to set the `gene_symbol`, i.e., standard name, to `tT(UAG)Q2`
        df.loc[df['sys_gene_id'] == 'tT(XXX)Q2', 'gene_symbol'] = "tT(UAG)Q2"

        # Q0130 is ATP9 in paper but standard name is OLI1 but SGD also lists 
        # conflicts with OLI1 and so I m going to compromise with both
        df.replace("OLI1","OLI1(ATP9)", inplace=True)




        # Sort dataframe, depnding on `sort_on_midpoint` setting. 
        if sort_on_midpoint:
            df = df.sort_values('midpoint', ascending=True)
            df = df.reset_index(drop=True) # res the index to reflect re-order




        # Pickle the dataframe unless `pickle_df` is False
        if pickle_df:
            out_name = make_sure_has_extension(output_file_name) # I 
            # found problems deleting files without an extension on Jupyter 
            #instances and so might as well specify if none provide
            df.to_pickle(out_name)
            # Let user know
            sys.stderr.write( "\n\nThe dataframe has been saved as a file in a "
                "manner where other Python programs can access\nthe created "
                "dataframe (pickled).\n"
                "The dataframe is stored as '{}'".format(out_name))

        return df



def generate_output_file_name(file_name):
    '''
    Takes a file name as an argument and returns string for the name of the
    output file. The generated name is based on the original file
    name.

    Specific examples
    =================
    Calling function with
        ("test.pk1")
    returns
        "test.pkl"

    Calling function with
        ("test")
    returns
        "test.pkl"
    '''
    main_part_of_name, file_extension = os.path.splitext(
        file_name) #from 
    # http://stackoverflow.com/questions/541390/extracting-extension-from-filename-in-python
    if file_extension:
        return file_name
    else:
        return main_part_of_name + ".pkl"

def make_sure_has_extension(file_name):
    '''
    function to handle generating file name because want to be sure it has an
    extension
    '''
    return generate_output_file_name(file_name) #this is my function to normally
    # handle file names and so migth as well be consistent


###--------------------------END OF HELPER FUNCTIONS---------------------------###
###--------------------------END OF HELPER FUNCTIONS---------------------------###















#*******************************************************************************
###-----------------Actual Main portion of script---------------------------###

def main():
    """ Main entry point of the script """
    if output_file_name == 'no_output':
        df = mine_mito_features(tx_file_name, pickle_df = False)
        sys.stderr.write("\n\nThe dataframe was not stored for use elsewhere because `no_output` was specified in place of the output file name.")
    else:
        df = mine_mito_features(tx_file_name)
        







if __name__ == "__main__" and '__file__' in globals():
    """ This is executed when run from the command line """
    # Code with just `if __name__ == "__main__":` alone will be run if pasted
    # into a notebook. The addition of ` and '__file__' in globals()` is based
    # on https://stackoverflow.com/a/22923872/8508004
    # See also https://stackoverflow.com/a/22424821/8508004 for an option to 
    # provide arguments when prototyping a full script in the notebook.
    ###-----------------for parsing command line arguments-----------------------###
    import argparse
    parser = argparse.ArgumentParser(prog='mine_mito_features_from_transcriptome.py',
        description="mine_mito_features_from_transcriptome.py  takes a \
        transcriptome and make a datafram of the gene information for the \
        mitochondrial chromosome.      \
        **** Script by Wayne Decatur   \
        (fomightez @ github) ***")

    parser.add_argument("tx_file_name", help="Name of the transcriptome file. \
        ", metavar="INPUT_FILE")
    parser.add_argument("output", nargs='?', help="**OPTIONAL**Name of file \
        for storing the dataframe. If none is provided, the dataframe will be \
        stored as '"+output_file_name+"'. To force nothing to be saved, enter \
        `no_output` without quotes as output file.", 
        default=output_file_name , metavar="OUTPUT_FILE")
    # Note see https://stackoverflow.com/questions/18862836/how-to-open-file-using-argparse#comment35222484_18863004
    # for why not using `argparse.FileType` approach here.
    # See
    # https://stackoverflow.com/questions/4480075/argparse-optional-positional-arguments 
    # and 
    # https://docs.python.org/2/library/argparse.html#nargs for use of `nargs='?'` 
    # to make input and output file names optional. Note that the square brackets
    # shown in the usage out signify optional according to 
    # https://stackoverflow.com/questions/4480075/argparse-optional-positional-arguments#comment40460395_4480202
    # , but because placed under positional I added clarifying text to help 
    # description.
    # IF MODIFYING THIS SCRIPT FOR USE ELSEWHERE AND DON'T NEED/WANT THE OUTPUT 
    # FILE TO BE OPTIONAL, remove `nargs` (& default?) BUT KEEP WHERE NOT
    # USING `argparse.FileType` AND USING `with open` AS CONISDERED MORE PYTHONIC.


    #I would also like trigger help to display if no arguments provided because need at least one input file
    if len(sys.argv)==1:    #from http://stackoverflow.com/questions/4042452/display-help-message-with-python-argparse-when-script-is-called-without-any-argu
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    # Note see https://stackoverflow.com/questions/18862836/how-to-open-file-using-argparse#comment35222484_18863004
    # for why not using `argparse.FileType` approach here.
    tx_file_name= args.tx_file_name
    output_file_name = args.output


    main()

#*******************************************************************************
###-***********************END MAIN PORTION OF SCRIPT***********************-###
#*******************************************************************************
