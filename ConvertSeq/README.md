Sequence Converters
===================

**My collection of files to deal with converting sequences and sequence-related files.**

---

**Description of each script**

- convert_fasta_to_reverse_complement.py
> DNA -->  reverse complement DNA sequence

`convert_fasta_to_reverse_complement.py` takes a sequence file in FASTA format and converts the sequence (or sequences) in it to the reverse complement.

- ConvertFASTAdnaSEQtoRNA.py     

> DNA --> mRNA  
`ConvertFASTAdnaSEQtoRNA.py` takes as input a file of nucleotide sequences in FASTA format from NCBI and changes all 'T's to 'U's in the sequence, effectively changing it to mRNA sequence. Importantly, it doesn't alter the [description line](http://blast.ncbi.nlm.nih.gov/blastcgihelp.shtml) of the FASTA entries, i.e., the line beginning '>' before the sequence. Note there is no check if these are protein or nucleic sequences. Thus it will change T to U in protein if you try to give it protein sequences.

**example of input and output for `ConvertFASTAdnaSEQtoRNA.py`:**

original input:
```
>gi|392920960|ref|NM_073597.2| Caenorhabditis elegans Protein FLR-2 (flr-2) mRNA, complete cds
ATGGGCTCCAAAGCACGAGCACGACGACGTTTAAGTTGTTTTTTAAGCGTTTTTGTTGTGACATGCTTAT
TACAGTACTGCACAGCAGGTGTTACTAAGAATAATAGTTGCAAAAAAGTTGGAGTGGAGGAACTTATAGA
TGAAGAAGGCTGTGATTTGATGATAATTCGAATCAATCGATGCAGTGGGCATTGCTTCTCATTTACATTT
CCTAATCCCTTAACGAAAAAATATTCAGTGCATGCGAAGTGCTGCCGGATGGTTGAATGGGAAATGCTTG
AAACAGAATTAAAATGTTCCAAAGGAAACCGAAATCTTCGAATACCATCTGCAACACAATGTGAATGTTT
TGATTGTCTTGTTCGATAG

>gi|XM_003122551.1:0-390(+) PREDICTED: glycoprotein hormone alpha-2-like [Sus scrofa] mRNA cds , corresponds to cds for protein sequence GI_number 311247365
ATGCCCATGGCCTCCCCCCAAACCCTGCTCCTCTGCCTGCTGGTCCTGGCAATCCCTGAA
GGCCAGGGTCAGCAGGCAGCCATCCCAGGCTGCCACTTGCACCCCTTCAACGTGACCGTG
CGAAGTGACCGCCAAGGCACCTGCCAGGGCTCCCATGTGGCACAGGCCTGTGTGGGCCAC
TGTGAGTCCAGTGCCTTCCCATCCCGGTACTCCGTGCTGGTGGCCAGCGGCTATCGACAC
AACATCACCTCCGTCTCTCAGTGCTGCACCATCAGCAGCCTGAGGAAGGTGAAGGTGCAG
CTGCACTGTGGGGGGGACCGGAGGGAGGAGCTGGAGATCTTCACGGCCAGGGCCTGCCAG
TGCGACATGTGTCGCCTCTCACGCTACTAG
```
example command to run:
`python ConvertFASTAdnaSEQtoRNA.py fasta_seq.txt`

stderr will report information as script runs.

actual output to be found in new output file after running:
```
>gi|392920960|ref|NM_073597.2| Caenorhabditis elegans Protein FLR-2 (flr-2) mRNA, complete cds
AUGGGCUCCAAAGCACGAGCACGACGACGUUUAAGUUGUUUUUUAAGCGUUUUUGUUGUGACAUGCUUAU
UACAGUACUGCACAGCAGGUGUUACUAAGAAUAAUAGUUGCAAAAAAGUUGGAGUGGAGGAACUUAUAGA
UGAAGAAGGCUGUGAUUUGAUGAUAAUUCGAAUCAAUCGAUGCAGUGGGCAUUGCUUCUCAUUUACAUUU
CCUAAUCCCUUAACGAAAAAAUAUUCAGUGCAUGCGAAGUGCUGCCGGAUGGUUGAAUGGGAAAUGCUUG
AAACAGAAUUAAAAUGUUCCAAAGGAAACCGAAAUCUUCGAAUACCAUCUGCAACACAAUGUGAAUGUUU
UGAUUGUCUUGUUCGAUAG

>gi|XM_003122551.1:0-390(+) PREDICTED: glycoprotein hormone alpha-2-like [Sus scrofa] mRNA cds , corresponds to cds for protein sequence GI_number 311247365
AUGCCCAUGGCCUCCCCCCAAACCCUGCUCCUCUGCCUGCUGGUCCUGGCAAUCCCUGAA
GGCCAGGGUCAGCAGGCAGCCAUCCCAGGCUGCCACUUGCACCCCUUCAACGUGACCGUG
CGAAGUGACCGCCAAGGCACCUGCCAGGGCUCCCAUGUGGCACAGGCCUGUGUGGGCCAC
UGUGAGUCCAGUGCCUUCCCAUCCCGGUACUCCGUGCUGGUGGCCAGCGGCUAUCGACAC
AACAUCACCUCCGUCUCUCAGUGCUGCACCAUCAGCAGCCUGAGGAAGGUGAAGGUGCAG
CUGCACUGUGGGGGGGACCGGAGGGAGGAGCUGGAGAUCUUCACGGCCAGGGCCUGCCAG
UGCGACAUGUGUCGCCUCUCACGCUACUAG


```
Using the example command above, the output produced would be in the file `fasta_seq_mRNAconv.txt`.


- collapse_large_unknown_blocks_in_DNA_sequence.py
> DNA with a large block or blocks of uknown nucleotide repeats -->  DNA with small 'blocks' of uknown nucleotide repeats of a set size 

`collapse_large_unknown_blocks_in_DNA_sequence.py` takes a DNA sequence in FASTA format and reduces the large blocks of unknown nucleotides, represented by repeats of uninterrupted `N`s, to be shorter. This is meant to be used for sequences that appear to have arbitrarily-sized blocks of unknown nucleotides so the blocks are short and easier to align / or visualize how the should align in multiple sequence alignment in cases where it looks like the block may not be as large as arbitarily made. (Or to help assess if that may be the case.)

The output sequence remains in FASTA format. Also works if the provided FASTA file is a multi-FASTA, i.e., contains multiple sequences in FASTA format in the one file. All sequences in the file will have the large blocks of repeated Ns collapsed to a small size. 

Typcially the output from this script would then be submitted to software that performs a MULTIPLE SEQUENCE ALIGNMENT to see if some conerved elements are better aligned now, indicating the perhaps arbitrary blocks were making it hard to see the adjacent placement while keeping in mind these new small blocks of `N`s are indeed 'fuzzy'.  
More explanation on that:  
In order to better judge some sequences extracted from alignments against a chromosome/genome, I've found that reducing apparently arbitrarily sized repeats of unknown nucleotides, represented as uninterrupted repeats of `N`s, such as "NNNNNNNNNNNNNNN" for example, can help in assessing the sequences extracted from a much larger alignment or help prepare them for aligning again, especially when additional knowledge, such as matches in flanking sequence seem to suggest the spacing doesn't match the number of unknown nucleotides shown in the represtation. In other words, often in assemblies generated from Illumina pipelines 'educated guesses' or abtirary numbers (such as 50) `N`s in a row will be introduced for smallish 'gaps' (i.e., ones where software has concluded an intact scaffold is still indicated) in the sequence assembly and reducing these to a smaller size can sometimes make alignments of sub-elements in different rows more obvious.  This script, `collapse_large_unknown_blocks_in_DNA_sequence.py,` does this reduction in size of the blocks. It is important to be clear this was done and keep both sets of data.





Related scripts
---------------

Since sequence manipulations are at the heart of many of my computational endeavors, several other of my code repositories hold code that is related. Here are some:


 - In the `RetrieveSeq` folder is a script that takes a list of protein sequence records in FASTA format and gets the mRNA sequence correspoding to each. See `GetmRNAforProtein.py`.

- In the `RetrieveSeq` folder is a script that takes a list of protein sequence records in FASTA format and gets the mRNA sequence or at least the coding sequence corresponding to each one. See `GetmRNAorCDSforProtein.py`.

- In the `AdjustFASTA_or_FASTQ` folder is a script that takes a list of sequence records in FASTA format and adds the source organism information to the description line for each entry. See `add_source_organism_info_to_FASTA.py`.

- In the `CompareFASTA_or_FASTQ` folder is a script that takes two lists of sequence records in FASTA format and checks if the source organisms overlap. See `compare_organisms_in_two_files_of_fasta_entries.py`.

- In the [`sequencework/alignment-utilities`](https://github.com/fomightez/sequencework/tree/master/alignment-utilities) folder is a script that takes an alignment and generates that alignment with the reverse complement sequences. It is called `reverse_complement_of_clustal_alignment.py`. One of the core actions it does is use Biopython's sequence object with the `.reverse_complement()` method to make the reverse complementats.

- In the `Adjust_Annotation` folder are some scripts for dealing with converting sequence-associated annotation files to be more useful.

- In the `SOLiD` folder are some scripts for dealing with converting Applied Biosystems SOLiD data from NCBI's Short Read Archive into a useful form.

- In the `UGENE_help` folder in a different repository is some guidance and one script for converting between Unipro's UGENE sequence anysis software and Serial Cloner. See [here](https://github.com/fomightez/UGENE_help) for `UGENE help`.

- In the `yeastmine` folder in a different repository are some scripts, such as `geneID_list_to_systematic_names.py`, useful for converting identifiers of gene and protein identifiers to other useful forms. See [here](https://github.com/fomightez/yeastmine) for `yeastmine` work.


Related scripts by others
------------------------

I tend to use Biopython's sequence object with the `.reverse_complement()` method to generate reverse complements of sequences (see above) but these posts are good to know about:

- [python code for getting the reverse complement DNA strand](http://crazyhottommy.blogspot.com/2013/10/python-code-for-getting-reverse.html)

- [ranking code for making reverse complements - "What is the fastest way to get the reverse complement of a DNA sequence in python?"](https://bioinformatics.stackexchange.com/questions/3583/what-is-the-fastest-way-to-get-the-reverse-complement-of-a-dna-sequence-in-pytho?rq=1), leads to an nice repo with the data [here](https://github.com/conchoecia/fastest_rc_python)
