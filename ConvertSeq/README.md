Sequence Converters
===================

**My collection of files to deal with converting sequences and sequence-related files.**

---

**Description of each script**

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

- In the `UGENE_help` folder in a different repository is some guidance and one script for converting between Unipro's UGENE sequence anysis software and Serial CLoner. See [here](https://github.com/fomightez/UGENE_help) for `UGENE help`.

- In the `yeastmine` folder in a different repository are some scripts, such as `geneID_list_to_systematic_names.py`, useful for converting identifiers of gene and protein identifiers to other useful forms. See [here](https://github.com/fomightez/yeastmine) for `yeastmine` work.


Related scripts by others
------------------------

I tend to use Biopython's sequence object with the `.reverse_complement()` method to generate reverse complements of sequences (see above) but this post is good to know about:

[python code for getting the reverse complement DNA strand](http://crazyhottommy.blogspot.com/2013/10/python-code-for-getting-reverse.html)
