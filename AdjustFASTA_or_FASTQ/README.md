Sequence Adjusters
===================

**My collection of files to deal with minor adjustments to sequences and sequence-related files.**

---

**Description of each script**

- add_source_organism_info_to_FASTA.py
> FASTA file -->  FASTA file with organism name info injected into description line / header


Purpose:  
Takes a sequence file in FASTA format anddds source organism information in brackets to end of descrption line of each
FASTA-formatted record for protein or nucleic acid in a file. The source
organism information added will be the genus and species plus any other
information such as strain and variety associated with that sequence entry.

**The sequences in the provided file must be all of the same kind, either all
nucleic or all protein.** The script will try to detect which one
based on the first sequence in the provided file. **Subsequently, the script
will fail if the sequences in the file are a mix of protein and nucleic.**


The script will leave the description line untouched if it already has some text
occuring between brackets as it assumes that is the source organism designation
as is typical in many protein entry description lines of FASTA-formatted
sequences at NCBI.

An typical use case is to prepare sequence files to use with my associated script, `compare_organisms_in_two_files_of_fasta_entries.py`, that is found in [my CompareFASTA_or_FASTQ sub-repo](https://github.com/fomightez/sequencework/tree/master/CompareFASTA_or_FASTQ).

There is a notebook demonstrating this script available in sessions launched from [here](https://github.com/fomightez/cl_sq_demo-binder) by clicking the `launch binder` button. After the session launches select the notebook 'Demo of script to add organism names to records in FASTA file' from the list of available notebooks. Static version of the notebooks is [nicely rendered here](https://github.com/fomightez/cl_sq_demo-binder/blob/master/notebooks/demo%20add_source_organism_info_to_FASTA.ipynb).



To GET HELP/MANUAL, enter on command line:

```python
python add_source_organism_info_to_FASTA.py  --help
```

-or- if using in a Jupyter notebook:

```python
%run add_source_organism_info_to_FASTA.py  --help
```

Won't work once GI numbers phased out as detailed [here](https://ftp.ncbi.nlm.nih.gov/genbank/gbrel.txt).

Originally written in Python 2.7; however, tested & works in 3.


&nbsp;<p></p>


- delete_seq_following_pattern_within_FASTA.py
> sequence(s) in FASTA file and a pattern to match   --->  sequence deleted after match in derived FASTA file 

There is a [demo notebook for this script in this repo](https://github.com/fomightez/sequencework/blob/master/AdjustFASTA_or_FASTQ/demo%20delete_seq_following_pattern_within_multiFASTA.ipynb) that can be viewed nicely displayed [here](https://nbviewer.org/github/fomightez/sequencework/blob/master/AdjustFASTA_or_FASTQ/demo%20delete_seq_following_pattern_within_multiFASTA.ipynb).

`delete_seq_following_pattern_within_FASTA.py` takes a sequence pattern string, a sequence file (FASTA-format), and a record id, and deletes any sequence following the sequence pattern. In other, words it trims the specified sequence, to make the first match to the pattern the new end. (The FASTA-formatted sequence file is assumed by default to be a multi-FASTA, i.e., multiple sequences in the provided file, although it definitely doesn't have to be. In case it is only a single sequence, the record id becomes moot and you can enter any nonsense for that argument. Nothing will be returned but a copy of the FASTA sequence file with the truncated sequence will be produced.

The provided sequence pattern will be matched regardless of case, as both the input sequence and pattern to search will be converted to lowercase. Beyond being insensitive of the case, REGULAR EXPRESSION SEARCH TERM SYNTAX IS ACCEPTABLE in the provided sequence pattern.

Note that if there is only one record in the specified sequence file, the record id is moot and you can instead provide any string for that parameter as it will be ignored. This makes the script more flexible in cases where sequence files aren't complex as the user doesn't need to provide an actual record id. 
This script is meant to be used after you have performed a large alignment, say of an entire chromosome, in order to have individual occurrences of related segments fall linearly with where they match up along the span of the sequence. Often due to large (seeming-to-be) arbitratrily-sized blocks of repeated unknown nucleotides (which are often good to 'collapse', see `collapse_large_unknown_blocks_in_DNA_sequence.py`) the 'ends' of regions often fail to get extracted exactly right and you can end up with some sequences that trail on for longer than they should.

It is designed to handle/filter gaps ('dashes') in the provided sequence patterns. The idea being that the known sequence ends may be manually extracted from sequence alignments. This way the user is not wasting time removing the gap indications / dashes from the collected text lines. The default handling of removing the gaps to ignore them can be overriden. The idea is that maybe you'll have a multiple sequence alignment file saved as FASTA with dashes, i.e., aligned FASTA file format and may want to use this script. 


There is a [demo notebook for this script in this repo](https://github.com/fomightez/sequencework/blob/master/AdjustFASTA_or_FASTQ/demo%20delete_seq_following_pattern_within_multiFASTA.ipynb) that can be viewed nicely displayed [here](https://nbviewer.org/github/fomightez/sequencework/blob/master/AdjustFASTA_or_FASTQ/demo%20delete_seq_following_pattern_within_FASTA.ipynb).

&nbsp;<p></p>

- remove_seq_from_multiFASTA_with_match_in_description.py
> sequence(s) in FASTA file and a pattern to match   --->  sequences without any that have match in description line 
remove_seq_from_multiFASTA_with_match_in_description.py takes a sequence file in FASTA format and a text string and then ...  

***SCRIPT NEEDS TO BE TESTED***  
***DOCUMENTATION NEEDS TO BE COMPLETED***

- permute_seq_within_FASTA_an_amount.py
> FASTA file -->  circular sequence in FASTA file start breakpoint permuted to position specified

Takes a sequence file (FASTA-format) & a record id (unless single 
sequence in file) of a circular sequence such as a plasmid or mitochondrial 
genome, and moves the specified distance in bps the start/end (breakppoint) of
that sequence within the specified FASTA record. In other words, the so the 
position that is the specified amount from the original 'start' becomes the 
new 'start' position.
(The FASTA-formatted sequence file in which to act is assumed by default 
to be a multi-FASTA, i.e., multiple sequences in the provided file, although it 
definitely doesn't have to be. In case it is only a single sequence, the
record id becomes moot, see below.) 
Nothing is returned when using this on the command line. A file will be saved
of the permuted sequence. In the case of importing and using the main funcion,
by default the name of the saved file is returned. Optionally, when calling 
the function, returning the name of the permuted file can be disabled.
Note that if there is only one record in the specified sequence file, the 
record id is moot and you can instead provide any string for that parameter 
as it will be ignored. This makes the script more flexible in cases where 
sequence files aren't complex as the user doesn't need to provide an actual 
record id.

A good reference for the use of the term 'breakpoint' for the arbitrary and 
artificial location of the 'start'/'end' of circular sequence is at 
https://software.broadinstitute.org/gatk/blog?id=23598 .

Typical use cases for this script include:
1) Permuting a circular sequence to match the published or 'standard' sequence 
for ease in comparison. Or vice versa. For example if you have 
experimentally-determined sequence that spans the artifical start/end 
breakpoint and you want to compare.
2) For a circular sequence, if you don't find an expected sequence, you may 
wish to permute to see if possible the sequence happens to fall across the
artificial start/end (end/start) breakpoint.




- replace_unusual_nts_within_FASTA.py
> sequence(s) in FASTA file --->  sequences with unusual nts replaced by a single character

As it stands, you have to edit the script itself to change the character used in the substitution. Search for `character_for_subbing` in the 'USER ADJUSTABLE VALUES' section of the script.

Takes a sequence file (FASTA-format) and replaces the unusual nts with a single character specified. It also summarizes the nts in the sequence(s). Assumes multi-FASTA, but single sequence entry is fine, too. When running on the command line, it will print out a summary table of counts of nucleotides and other character in each sequence and totals. When #calling the main function it will, by default, return a dataframe with this information.

Only valid for DNA sequences; script has no step checking for data type, and so you are responsible for verifying appropriate input.

There's **a Jupyter notebook demonstrating this script** accessible by [first clicking here](https://mybinder.org/v2/gh/fomightez/cl_sq_demo-binder/master?filepath=index.ipynb), and then when the session starts select 'Demo of script to replace unusual nts appearing in FASTA file' from the available notebooks listed.    
[A static version of the demo can be seen here](https://nbviewer.org/github/fomightez/cl_sq_demo-binder/blob/master/notebooks/demo%20replace_unusual_nts_within_FASTA.ipynb).

This script makes use of the main function in my script `summarize_all_nts_even_ambiguous_present_in_FASTA.py` to summarize the nucleotdies to provide a better sense of the amount of unusual vs. normal nucleotdies. If it all goes right, it handles this all behind the scenes. I'm just pointing it out mainly so anyone reading about this script will see that summarizing functionality is available separate from the replacement steps. (Also, if fetching of `summarize_all_nts_even_ambiguous_present_in_FASTA.py` fails, it will ask to place that script in the same location with `replace_unusual_nts_within_FASTA.py` , and so it may be nice to know why it is asking that.)

Related code and suggestions by others [here](https://www.biostars.org/p/9479894/#9479894) and [here](https://www.biostars.org/p/9538916/#9538947).


&nbsp;<p></p>


Related scripts
---------------

Since sequence manipulations are at the heart of many of my computational endeavors, several other of my code repositories hold code that is related. Here are some:


 - In the `ConvertSeq` folder is a script that converts a sequence (or sequences) to the reverse complement. See `convert_fasta_to_reverse_complement.py`.



Related scripts by others
------------------------

[SeqKit - a cross-platform and ultrafast toolkit for FASTA/Q file manipulation](https://github.com/shenwei356/seqkit) by shenwei356 has categories 'edit' and 'ordering' that seem to allow adjusting FASTA or FASTQ files. See [here](https://github.com/shenwei356/seqkit#subcommands) for subcommands listing.

[Python code to remove sequences with non-amino acid residues or ambiguous residues from a multi-FASTA file, a.k.a. remove imperfect protein sequences](https://www.biostars.org/p/9531892/#9531906) and the related, [How to identify DNA sequences with ambiguous nucleotides such as N, Y, R, W.. in a multifasta file and then remove these sequences with Biopython](https://www.biostars.org/p/486742/).
