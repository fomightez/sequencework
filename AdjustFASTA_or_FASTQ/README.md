Sequence Adjusters
===================

**My collection of files to deal with minor adjustments to sequences and sequence-related files.**

---

**Description of each script**


- delete_seq_following_pattern_within_FASTA.py
> sequence(s) in FASTA file and a pattern to match   --->  sequence deleted after match in derived FASTA file 

There is a [demo notebook for this script in this repo](https://github.com/fomightez/sequencework/blob/master/AdjustFASTA_or_FASTQ/demo%20delete_seq_following_pattern_within_FASTA.ipynb) that can be viewed nicely displayed [here](https://nbviewer.jupyter.org/github/fomightez/sequencework/blob/master/AdjustFASTA_or_FASTQ/demo%20delete_seq_following_pattern_within_FASTA.ipynb).

`delete_seq_following_pattern_within_FASTA.py` takes a sequence pattern string, a sequence file (FASTA-format), and a record id, and deletes any sequence following the sequence pattern. In other, words it trims the specified sequence, to make the first match to the pattern the new end. (The FASTA-formatted sequence file is assumed by default to be a multi-FASTA, i.e., multiple sequences in the provided file, although it definitely doesn't have to be. In case it is only a single sequence, the record id becomes moot and you can enter any nonsense for that argument. Nothing will be returned but a copy of the FASTA sequence file with the truncated sequence will be produced.

The provided sequence pattern will be matched regardless of case, as both the input sequence and pattern to search will be converted to lowercase. Beyond being insensitive of the case, REGULAR EXPRESSION SEARCH TERM SYNTAX IS ACCEPTABLE in the provided sequence pattern.

Note that if there is only one record in the specified sequence file, the record id is moot and you can instead provide any string for that parameter as it will be ignored. This makes the script more flexible in cases where sequence files aren't complex as the user doesn't need to provide an actual record id. 
This script is meant to be used after you have performed a large alignment, say of an entire chromosome, in order to have individual occurrences of related segments fall linearly with where they match up along the span of the sequence. Often due to large (seeming-to-be) arbitratrily-sized blocks of repeated unknown nucleotides (which are often good to 'collapse', see `collapse_large_unknown_blocks_in_DNA_sequence.py`) the 'ends' of regions often fail to get extracted exactly right and you can end up with some sequences that trail on for longer than they should.

It is designed to handle/filter gaps ('dashes') in the provided sequence patterns. The idea being that the known sequence ends may be manually extracted from sequence alignments. This way the user is not wasting time removing the gap indications / dashes from the collected text lines. The default handling of removing the gaps to ignore them can be overriden. The idea is that maybe you'll have a multiple sequence alignment file saved as FASTA with dashes, i.e., aligned FASTA file format and may want to use this script. 


There is a [demo notebook for this script in this repo](https://github.com/fomightez/sequencework/blob/master/AdjustFASTA_or_FASTQ/demo%20delete_seq_following_pattern_within_FASTA.ipynb) that can be viewed nicely displayed [here](https://nbviewer.jupyter.org/github/fomightez/sequencework/blob/master/AdjustFASTA_or_FASTQ/demo%20delete_seq_following_pattern_within_FASTA.ipynb).

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


- add_source_organism_info_to_FASTA.py
> FASTA file -->  FASTA file with organism info injected in text

`add_source_organism_info_to_FASTA.py` takes a sequence file in FASTA format and ...  

***DOCUMENTATION NEEDS TO BE COMPLETED***



Related scripts
---------------

Since sequence manipulations are at the heart of many of my computational endeavors, several other of my code repositories hold code that is related. Here are some:


 - In the `ConvertSeq` folder is a script that converts a sequence (or sequences) to the reverse complement. See `convert_fasta_to_reverse_complement.py`.



Related scripts by others
------------------------

[SeqKit - a cross-platform and ultrafast toolkit for FASTA/Q file manipulation](https://github.com/shenwei356/seqkit) by shenwei356 has categories 'edit' and 'ordering' that seem to allow adjusting FASTA or FASTQ files. See [here](https://github.com/shenwei356/seqkit#subcommands) for subcommands listing.
