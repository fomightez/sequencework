# Extracting sequences from a FASTA file

Repo for my computational resources dealing with extracting sequence information from FASTA-formatted sequence files.

## My scripts

* get_seq_following_seq_from_multiFASTA.py
> sequence(s) and a pattern to match   ---> first sequence after the match in specific sequence

At the time this was developed, it was intended to help examine sequences just downstream of sequences extracted from much larger alignments to check if all pertinent sequence included already.

Takes a sequence pattern string, a sequence file (FASTA-format), and  record id, and extracts a sequence of specified size following the sequence pattern. (The FASTA-formatted sequence file is assumed by default to be a multi-FASTA, i.e., multiple sequences in the provided file, although it definitely doesn't have to be. In case it is only a single sequence, the record id becomes moot, see below.

A sequence string of the specified length will be returned. Redirect the output to a file if that is what is needed. (Unless using main function as coverd in the demo notebook, see below.)

The provided sequence pattern will be matched regardless of case, as both the input sequence and pattern to search will be converted to lowercase. Beyond being insensitive of the case, **REGULAR EXPRESSION SEARCH TERM SYNTAX IS ACCEPTABLE in the provided sequence pattern**.

Note that if there is only one record in the specified sequence file, the record id is moot and you can instead provide any string for that parameter as it will be ignored. This makes the script more flexible in cases where sequence files aren't complex as the user doesn't need to provide the record id.

This script is meant to be used after you have performed a large alignment, say of an entire chromosome, in order to have individual occurrences of related segments fall linearly with where they match up along the span of the sequence. Often due to  large (seeming-to-be) arbitratrily-sized blocks of repeated unknown nucleotides (which are often good to 'collapse', see `collapse_large_unknown_blocks_in_DNA_sequence.py`) the 'ends' of region often fail to get extracted exactly right. This script is meant to help in determining how best to clean up such instances. For example, in the obtained sequence, is there an 'end' that matches up better with the pattern of known 'ends' and should be added?


It is designed to handle/filter gaps ('dashes') in the provided sequence patterns. The idea being that the known sequence ends may be manually extracted from sequence alignments. This way the user is not wasting time removing the gap indications / dashes from the collected text lines. The default handling of removing the gaps to ignore them can be overriden. The idea is that maybe you'll have a multiple sequence alignment file saved as FASTA with dashes, i.e., aligned FASTA file format and may want to use this script.  (The caveat is that number of residues to get will then be counting the gaps / dashes too.)

There is a [demo notebook for this script in this repo](UPDATE THIS URL HERE)<---UPDATE THIS URL HERE that can be viewed nicely displayed [here](UPDATE THIS URL HERE)<----UPDATE THIS URL HERE.

## Related scripts / utilities in my other repositories

To extract sequences from a Clustal-formatted alignment file, see `extract_regions_from_clustal_alignment.py` in [my Alignment utilities sub-repository](https://github.com/fomightez/sequencework/tree/master/alignment-utilities).

## Scripts / Resources by others

- >"Fastest way to extract some seqs from a BIG fasta file based on header?" - SOURCE: https://twitter.com/macmanes/status/1045659290919481344

  Good folks with lots of ideas in that response thread.
