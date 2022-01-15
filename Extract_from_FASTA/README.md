# Extracting sequences from a FASTA file

Repo for my computational resources dealing with extracting sequence from FASTA-formatted sequence files. If you need to gather some details or other information related to the sequences, you should probably look into my ['Sequencework/Extract_Details_or_Annotation' code repository](https://github.com/fomightez/sequencework/tree/master/Extract_Details_or_Annotation/).


## My scripts

* extract_subsequence_from_FASTA.py
> sequence file ---> subsequence

Quick-take: Takes a sequence file and extracts a specified subsequence specified by providing positions to span.

There is a [demo notebook for this script in this repo](https://github.com/fomightez/cl_sq_demo-binder). Launch a binder session from there and select 'Demo of script to specific subsequence from FASTA file' to run it actively.  The particular notebook can be viewed statically, nicely displayed [here](https://nbviewer.jupyter.org/github/fomightez/cl_sq_demo-binder/blob/master/notebooks/Demo%20of%20script%20to%20get%20specific%20subsequence%20from%20FASTA%20file.ipynb) (although that version does cut off the full sequence in last cells). There is a more typical use case illustrated in a workflow [here](https://nbviewer.jupyter.org/github/fomightez/blast-binder/blob/master/notebooks/Searching%20for%20coding%20sequences%20in%20genomes%20using%20BLAST%20and%20Python.ipynb). It can be launched as an active Jupyter session by clicking `launch binder` [here](https://github.com/fomightez/blast-binder) and then selecting 'Searching for coding sequences in genomes using BLAST and Python' from the index page once the session launches.

Takes a sequence file and extracts a specified subsequence specified by providing positions to span. If there is more than one sequence, it will know which one to extract those positions form using the `record identifier`, a.k.a. the text before the first space on the description line of a FASTA entry, not including the caret symbol (less-than-symbol).

Assumes multiFASTA, but it can be a FASTA formatted file with a single entry. In that case, the `record identifier` parameter is moot and nonsensical text can be provided in its place.

Typical example command:

```
python extract_subsequence_from_FASTA.py SK1.genome.fa chrIII 101-200
```

* get_seq_from_multiFASTA_with_match_in_description.py
> sequences and a text pattern to match   ---> first sequence that has a match in description

There is a [demo notebook for this script in this repo](https://github.com/fomightez/cl_sq_demo-binder). Launch a binder session from there and select 'Demo of script to get sequence from multiFASTA file when description contains matching text' to run it actively.  The particular notebook can be viewed statically, nicely displayed [here](https://nbviewer.jupyter.org/github/fomightez/cl_sq_demo-binder/blob/master/notebooks/Demo%20of%20script%20to%20get%20sequence%20from%20multiFASTA%20file%20when%20description%20contains%20matching%20text.ipynb) (although that version does caught off the full sequence in last cells).

Takes any sequences in FASTA format and gets the first sequence with a description line containing a match to provided text string. For example, if provided a multi-sequence FASTA file and a gene identifier, such as `YDL140C`, it will pull out the first sequence matching that anywhere in the description line. Defaults to ignoring case.

Typical example command:

```
python get_seq_from_multiFASTA_with_match_in_description.py DBVPG6044.mt.pep.fa cox1
```



* get_seq_following_seq_from_FASTA.py
> sequence(s) and a pattern to match   ---> first sequence after the match in specific sequence

There is a [demo notebook for this script in this repo](https://github.com/fomightez/sequencework/blob/master/Extract_from_FASTA/demo%20get_seq_following_seq_from_FASTA.ipynb) that can be viewed nicely displayed [here](https://nbviewer.jupyter.org/github/fomightez/sequencework/blob/master/Extract_from_FASTA/demo%20get_seq_following_seq_from_FASTA.ipynb).

At the time this was developed, it was intended to help examine sequences just downstream of sequences extracted from much larger alignments to check if all pertinent sequence included already.

Takes a sequence pattern string, a sequence file (FASTA-format), and  record id, and extracts a sequence of specified size following the sequence pattern. (The FASTA-formatted sequence file is assumed by default to be a multi-FASTA, i.e., multiple sequences in the provided file, although it definitely doesn't have to be. In case it is only a single sequence, the record id becomes moot, see below.

A sequence string of the specified length will be returned. Redirect the output to a file if that is what is needed. (Unless using main function as coverd in the demo notebook, see below.)

The provided sequence pattern will be matched regardless of case, as both the input sequence and pattern to search will be converted to lowercase. Beyond being insensitive of the case, **REGULAR EXPRESSION SEARCH TERM SYNTAX IS ACCEPTABLE in the provided sequence pattern**.

Note that if there is only one record in the specified sequence file, the record id is moot and you can instead provide any string for that parameter as it will be ignored. This makes the script more flexible in cases where sequence files aren't complex as the user doesn't need to provide an actual record id.

This script is meant to be used after you have performed a large alignment, say of an entire chromosome, in order to have individual occurrences of related segments fall linearly with where they match up along the span of the sequence. Often due to  large (seeming-to-be) arbitratrily-sized blocks of repeated unknown nucleotides (which are often good to 'collapse', see `collapse_large_unknown_blocks_in_DNA_sequence.py`) the 'ends' of region often fail to get extracted exactly right. This script is meant to help in determining how best to clean up such instances. For example, in the obtained sequence, is there an 'end' that matches up better with the pattern of known 'ends' and should be added?


It is designed to handle/filter gaps ('dashes') in the provided sequence patterns. The idea being that the known sequence ends may be manually extracted from sequence alignments. This way the user is not wasting time removing the gap indications / dashes from the collected text lines. The default handling of removing the gaps to ignore them can be overriden. The idea is that maybe you'll have a multiple sequence alignment file saved as FASTA with dashes, i.e., aligned FASTA file format and may want to use this script.  (The caveat is that number of residues to get will then be counting the gaps / dashes too.)

There is a [demo notebook for this script in this repo](https://github.com/fomightez/sequencework/blob/master/Extract_from_FASTA/demo%20get_seq_following_seq_from_FASTA.ipynb) that can be viewed nicely displayed [here](https://nbviewer.jupyter.org/github/fomightez/sequencework/blob/master/Extract_from_FASTA/demo%20get_seq_following_seq_from_FASTA.ipynb).



* get_specified_length_of_end_of_seq_from_FASTA.py
> sequence(s) and a number of residues to get   ---> specified number of characters at the end of the sequence for the indicated record

There is a [demo notebook for this script in this repo](https://github.com/fomightez/cl_sq_demo-binder/tree/master/notebooks/demo%20get_specified_length_of_end_of_seq_from_FASTA.ipynb) that can be viewed nicely displayed [here](https://nbviewer.jupyter.org/github/fomightez/cl_sq_demo-binder/tree/master/notebooks/demo%20get_specified_length_of_end_of_seq_from_FASTA.ipynb). It can actively by run in your browser by going to [here](https://github.com/fomightez/cl_sq_demo-binder), clicking on the `launch binder` badge, and selecting the 'Demo of script to get specified length of sequence from end of a record in FASTA file' notebook from the listing.

At the time this was developed, it was intended to extract sequences coding for a C-terminal protein domain.

Takes a sequence file (FASTA-format), & a record id, and a number (integer), and extracts a sequence of specified length from the end of the indicated sequence. The number provided is what specifies the length extracted. (The FASTA-formatted sequence file is assumed by default to be a multi-FASTA, i.e., multiple sequences in the provided file, although it definitely doesn't have to be. In case it is only a single sequence, the record id becomes moot, see below.

A sequence string of the specified length will be returned. Redirect the output to a file if that is what is needed. (Unless using main function as coverd in the demo notebook, see below.)

Note that if there is only one record in the specified sequence file, the record id is moot and you can instead provide any string for that parameter as it will be ignored. This makes the script more flexible in cases where sequence files aren't complex as the user doesn't need to provide an actual record id.


There is a [demo notebook for this script in this repo](https://github.com/fomightez/cl_sq_demo-binder/tree/master/notebooks/demo%20get_specified_length_of_end_of_seq_from_FASTA.ipynb) that can be viewed nicely displayed [here](https://nbviewer.jupyter.org/github/fomightez/cl_sq_demo-binder/tree/master/notebooks/demo%20get_specified_length_of_end_of_seq_from_FASTA.ipynb). It can actively by run in your browser by going to [here](https://github.com/fomightez/cl_sq_demo-binder), clicking on the `launch binder` badge, and selecting the 'Demo of script to get specified length of sequence from end of a record in FASTA file' notebook from the listing.


## Related scripts / utilities in my other repositories

To extract sequences from a Clustal-formatted alignment file, see `extract_regions_from_clustal_alignment.py` in ['Sequencework/alignment-utilities' code repository](https://github.com/fomightez/sequencework/tree/master/alignment-utilities/).

If you are gathering some details or other information related to the sequences, you probably want to see my ['Sequencework/Extract_Details_or_Annotation' code repository](https://github.com/fomightez/sequencework/tree/master/Extract_Details_or_Annotation/).

See `drafting_extracting_mito_chr_function.py`. 

[`get_chromosomal_coordinates_as_FASTA.py` in my 'YeastMine' code repository](https://github.com/fomightez/yeastmine) actually gets subsequences from the full genomic sequence of S. cerevisiae S288C sequences. See also these where a subsequences is retrived provided id info for the gene/protein:

* [`get_protein_seq_as_FASTA.py` in my 'YeastMine' code repository](https://github.com/fomightez/yeastmine)
* [`get_gene_genomic_seq_as_FASTA.py` in my 'YeastMine' code repository](https://github.com/fomightez/yeastmine)
)

## Scripts / Resources by others

- >"Fastest way to extract some seqs from a BIG fasta file based on header?" - SOURCE: https://twitter.com/macmanes/status/1045659290919481344

  Good folks with lots of ideas in that response thread.
  
- [SeqKit - a cross-platform and ultrafast toolkit for FASTA/Q file manipulation](https://github.com/shenwei356/seqkit) by shenwei356 has under category 'basic', the commands 'subseq', 'sliding', and 'faidx' that look related. There's also 'Searching', 'amplicon' that is somewhat related. See [here](https://github.com/shenwei356/seqkit#subcommands) for subcommands listing.

- [Example uses Biopython to make a Pandas dataframe from FASTA sequences](https://www.biostars.org/p/9505933/#9505933)
