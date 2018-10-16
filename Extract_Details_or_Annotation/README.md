# Collection of scripts to Extract Details or Annotations

Repo for my computational resources dealing with gathering information or some details from sequence files and not the actual sequence. If you are seeking to include the sequence, you should probably look into my other repositories, see ['Related' below](https://github.com/fomightez/sequencework/tree/master/Extract_Details_or_Annotation#related). 


# The scripts

* mine_mito_features_from_transcriptome.py
> transcriptome in FASTA format --> dataframe of mitochdondrial chromsome genes/features

Takes a transcriptome file of entries in FASTA format and mines the 
mitochondrial details from definition lines like below:


```
>Q0065 cdna chromosome:R64-1-1:Mito:13818:21935:1 gene:Q0065 gene_biotype:protein_coding transcript_biotype:protein_coding gene_symbol:AI4 description:Endonuclease I-SceII; encoded by a mobile group I intron within the mitochondrial COX1 gene; intron is normally spliced by the BI4p maturase but AI4p can mutate to acquire the same maturase activity [Source:SGD;Acc:S000007264]
>Q0143 cdna chromosome:R64-1-1:Mito:51277:51429:1 gene:Q0143 gene_biotype:protein_coding transcript_biotype:protein_coding description:Dubious open reading frame; unlikely to encode a functional protein, based on available experimental and comparative sequence data [Source:SGD;Acc:S000007277]
```

Returns a dataframe with the details for the mito features.
Specifically, in this case systematic_id, start, end, midpoint, gene_symbol
(standard name) if there is one noted

* mine_SUTs_and merge_into_mito_transcripts_df.py 
> SUTs data and transcriptome data --->  dataframe of mitochdondrial chromsome genes/features

Intended to be used in conjuction with `mine_mito_features_from_transcriptome.py`.  
Returns a dataframe with the details for the mito features, now with yeast mitochondrial SUTs added.

* mine_XUTs_and merge_into_mito_transcripts_df.py 
> XUTs data and transcriptome/SUTs data --->  dataframe of mitochdondrial chromsome genes/features

Intended to be used in conjuction with `mine_mito_features_from_transcriptome.py` and `mine_SUT_and merge_into_mito_transcripts_df.py`, specifically after `mine_SUT_and merge_into_mito_transcripts_df.py`.  
Returns a dataframe with the details for the mito features, now with yeast mitochondrial XUTs added.

* report_coordinates_for_seq_within_multiFASTA.py 
> subsequence from a specific sequence in a FASTA file --->  coordinates of the match

There is a [demo notebook for this script in this repo](https://github.com/fomightez/sequencework/blob/master/Extract_Details_or_Annotation/demo%20report_coordinates_for_seq_within_multiFASTA.ipynb) that can be viewed nicely displayed [here](https://nbviewer.jupyter.org/github/fomightez/sequencework/blob/master/Extract_Details_or_Annotation/demo%20report_coordinates_for_seq_within_multiFASTA.ipynb).

Takes a sequence pattern string, a sequence file (FASTA-format), and a record identifier, and reports the start and end coordinates of that sequence within the specified FASTA record.
Because inherently the sequence to search can be provided as a pattern, I am saying 'pattern'; however, the impetus for this script was for when you have a known sequence, that will only occour once in that record and you need the coordinates. The FASTA-formatted sequence file in which to search is assumed by default to be a multi-FASTA, i.e., multiple sequences in the provided file, although it definitely doesn't have to be. More on that below. 

After running the script, two values, the start and the end will be returned, separated by a tab. Redirect the output to a file if using command line version and a file is needed.

Importantly, the coordinates numbering is in 'common' terms where the position numbered one corresponds to the first position. (In other words, thecoordinates returned are in the more conventional form & are not zero-indexed even though behind-the-scenes Python/Biopython is indeed zero-indexed.)

Note that if there is only one record in the specified sequence file, the record id is moot and you can instead provide any string for that parameter as it will be ignored. This makes the script more flexible in cases where sequence files aren't complex as the user doesn't need to provide an actual record id.

The provided sequence pattern will be matched regardless of case, as both the input sequence and pattern to search will be converted behind-the-scenes to lowercase for the comparison. Beyond being insensitive of the case, REGULAR EXPRESSION SEARCH TERM SYNTAX IS ACCEPTABLE in the provided sequence pattern.

It is designed to handle/filter gaps ('dashes') in the provided sequence patterns. The idea being that the known sequence ends may be manually extracted from sequence alignments. This way the user is not wasting time removing the gap indications / dashes from the collected text lines.

Note that why some aspects of this script may seem redundant with my scripts [`find_sequence_element_occurrences_in_sequence.py`](https://github.com/fomightez/sequencework/tree/master/FindSequence) and [`blast_to_df.py`](https://github.com/fomightez/sequencework/tree/master/blast-utilities) that include coordinates in the output tables, those were meant to mainly be used when looking for a table of information about matches or exploring matches that may occur multiples times in many sequences in the supplied file(s). This script is targeted at the case where you know record id and expect only one match to that record with that identifier.

There is a [demo notebook for this script in this repo](https://github.com/fomightez/sequencework/blob/master/Extract_Details_or_Annotation/demo%20report_coordinates_for_seq_within_multiFASTA.ipynb) that can be viewed nicely displayed [here](https://nbviewer.jupyter.org/github/fomightez/sequencework/blob/master/Extract_Details_or_Annotation/demo%20report_coordinates_for_seq_within_multiFASTA.ipynb).

# Related

If you are seeking to include the sequence in the data gathered, you should probably look into my other repositories, such as:

* ['Sequencework/Extract_from_FASTA' code repository](https://github.com/fomightez/sequencework/tree/master/Extract_from_FASTA/)
* ['Sequencework/FindSequence' code repository](https://github.com/fomightez/sequencework/tree/master/FindSequence/)
* ['Sequencework/alignment-utilities' code repository](https://github.com/fomightez/sequencework/tree/master/alignment-utilities/)
* ['Sequencework/RetrieveSeq' code repository](https://github.com/fomightez/sequencework/tree/master/RetrieveSeq/)

Specific, related scripts to note here:

  * `UCSC_chrom_sizes_2_circos_karyotype.py` extracts chromosome sizes from UCSC `chrom.sizes` files and can be found in [my collection of circos-related utility scripts (Python)](https://github.com/fomightez/sequencework/tree/master/circos-utilities).
