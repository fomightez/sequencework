CompareFASTA_or_FASTQ utilities
===============================

Repo for my own computational resources dealing with comparing sequences in FASTA or FASTQ format (both nucleic and protein), plus give me a place to reference other handy resources.

- `roughly_score_relationships_to_subject_seq_pairwise_premsa.py`

This script should probably be here but I was thinking about it in context of comparing sequences in a rough manner prior to making a full sequence alignment and put it there


Related utilities in my other repositories
------------------------------------------

- `roughly_score_relationships_to_subject_seq_pairwise_premsa.py`

This script should probably be here but I was thinking about it in context of comparing sequences in a rough manner prior to making a full sequence alignment and put it there


  Often when you use an alignment of entire chromosome to extract a collinear region to identify orthologs sometimes the ends aren't exactly perfect. You may wish to look at the sequence just beyond the edge of collected sequence, and so you can use the end of that extracted seqeunce to search that pattern agains the original source to get the downstream sequence. I made a script to do this post-alignment processing. It is called `get_seq_following_seq_from_multiFASTA.py` and it can be found in my ['Sequencework/Extract_from_FASTA' code repository](https://github.com/fomightez/sequencework/tree/master/Extract_from_FASTA).
  
