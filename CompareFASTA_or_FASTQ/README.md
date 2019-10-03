CompareFASTA_or_FASTQ utilities
===============================

Repo for my own computational resources dealing with comparing sequences in FASTA or FASTQ format (both nucleic and protein), plus give me a place to reference other handy resources.

- `compare_organisms_in_two_files_of_fasta_entries.py`

  See my comparing lists resources for related items.
  
- `score_differences_between_sequences_by_pairwise_alignment.py`

  This script takes sequences in FASTA format (single multi-FASTA file) and scores the different to each of the others. It makes a matrix of the differences.

  This is located [here](https://github.com/fomightez/sequencework/tree/master/alignment-utilities) in my [alignment-utlities sub-repo](https://github.com/fomightez/sequencework/tree/master/alignment-utilities).

  This script should probably be here in this sub-repo; however, I was thinking about it in context of comparing sequences in a rough manner prior to making a full sequence alignment and put it there.

- `roughly_score_relationships_to_subject_seq_pairwise_premsa.py`

  This script takes sequences in FASTA format (single multi-FASTA file) and makes quick assessment of similarity of first sequence to each of the others. It makes a matrix of the differences.

  This is located [here](https://github.com/fomightez/sequencework/tree/master/alignment-utilities) in my [alignment-utlities sub-repo](https://github.com/fomightez/sequencework/tree/master/alignment-utilities).

  This script should probably be here in this sub-repo; however, I was thinking about it in context of comparing sequences in a rough manner prior to making a full sequence alignment and put it there.


Related utilities in my other repositories
------------------------------------------
- `matches_a_patmatch_pattern.py`

  This script can take a sequence in FASTA format or as a text string and tell if it contains a match to a pattern in PatMatch syntaz. *It doesn't expect multi-FASTA file entries though, and will only use the first sequence if one is used as input.* Just reports if matches or not. Really meant to compare a general pattern in PatMatch syntax.  
   This is located [here](https://github.com/fomightez/sequencework/tree/master/patmatch-utilities) in my [patmatch-utlities sub-repo](https://github.com/fomightez/sequencework/tree/master/patmatchutilities).
  
  See my [patmatch-binder](https://github.com/fomightez/patmatch-binder) if you need to locate matches in FASTA sequences and learn the details.

- `score_differences_between_sequences_by_pairwise_alignment.py`

  This script takes sequences in FASTA format (single multi-FASTA file) and scores the different to each of the others. It makes a matrix of the differences.

  This is located [here](https://github.com/fomightez/sequencework/tree/master/alignment-utilities) in my [alignment-utlities sub-repo](https://github.com/fomightez/sequencework/tree/master/alignment-utilities).

  This script should probably be here in this sub-repo; however, I was thinking about it in context of comparing sequences in a rough manner prior to making a full sequence alignment and put it there.


- `roughly_score_relationships_to_subject_seq_pairwise_premsa.py`

  This script takes sequences in FASTA format (single multi-FASTA file) and makes quick assessment of similarity of first sequence to each of the others. It makes a matrix of the differences.

  This is located [here](https://github.com/fomightez/sequencework/tree/master/alignment-utilities) in my [alignment-utlities sub-repo](https://github.com/fomightez/sequencework/tree/master/alignment-utilities).

  This script should probably be here in this sub-repo; however, I was thinking about it in context of comparing sequences in a rough manner prior to making a full sequence alignment and put it there.
  
  
- `report_diff_between_two_seq_strings.py`

  This script directly compares two Python strings using an approach meant for aligning biological strings. This is like a single, direct version of my script `score_differences_between_sequences_by_pairwise_alignment.py`.
  
  This is located [here](https://github.com/fomightez/sequencework/tree/master/Compare_biological_seq_strings) in my [compare_biological_seq_strings utilities sub-repo](https://github.com/fomightez/sequencework/tree/master/Compare_biological_seq_strings).
  
  Put there because this sub-repo is meant for when sequences already in FASTA or FASTQ format and I had occaision to need to compare sequences that are strings in dataframe cells that were pulled out has hits by PatMatch.
  
