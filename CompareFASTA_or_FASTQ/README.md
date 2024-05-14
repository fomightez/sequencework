CompareFASTA_or_FASTQ utilities
===============================

Repo for my own computational resources dealing with comparing sequences in FASTA or FASTQ format (both nucleic and protein), plus give me a place to reference other handy resources.

- `compare_organisms_in_two_files_of_fasta_entries.py`

  See my resources for comparing lists among [my text_mining repository](https://github.com/fomightez/text_mining) for related items.
  
  &nbsp;<p></p>
  
- `score_differences_between_sequences_by_pairwise_alignment.py`

  This script takes sequences in FASTA format (single multi-FASTA file) and scores the different to each of the others. It makes a matrix of the differences.

  This is located [here](https://github.com/fomightez/sequencework/tree/master/alignment-utilities) in my [alignment-utlities sub-repo](https://github.com/fomightez/sequencework/tree/master/alignment-utilities).

  This script should probably be here in this sub-repo; however, I was thinking about it in context of comparing sequences in a rough manner prior to making a full sequence alignment and put it there.

&nbsp;<p></p>

- `roughly_score_relationships_to_subject_seq_pairwise_premsa.py`

  This script takes sequences in FASTA format (single multi-FASTA file) and makes quick assessment of similarity of first sequence to each of the others. It makes a matrix of the differences.

  This is located [here](https://github.com/fomightez/sequencework/tree/master/alignment-utilities) in my [alignment-utlities sub-repo](https://github.com/fomightez/sequencework/tree/master/alignment-utilities).

  This script should probably be here in this sub-repo; however, I was thinking about it in context of comparing sequences in a rough manner prior to making a full sequence alignment and put it there.

&nbsp;<p></p>

Related utilities in my other repositories
------------------------------------------
- `matches_a_patmatch_pattern.py`

  This script can take a sequence in FASTA format or as a text string and tell if it contains a match to a pattern in PatMatch syntaz. *It doesn't expect multi-FASTA file entries though, and will only use the first sequence if one is used as input.* Just reports if matches or not. Really meant to compare a general pattern in PatMatch syntax.  
   This is located [here](https://github.com/fomightez/sequencework/tree/master/patmatch-utilities) in my [patmatch-utlities sub-repo](https://github.com/fomightez/sequencework/tree/master/patmatchutilities).
  
  See my [patmatch-binder](https://github.com/fomightez/patmatch-binder) if you need to locate matches in FASTA sequences and learn the details.

&nbsp;<p></p>

- `score_differences_between_sequences_by_pairwise_alignment.py`

  This script takes sequences in FASTA format (single multi-FASTA file) and scores the different to each of the others. It makes a matrix of the differences.

  This is located [here](https://github.com/fomightez/sequencework/tree/master/alignment-utilities) in my [alignment-utlities sub-repo](https://github.com/fomightez/sequencework/tree/master/alignment-utilities).

  This script should probably be here in this sub-repo; however, I was thinking about it in context of comparing sequences in a rough manner prior to making a full sequence alignment and put it there.


- `roughly_score_relationships_to_subject_seq_pairwise_premsa.py`

  This script takes sequences in FASTA format (single multi-FASTA file) and makes quick assessment of similarity of first sequence to each of the others. It makes a matrix of the differences.

  This is located [here](https://github.com/fomightez/sequencework/tree/master/alignment-utilities) in my [alignment-utlities sub-repo](https://github.com/fomightez/sequencework/tree/master/alignment-utilities).

  This script should probably be here in this sub-repo; however, I was thinking about it in context of comparing sequences in a rough manner prior to making a full sequence alignment and put it there.
  
  &nbsp;<p></p>
  
- `report_diff_between_two_seq_strings.py`

  This script directly compares two Python strings using an approach meant for aligning biological strings. This is like a single, direct version of my script `score_differences_between_sequences_by_pairwise_alignment.py`.
  
  This is located [here](https://github.com/fomightez/sequencework/tree/master/Compare_biological_seq_strings) in my [compare_biological_seq_strings utilities sub-repo](https://github.com/fomightez/sequencework/tree/master/Compare_biological_seq_strings).
  
  Put there because this sub-repo is meant for when sequences already in FASTA or FASTQ format and I had occaision to need to compare sequences that are strings in dataframe cells that were pulled out has hits by PatMatch.
  
  
Related utilities by others
----------------------------

- [FastANI](https://github.com/ParBLiSS/FastANI) - FastANI is developed for fast alignment-free computation of whole-genome Average Nucleotide Identity (ANI)  
   >"FastANI rapidly compares fasta sequences and will output global identity. It is meant only for relatively similar genomes, and will not create any output if identity is below 75-80%." - from Mensur Dlakic's Biostars answer [here](https://www.biostars.org/p/9509271/#9509767)

- [degenotate - Annotate degeneracy of sites in coding regions of a genome](https://github.com/harvardinformatics/degenotate) Learned of it from [this Biostar's question seeking help to calculate degreee of degeneracy](https://www.biostars.org/p/9594758/#9594839).
  >"degenotate takes as input either a genome FASTA file and a corresponding annotation file (GFF or GTF) OR file or directory of files that contain coding sequences in FASTA format and outputs a bed-like file that contains the degeneracy score (0-, 2-, 3-, or 4-fold) of every coding site."
Vaguely related to my 'score' utilities featured on this page.
