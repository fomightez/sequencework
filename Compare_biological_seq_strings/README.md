Compare_biological_seq_strings utilities
===============================

Repo for my own computational resources dealing with comparing sequences (both nucleic and protein) as Python strings, plus give me a place to reference other handy resources.

- `report_diff_between_two_seq_strings.py`

  Takes two Python strings that correspond to biological sequences and reports differences between them as calculated from a global pairwise alignment and optionally reports the length difference of the two unaligned strings. 
  
  The idea for now is that this reported data can be used to gauge amount of changes between matches to a sequence pattern in order to summarize a lot of matches to the pattern in many genomes. However, I am writing it as a general function and so it may be useable in several other contexts.
  
  For now, I have added a very basic demo of this script to this sub-repo as a 'stub', `stub of demo for report_diff_between_two_seq_strings script.ipynb`.



Related utilities in my other repositories
------------------------------------------

See:  

- `score_differences_between_sequences_by_pairwise_alignment.py` in [my alignment-utilities sub-repo](https://github.com/fomightez/sequencework/tree/master/alignment-utilities).
It does what the script `report_diff_between_two_seq_strings.py` does in calculating differences but it does it between every FASTA sequence entry in a multi-entry FASTA sequence.

- the contents of [my CompareFASTA_or_FASTQ sub-repo](https://github.com/fomightez/sequencework/tree/master/CompareFASTA_or_FASTQ). In those cases, the input sequences are in FASTA/Q format. In this current sub-repo, I am just considering Python strings.


Related utilities by others
---------------------------

- There are other ways to score alignments. For now I have been using simplistic assessments based on differences. This can be done much fancier based on probability models, with different scores for sections in between the gaps as well. See [How sequence alignment scores correspond to probability models. Martin C Frith. Bioinformatics, btz576, https://doi.org/10.1093/bioinformatics/btz576](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btz576/5536873) for recent coverage of this, although I don't see an implementation linked to there on my cursory examination.
