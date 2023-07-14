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

- [FastANI](https://github.com/ParBLiSS/FastANI) - FastANI is developed for fast alignment-free computation of whole-genome Average Nucleotide Identity (ANI)  
>"FastANI rapidly compares fasta sequences and will output global identity. It is meant only for relatively similar genomes, and will not create any output if identity is below 75-80%." - from Mensur Dlakic's Biostars answer [here](https://www.biostars.org/p/9509271/#9509767)


- If you are doing massive scales to compare non-rendundant sequences to get a sense of the degree of difference or increase in content of a databases since last time you checked, see this discussion:  
https://twitter.com/BenLangmead/status/1679658259400077314    July 2023
>"Genome folks: saw a tweet like this today? Can we do better than gzip as a way of measuring non-redundant sequence content & similarity, including between huge sequence collections, without gzip's windowing?  Super efficiently? See this by Jessica Bonnie: https://biorxiv.org/content/10.1101/2023.02.02.526837v1 "  
https://twitter.com/jnalanko/status/1679843958229200896     July 2023  
>"Twitter makes this thing sound super novel and creative. Actually, the concept of normalized compression distance (NCD) has been known at least since 2005 [1]. 🧵"  
>"I played around with this for clustering metagenomic samples with Lempel-Ziv compression as a summer intern in 2015. While the initial results were okay, we abandoned the project because it was looking like just a hack to approximate more principled substring distances."
