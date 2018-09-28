# Alignment utilities

Repo for my own computational resources dealing with aligning and alignments (both nucleic and protein), plus give me a place to reference other handy resources.

* add_actual_position_numbering_to_sequence_in_top_line_of_clustal_output.py
> alignment text --> alignment text with top sequence in alignment blocks annotated with numbering according to actual contiguous sequence (i.e., excluding gaps)

info here

* BIOALIGNAMER (* NOT HERE YET BUT BELONGS HERE*)
>alignment text with cryptic identifiers for organisms --> alignment text with common names as ids

I'd need to see if it needs updating/fixing (for example, I think it is purely 2.7) but I have one of the more fleshed out Python scripts I wrote that corrects names in alignments to add the common name automaigically. Created for facilating pipelines for submitting tasks to the GARLI webserver.
I intended to put up a web-served version of it [here](http://fomightez.pythonanywhere.com/) but I recall I hadn't figured out how to get the toggles to work at the bottom at the time for controlling the command line parameters. (I know I got some text entry forms working later, for example see [here](http://fomightez.pythonanywhere.com/ammonium_screen/) and [here](http://fomightez.pythonanywhere.com/spartan_fixer/) and [here](http://fomightez.pythonanywhere.com/triplet_freq/), and so I probably could figure it out now.)

* pretty_msa_maker_from_clustal_nucleic.py (* NOT HERE YET BUT SHOULD BE HERE*)
> alignment in text form --> alignment fully annotated with features and numbers in vector graphics form

info here

* extract_regions_from_clustal_alignment.py
> alignment of large region ---> alignment of a portion, plus saves sequence of ungapped corresponding block as well

Intended to be used then aligning the entire chromosome along full length, end-to-end to insure corresponding sub-regions best aligned (assuming no rearrangements). Under 'Whole genome alignment', about 4/5 through set of slides, [here](http://www.utdallas.edu/~pradiptaray/teaching/1_mult_seq_align.pdf) refers to this method as '• Identify “collinear” (orthologous) regions or blocks and perform piecewise alignment'.)

info here.


* reverse_complement_of_clustal_alignment.py
> alignment ---> alignment with same sequences in reverse complement form

Intended to be used after `extract_regions_from_clustal_alignment.py` for sequences on the other (a.k.a. reverse) strand because some might be on other strand when extract from entire chromosome and sometimes still it would make more sense to view the extracted sequence in the other direction. For example, if typically used to looking at the regions in one direction, etc.

info here.

* score_columns_in_clustal_msa.py
>alignment ---> two scores assessing conservation of an MSA overall in a quick-n-dirty manner

Intended to be used when comparing separate alignments of related but distinct genetic elements or protein paralogs/orthologs/gene family members/gene cluster members in order to objectively assess which of the related elements are more highly conserved overall.

info here

* calculate_cons_for_clustal_nucleic.py
> nucleic alignment without consensus line ---> alignment with consensus line

info here.


* roughly_score_relationships_to_subject_seq_pairwise_premsa.py
> sequences in FASTA format (single multi-FASTA file)  ---> quick assessment of similarity of first sequence to each of the others

Meant to be run prior a multi-sequence alignment task in order to tell quickly what is similar to one another and maybe leave out redundant sequences if you are looking to reduce the submitted alignment job.

Produces a ranked list where the highest scoring one is closest to the first sequence in the provided sequence.

The core of this is based on [this script](https://github.com/berrisfordjohn/adding_stats_to_mmcif/blob/master/adding_stats_to_mmcif/pairwise_align.py) by [John Berrisford](https://github.com/berrisfordjohn).

There is a [demo notebook for this script in this repo](https://github.com/fomightez/sequencework/blob/master/alignment-utilities/demo%20roughly_score_relationships_to_subject_seq_pairwise_premsa.ipynb) that can be viewed nicely displayed [here](https://nbviewer.jupyter.org/github/fomightez/sequencework/blob/master/alignment-utilities/demo%20roughly_score_relationships_to_subject_seq_pairwise_premsa.ipynb).

CAVEAT:  
The original intention was to align the entire length of provided sequences but for moderate- or large-sized sequences (>5 kb) this is not possible unless you have a substantial amount of computing power due to demanding memory resources from Biopython's aligmnet algorithmn (see [here](https://github.com/biopython/biopython/pull/1655)). Thus, it was modified to sample what seems to work when using the script in Jupyter sessions launched from MyBinder.org via [here](https://github.com/fomightez/blast-binder). Therefore, the reported score may not represent the entire range of the provided sequences but that is reported back if that is the case.

Related utilities in my other repositories
------------------------------------------

- `collapse_large_unknown_blocks_in_DNA_sequence.py`

  In order to better judge some sequences extracted from alignments against a chromosome/genome, I've found that reducing apparently arbitrarily sized repeats of unknown nucleotides, represented as uninterrupted repeats of `N`s, such as "NNNNNNNNNNNNNNN" for example, can help in assessing the sequences extracted from a much larger alignment or help prepare them for aligning again, especially when additional knowledge, such as matches in flanking sequence seem to suggest the spacing doesn't match the number of unknown nucleotides shown in the represtation. In other words, often in assemblies generated from Illumina pipelines 'educated guesses' or abtirary numbers (such as 50) `N`s in a row will be introduced for smallish 'gaps' (i.e., ones where software has concluded an intact scaffold is still indicated) in the sequence assembly and reducing these to a smaller size can sometimes make alignments of sub-elements in different rows more obvious. It is important to be clear this was done and keep both sets of data. I have made a script that does this reduction in the size of the blocks; it is called `collapse_large_unknown_blocks_in_DNA_sequence.py` and it can be found in my ['ConvertSeq' code repository](https://github.com/fomightez/sequencework/tree/master/ConvertSeq).


Utilities by others
------------------

- [Adding unaligned sequences into an existing alignment using MAFFT and LAST](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3516148/)

I have several listed [here](http://proteopedia.org/wiki/index.php/User:Wayne_Decatur/Sequence_analysis_tools) such as those that interconvert/reformat alignments in various forms. Just search `align` on that page to step through examining them. In particular note:

- Under 'Installable software for fine-tuning sequence alignments' there, I discuss using SEQOTRON on my work Mac and  `cons` alignment consensus program and many others at EMBOSS explorer website to put in conservation. Along the way it looks like it adds numbering similar to Mview where it is just numbering for total number on line INCLUDING GAPS. (Not to self, see my `edited t-coffee alignment for true XXXXX XXXXX.md` for example.)
