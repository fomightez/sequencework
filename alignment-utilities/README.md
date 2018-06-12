# Alignment utilities

Repo for my own computational resources dealing with aligning and alignments (both nucleic and protein), plus give me a place to reference other handy resources.

* add_actual_position_numbering_to_sequence_in_top_line_of_clustal_output.py
> alignment text --> alignment text with top sequence in alignment blocks annotated with numbering according to actual contiguous sequence (i.e., excluding gaps)

info here

* BIOALIGNAMER (* NOT HERE YET BUT BELONGS HERE*)
>alignment text with cryptic identifiers for organisms --> alignment text with common names as ids

I'd need to see if it needs updating/fixing (for example, I think it is purely 2.7) but I have one of the more fleshed out Python scripts I wrote that corrects names in alignments to add the common name automaigically. Created for facilating pipelines for submitting tasks to the GARLI webserver.
I intended to put up a web-served version of it [here](http://fomightez.pythonanywhere.com/) but I recall I hadn't figured out how to get the toggles to work at the bottom at the time for controlling the command line parameters. (I knot I got some text entry forms working later, for example see [here](http://fomightez.pythonanywhere.com/ammonium_screen/) and [here](http://fomightez.pythonanywhere.com/spartan_fixer/) and [here](http://fomightez.pythonanywhere.com/triplet_freq/), and so I probably could figure it out now.)

* pretty_msa_maker_from_clustal_nucleic.py (* NOT HERE YET BUT SHOULD BE HERE*)
> alignment in text form --> alignment fully annotated with features and numbers in vector graphics form

info here

* extract_regions_from_clustal_alignment.py
> alignment of large region ---> alignment of a portion plus sequence of ungapped corresponding block/

Intended to be used then aligning the entire chromosome along full length to insure corresponding sub-regions best aligned (assuming no rearramgements).

info here.


* reverse_complement_of_clustal_alignment.py
> alignment ---> alignment with same sequences in reverse complement form

Intended to be used after `extract_regions_from_clustal_alignment.py` for sequences on the other (a.k.a. reverse) strand because some might be on other strand when extract from entire chromosome and sometimes still it would make more sense to view the extracted sequence in the other direction. For example, if typically used to looking at the regions in one direction, etc.

info here.


* calculate_cons_for_clustal_nucleic.py
> nucleic alignment without consensus line ---> alignment with consensus line

info here.


Utilities by others
------------------

- [Adding unaligned sequences into an existing alignment using MAFFT and LAST](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3516148/)

I have several listed [here](http://proteopedia.org/wiki/index.php/User:Wayne_Decatur/Sequence_analysis_tools) such as those that interconvert/reformat alignments in various forms. Just search `align` on that page to step through examining them. In particular note:

- Under 'Installable software for fine-tuning sequence alignments' there, I discuss using SEQOTRON on my work Mac and  `cons` alignment consensus program and many others at EMBOSS explorer website to put in conservation. Along the way it looks like it adds numbering similar to Mview where it is just numbering for total number on line INCLUDING GAPS. (Not to self, see my `edited t-coffee alignment for true XXXXX XXXXX.md` for example.)
