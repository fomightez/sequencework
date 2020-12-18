# Annotation utilities

Scripts
-----------

- `measure_intergenic_regions_in_mito_annotations.py`

>annotation file and corresponding sequence file for a fungal mitochondrial genome > list of intergenic gap sizes

There is a [demo notebook for this script in this repo](https://github.com/fomightez/cl_sq_demo-binder). Launch a binder session from there and select 'Demo of script to get intergenic gap sizes from annotation file' to run it actively.  The particular notebook can be viewed statically, nicely displayed [here](https://nbviewer.jupyter.org/github/fomightez/cl_sq_demo-binder/blob/master/notebooks/Demo%20of%20script%20to%20get%20intergenic%20gap%20sizes%20from%20annotation%20file.ipynb).

Utilities by Others
-------------------

[gencode_regions - Extract 3'UTR, 5'UTR, CDS, Promoter, Genes, Introns, Exons from GTF files](https://github.com/saketkc/gencode_regions#confused-about-exons-and-utrs) - Even has several yeast examples, such as [sacCerR64.ipynb](https://github.com/saketkc/gencode_regions/blob/master/notebooks/sacCerR64.ipynb) [S_cerevisiae-BY4741-stanford.ipynb](https://github.com/saketkc/gencode_regions/blob/master/notebooks/S_cerevisiae-BY4741-stanford.ipynb), besides other organisms.

[Cool R script to colorfully annotate adapter sequences in a sequence - could be useful to adapt for displaying PCR primers, etc](https://twitter.com/clintcodesbio/status/1339947174239612929); R script [here](https://gitlab.com/gringer/bioinfscripts/-/blob/master/read_annotator.r)

See also
-------

The contents of [the Adjust_Annotation sub-repositor](https://github.com/fomightez/sequencework/tree/master/Adjust_Annotation) would probably be best put under this topic but predates this sub-directory. Currently it houses two scripts and points at some other resources:
* fix_lsu_rRNA_annotation_in_gff_resulting_from_mfannot.py
* makes_length_annotation_file.py

See `check_for_omega_intron.py` under [the omega-presence sub-repository](https://github.com/fomightez/sequencework/tree/master/omega-presence). Technically that could be used to make annotations but it isn't presently, and doesn't rely on gff3 annoation file, and so placed in its own sub-directory for now.
