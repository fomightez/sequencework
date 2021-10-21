# Collection of scripts to find sequence elements in sequence data

# The scripts

* find_sequence_element_occurrences_in_sequence.py
> smaller sequence and large sequence --> occurrences of EXACT matches of small sequence in larger one

See [the demo notebook](https://nbviewer.jupyter.org/github/fomightez/sequencework/blob/master/FindSequence/demo%20find_sequence_element_occurrences_in_sequence%20script.ipynb) for use. (That notebook can easily be uploaded and actively run once you launch a binder instance by pressing the `launch binder` button [here](https://github.com/fomightez/qgrid-notebooks).)

For searching for ambiguous sequences/motifs or mismatches, see [below](#related-items-by-others-where-i-have-added-some-utilities-for-associating-the-data-with-python) about [PatMatch](https://github.com/fomightez/patmatch-binder). However, if you know regular expression format you can indeed use the script `find_sequence_element_occurrences_in_sequence.py` by supplying it text strings formatted in that form to search in the sequence.

Additional **tips for using `find_sequence_element_occurrences_in_sequence.py`**:
- If you are using this with protein sequences and thus the reverse strand search is moot/confusing, you can run the following where you downloaded the script to comment out the lines running the reverse strand processing: 

  !sed -i '371s/.*/    #/' find_sequence_element_occurrences_in_sequence.py  
  !sed -i '372s/.*/    #/' find_sequence_element_occurrences_in_sequence.py

(That replacement works with [this version](https://github.com/fomightez/sequencework/commit/effaf12354468c9b0288f5c3eca129192f70d350).Unfortunately once I change the lines this will break!!)
- If you know regular expression format you can indeed use the script `find_sequence_element_occurrences_in_sequence.py` by supplying it text strings formatted in that form to search in the sequence. For example, I was able to use `element = "F[ST]P.*?G.[LF].*?[ST].Y"` for using the imported main function

# Related

- If looking to retrieve a sequence from a database, you want the [RetrieveSeq](https://github.com/fomightez/sequencework/tree/master/RetrieveSeq) directory in this repo.

- If you want to retrieve the sequence following a match to a sequence, you want the [Extract_from_FASTA](https://github.com/fomightez/sequencework/tree/master/Extract_from_FASTA) directory in this repo

- If you want to delete sequence following a match to a sequence, you want the [AdjustFASTA_or_FASTQ](https://github.com/fomightez/sequencework/tree/master/AdjustFASTA_or_FASTQ) directory in this repo

- If you just want the coordinates spanned / position of a matched sequence, you may want to use `report_coordinates_for_seq_within_FASTA.py` script in my [collection of scripts to Extract Details or Annotations](https://github.com/fomightez/sequencework/tree/master/Extract_Details_or_Annotation) sub-repo. This potentially complements `find_sequence_element_occurrences_in_sequence.py`.


# Related items by others where I have added some utilities for associating the data with Python

- The venerable PatMatch allows much more flexibility beyond exact matches to nucleotide and protein sequences and motifs (and is easier to use than regular expression-based search terms because tailored to residues), and although it doesn't directly play well with Python, I have made a utility script that helps with that.  
See about the command line-based version and links to online implementations at [my repo, patmatch-binder](https://github.com/fomightez/patmatch-binder). It includes active demonstrations in launchable Jupyter notebooks of running the command line-based version and getting the data into Python. Go [there](https://github.com/fomightez/patmatch-binder) and click `launch binder` for that.  
Related to this, are the scripts I have placed in my [patmatch-utilities](https://github.com/fomightez/sequencework/tree/master/patmatch-utilities) sub-repo.


# Related items by others

- Came accross [this barebones, related example](https://twitter.com/strnr/status/986215833453113344) (see the extended discussion [here](https://twitter.com/strnr/status/986167127941042177)) after I had made my `find_sequence_element_occurrences_in_sequence.py` script for finding exact matches:

```
# "essentially find the location of a particular sequence." - by Stephen Turner
# source: https://twitter.com/strnr/status/986215833453113344
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)
vmatchPattern("GCGATCGC", Hsapiens)
```
https://twitter.com/BQPMalfoy/status/1202902187174707200     December 2019
>"Finaly reverse complement grep exist !"
[rcgrep: Search text files for DNA sequences and their reverse complements](https://github.com/dib-lab/rcgrep)  <--- this is from DIBs lab when traced back and must be what was I was thinking of when I thought screen migth have it

[SeqKit - a cross-platform and ultrafast toolkit for FASTA/Q file manipulation](https://github.com/shenwei356/seqkit) by shenwei356 has under category 'searching', the commands 'grep', 'locate', 'amplicon', and 'fish' that look related. See [here](https://github.com/shenwei356/seqkit#subcommands) for subcommands listing.


[Cool R script to colorfully annotate adapter sequences in a sequence - could be useful to adapt for displaying PCR primers, etc](https://twitter.com/clintcodesbio/status/1339947174239612929); R script [here](https://gitlab.com/gringer/bioinfscripts/-/blob/master/read_annotator.r) that highlights subseqeunces

