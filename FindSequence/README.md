# Collection of scripts to find sequence elements in sequence data

# The scripts

* find_sequence_element_occurrences_in_sequence.py
> smaller sequence and large sequence --> occurrences of EXACT matches of small sequence in larger one

See [the demo notebook](https://nbviewer.jupyter.org/github/fomightez/sequencework/blob/master/FindSequence/demo%20find_sequence_element_occurrences_in_sequence%20script.ipynb) for use. (That notebook can easily be uploaded and actively run once you launch a binder instance by pressing the `launch binder` button [here](https://github.com/fomightez/qgrid-notebooks).)

For searching for ambiguous motifs or mismatches, see below about [PatMatch](https://github.com/fomightez/patmatch-binder).



# Related

- If looking to retrieve a sequence from a database, you want the [RetrieveSeq](https://github.com/fomightez/sequencework/tree/master/RetrieveSeq) directory in this repo.


# Related items by others where I have added some utilities for associating the data with Python

- The venerable PatMatch allows much more flexibility beyond exact matches to nucleotide and protein sequences and motifs, and although it doesn't directly play well with Python, I have added a utility script that helps with that.

See about the command line-based version and links to online implementations at [my repo, patmatch-binder](https://github.com/fomightez/patmatch-binder). It includes active demonstrations in launchable Jupyter notebooks of running the command line-based version and getting the data into Python. Go [there](https://github.com/fomightez/patmatch-binder) and click `launch binder` for that.


# Related items by others

- Came accross [this barebones, related example](https://twitter.com/strnr/status/986215833453113344) (see the extended discussion [here](https://twitter.com/strnr/status/986167127941042177)) after I had made my `find_sequence_element_occurrences_in_sequence.py` script for finding exact matches:

```
# "essentially find the location of a particular sequence." - by Stephen Turner
# source: https://twitter.com/strnr/status/986215833453113344
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)
vmatchPattern("GCGATCGC", Hsapiens)
```
