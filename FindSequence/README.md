# Collection of scripts to find sequences in other sequences

# The scripts

* find_sequence_element_occurrences_in_sequence.py
> smaller sequence and large sequence --> occurrences of exact matches of small sequence in larger one

See [the demo notebook](https://nbviewer.jupyter.org/github/fomightez/sequencework/blob/master/FindSequence/demo%20find_sequence_element_occurrences_in_sequence%20script.ipynb) for use.




# Related

- If looking to retrieve a sequence from a database, you want the [RetrieveSeq](https://github.com/fomightez/sequencework/tree/master/RetrieveSeq) directory in this repo.

# Related items by others

- Came accross [this barebones, related example](https://twitter.com/strnr/status/986215833453113344) after I had made my `find_sequence_element_occurrences_in_sequence.py` script for finding exact matches:

```
# "essentially find the location of a particular sequence." - by Stephen Turner
# source: https://twitter.com/strnr/status/986215833453113344
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)
vmatchPattern("GCGATCGC", Hsapiens)
```

