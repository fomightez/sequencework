# patmach-utilities

Utility scripts for working with command line-based PatMatch, which is launchable [here](https://github.com/fomightez/patmatch-binder).

# The scripts

* patmatch_results_to_df.py
> PatMatch results --> dataframe of data for use in Python

Takes output from command line PatMatch and makes a dataframe from it for use with Python.

Verified compatible with both Python 2.7 and Python 3.6.

Written to run from command line or pasted/loaded inside a Jupyter notebook cell.  
The main ways to run the script are demonstrated in [a series of notebooks illustrating using PatMatch with Python](https://github.com/fomightez/patmatch-binder).

Assumes PatMatch is run with the `-c` flag to search on both strands.

Example calls to run script from command line:
```
python patmatch_results_to_df.py results_file.txt

python patmatch_results_to_df.py test.out --pattern DDWDWTAWAAGTARTADDDD -name promoter

python patmatch_results_to_df.py protein.out --pattern TYEETGLQGHPS -name motif -p
```

(Alternatively, upload the script to a Jupyter environment and use `%run patmatch_results_to_df.py results_file.txt` in a Python-backed notebook to run the example.)



**The pickle format ending in extension `.pkl` was chosen as the default output, despite not being human-readable, in order to store data efficiently since searches of entire gemomes has the potential to generate a lot of hits.** If you need to convert to a text-readable form, you can do the following with Jupyter or IPython where the pickled data is located to add a form saved as tab separated text:

    import pandas as pd
    df = pd.read_pickle("patmatch_pickled_df.pkl")
    df.to_csv('patmatch_data.tsv', sep='\t',index = False) 



#### For running in a Jupyter notebook:

To use this script after pasting or loading into a cell in a Jupyter notebook, in the next cell define the URL and then call the main function similar to below:
```
my_pattern= "GAATTC"
df = patmatch_results_to_df("results_file.txt", pattern=my_pattern, name="EcoRI")
```
See [here](https://github.com/fomightez/patmatch-binder) for notebooks demonstrating use within a Jupyter notebook.


#### Using patmatch_results_to_df with large genomes
When I used PatMatch with the megagenome of pine sugar on Cyverse, it was including in its output several warnings like below:

`Warning: recSearchFile: Record longer than buffer size (10000000) has been split\n`

It probably is issuing a warning so you know that it may not show a match if it happens to match where the split was chosen. They way I dealt with this was to remove the warnings from the input to `patmatch_results_to_df.py` or `patmatch_results_to_df()`, see example [here](https://github.com/fomightez/ptmbr-accompmatz/blob/master/notebooks/PatMatch%20use%20on%20the%20largest%20genome%20project%20to%20date.ipynb), because otherwise the warnings cause errors when the lines are processed as typical output. There error is caused because there are no start or end coordinates to extract. However, instead of having `patmatch_results_to_df.py` recognize the warning and handle it, I thought it best the user has to deal with the warnings so they are aware.


* matches_a_patmatch_pattern.py
> pattern in PatMatch syntax and a sequence --> report if they match

Takes pattern in PatMatch syntax and a sequence in text or FASTA format and reports if they match.

There is a [demo notebook for this script that runs in my patmatch-binder](https://github.com/fomightez/patmatch-binder). You can launch the series from [here](https://github.com/fomightez/patmatch-binder) and then selecting from the list at the bottom of the index page undert 'Additional topics: Technical' to go to the 'Demo of script to check for a match to a sequence pattern in PatMatch syntax' page. The direct link to a nicely-rendered, static version of that page is [here](https://nbviewer.jupyter.org/github/fomightez/ptmbr-accompmatz/blob/master/notebooks/Demo%20of%20script%20to%20check%20for%20a%20match%20to%20a%20sequence%20pattern%20in%20PatMatch%20syntax.ipynb).

*Details*: Takes a sequence pattern in PatMatch syntax, and checks if a provided sequence (file name for sequence in FASTA format or string text can be  provided) contains a match. It reports True or False depending on that
assessment. Optionally, it can be restricted to checking if the provided 
sequence is a match to a sequence pattern. Note that the point, i.e., 
checking that the provided sequence matches entirely to the pattern and is 
therefore a specific instance of the more general pattern, was the 
original impetus for writing this script. However, it seemed I might write a
more general script and then include an option to do that.
When the residue type is nucleic, both strands will be searched.
Only concerned with first sequence if a multi-sequence file is provided. Loop
over sequences and then call script with each if you need to do multiple.

Currently, designed with comparing short sequence strings in mind. For 
PatMatch syntax, see [here](https://www.yeastgenome.org/nph-patmatch#examples).
For those familiar with prepration of files for use with PatMatch, you don't 
have to apply `unjustify.pl` yourself to your sequences. That preparation of 
removing line endings will all be handled internally meaning 
STANDARD FASTA FILE FORMAT IS ACCEPTABLE. 

Note that with 'match_over_entirety', the pattern itself must be same length 
in characters as the expected matching sequence, and therefore it is 
recommended to avoid using it with fancier patterns using symbols as such 
patterns will likely cause this not to work correctly. If using such fancier 
search patterns, you'll want to screen the length of the sequence realtive the 
possible outcomes prior and combine the results of the assessment of match 
into your interpretation of what the report means. Or you'd need to edit the 
script.

Meant to be a 'yes' or 'no' answer. If you want to know about the location of 
the match or matches you'll want to use `patmatch_results_to_df.py` from [this same sub-repo](
https://github.com/fomightez/sequencework/tree/master/patmatch-utilities).
It actually uses that script to do the work here.

Note this is conceptually, vaguely related to my script 
`find_sequence_element_occurrences_in_sequence.py` found [here]( 
https://github.com/fomightez/sequencework/tree/master/FindSequence) except that
handles Regular Expression (REGEX) syntax, and here I want to be able to 
handle more biological contexts without having to add many sets to the regex 
code. By using PatMatch, I get the biological context without having to 
re-implement it.

Note on PatMatch syntax:  
The letters themselves match IUPAC ambiguity codes, see 'Nucleotide ambiguity code' list [here](https://www.dnabaser.com/articles/IUPAC%20ambiguity%20codes.html). However, pattern syntax for PatMatch goes beyond that as it can handle broader constraints than just the letters, as shown [here](https://www.yeastgenome.org/nph-patmatch#examples). See [here](https://www.biostars.org/p/264212/#264218) for a discussion of alignment tools that also allow the IUPAC codes. 



Related
-------

[patmatch-binder](https://github.com/fomightez/patmatch-binder) - for running command line-based PatMatch on sequences of your choice in your browser without need for downloads, installations, or maintenance.


Related Utilities by Others
---------------------------

Go the other direction .... [Biostars post 'Python3 or Ubuntu, not perl: Have Primers with Degenerate Bases, Need tool or way to List all Possible Nucleotide Sequences'](https://www.biostars.org/p/9589753/) containes resources to do the reverse. -- One of them is [Generate all the possible combinations of a Degenerate DNA/RNA sequence in FASTA format using Python - Degenerecy..py by Chamara Janaka Bandara](https://www.researchgate.net/publication/354451685_Generate_all_the_possible_combinations_of_a_Degenerate_DNARNA_sequence_in_FASTA_format_using_Python) (I haven't vetted the script yet myself by taking what it produces and seeing if all of them match according to PatMatch.
