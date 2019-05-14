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

It probably is warning so you know that it may not show a match if it happens to match where the split was chosen. They way I dealt with this was to remove the warnings from the input to `patmatch_results_to_df.py` or `patmatch_results_to_df()` because otherwise they cause errors when the lines are processed as typical output and there are no start or end coordinates to extract. However, instead of having `patmatch_results_to_df.py` handle the issues, I thought it best the user has to deal with the warnings so they are aware.


Related
-------

[patmatch-binder](https://github.com/fomightez/patmatch-binder) - for running command line-based PatMatch on sequences of your choice in your browser without need for downloads, installations, or maintenance.
