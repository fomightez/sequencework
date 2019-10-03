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



**The pickle format ending in extension `.pkl` was chosen as the default output, despite not being human-readable, in order to store data efficiently since searches of entire gemomes has the potential to generate a lot of hits.** If you need to convert to a text-readable form, you can do the following with Jupyter or IPython where the pickled data has been saved to save it as tab separated text:

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

Takes pattern in PatMatch syntax and a sequence in text or FAST format and reports if they match.

There is a [demo notebook for this script that runs in my patmatch-binder](https://github.com/fomightez/patmatch-binder). You can launch the series from [here](https://github.com/fomightez/patmatch-binder) and then selecting from the list at the bottom of the index page undert 'Additional topics: Technicall' to go to the 'Demo of script to check for a match to a sequence pattern in PatMatch syntax' page. The direct link to a nicely-rendered, static version of that page is [here](https://nbviewer.jupyter.org/github/fomightez/patmatch-binder/blob/master/notebooks/Demo%20of%20script%20to%20check%20for%20a%20match%20to%20a%20sequence%20pattern%20in%20PatMatch%20syntax.ipynb).


Related
-------

[patmatch-binder](https://github.com/fomightez/patmatch-binder) - for running command line-based PatMatch on sequences of your choice in your browser without need for downloads, installations, or maintenance.
