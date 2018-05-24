# BLAST-utilities

Utility scripts for working with command line-based BLAST, which is already installed and working in the Jupyter instances launchable [here](https://github.com/fomightez/qgrid-notebooks).

# The scripts

* BLAST_to_df.py
> BLAST results --> dataframe of data for use in Python

Takes output from command line BLAST and makes a dataframe from it for use with Python.

Verified compatible with both Python 2.7(??) and Python 3.6.

Written to run from command line or pasted/loaded inside a Jupyter notebook cell.  
The main ways to run the script are demonstrated in [a series of notebooks???](???).

Assumes BLAST is run with the flag below:
```
-outfmt "6 qseqid sseqid stitle pident qcovs length mismatch gapopen qstart qend sstart send qframe sframe frames evalue bitscore qseq sseq"
```

See [here](https://medium.com/@auguste.dutcher/turn-blast-results-into-a-presence-absence-matrix-cc44429c814) for illustration of that approach, and [here](https://blastedbio.blogspot.com/2014/11/column-headers-in-blast-tabular-and-csv.html) for a listing of explanations for those text codes.  
For drafting the actual BLAST command, I also found these helpful:

- [Running command-line BLAST](https://angus.readthedocs.io/en/2017/running-command-line-blast.html) 
- [Shawn  T. O’Neil's A Primer for Computational Biology: Chapter 7 Command Line BLAST](http://library.open.oregonstate.edu/computationalbiology/chapter/command-line-blast/)


Example calls to run script from command line:
```
python BLAST_to_df.py results_file.txt
```

(Alternatively, upload the script to a Jupyter environment and use `%run BLAST_to_df.py results_file.txt` in a Python-backed notebook to run the example.)




#### For running in a Jupyter notebook:

To use this script after pasting or loading into a cell in a Jupyter notebook, in the next cell define the URL and then call the main function similar to below:
```
df = BLAST_to_df("results_file.txt")
```
See [here????](???????) for notebooks demonstrating use within a Jupyter notebook.


Related
-------

???
