# hhsuite3-utilities

Utility scripts for working with command line-based [HH-suite3](https://github.com/soedinglab/hh-suite/wiki#generating-a-multiple-sequence-alignment-using-hhblits), which is already installed and working in the Jupyter instances launchable [here](https://github.com/fomightez/hhsuite3-binder).

# The scripts

* hhr_results_to_df.py
> HH-suite3 hhr results file format --> dataframe of data for use in Python

Takes output from command line BLAST and makes a dataframe from it for use with Python.

Verified compatible with both Python 2.7 and Python 3.6.

Written to run from command line or pasted/loaded inside a Jupyter notebook cell.  
The main ways to run the script are demonstrated in [a series of notebooks found here](https://github.com/fomightez/blast-binder). They can be run actively by pressing the `launch binder` button [there](https://github.com/fomightez/blast-binder). (https://nbviewer.jupyter.org/github/fomightez/sequencework/blob/master/blast-utilities/demo_blast2df_in_python2.ipynb) in Python 2.)

Assumes BLAST is run with the flag below:
```
-outfmt "6 qseqid sseqid stitle pident qcovs length mismatch gapopen qstart qend sstart send qframe sframe frames evalue bitscore qseq sseq"
```
Example BLAST command to make the results needed by `blast_to_df.py`:
```
blastn -query sequence.fa -db db_seq.fa -outfmt "6 qseqid sseqid stitle pident qcovs length mismatch gapopen qstart qend sstart send qframe sframe frames evalue bitscore qseq sseq" -out seq.x.db_seq.comp.txt
```

See [here](https://medium.com/@auguste.dutcher/turn-blast-results-into-a-presence-absence-matrix-cc44429c814) for illustration of that approach, and [here](https://blastedbio.blogspot.com/2014/11/column-headers-in-blast-tabular-and-csv.html) for a listing of explanations for those text codes.  
For drafting the actual BLAST command, I also found these helpful:

- [Running command-line BLAST](https://angus.readthedocs.io/en/2017/running-command-line-blast.html) 
- [Shawn  T. O’Neil's A Primer for Computational Biology: Chapter 7 Command Line BLAST](http://library.open.oregonstate.edu/computationalbiology/chapter/command-line-blast/)


Example calls to run the `blast_to_df.py` script from command line:
```
python blast_to_df.py results_file.txt
```

(Alternatively, upload the script to a Jupyter environment and use `%run blast_to_df.py results_file.txt` in a Python-backed notebook to run the example.)


**The pickle format ending in extension `.pkl` was chosen as the default output, despite not being human-readable, in order to store data efficiently since searches of entire gemomes has the potential to generate a lot of hits.** If you need to convert to a text-readable form, you can do the following with Jupyter or IPython where the pickled data has been saved to save it as tab separated text:

    import pandas as pd
    df = pd.read_pickle("patmatch_pickled_df.pkl")
    df.to_csv('patmatch_data.tsv', sep='\t',index = False) 



#### For running in a Jupyter notebook:

To use this script after pasting or loading into a cell in a Jupyter notebook, in the next cell define the URL and then call the main function similar to below:
```
df = blast_to_df("results_file.txt")
```
See [here](https://git.io/vh8M7) and for notebooks demonstrating use within a Jupyter notebook; click `luanch binder` to launch the notebooks from there.


Related
-------

- ?
