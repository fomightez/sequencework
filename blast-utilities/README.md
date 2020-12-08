# BLAST-utilities

Utility scripts for working with command line-based BLAST, which is already installed and working in the Jupyter instances launchable [here](https://github.com/fomightez/blast-binder).

# The scripts

* blast_to_df.py
> BLAST results --> dataframe of data for use in Python

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
See [here](https://git.io/vh8M7) for notebooks demonstrating use within a Jupyter notebook; click `launch binder` to launch a session that will allow you to use the notebooks from there.


Related
-------

- About a month before I made my `blast_to_df.py`, Noam Harel in Adi Stern's Labarotary made one shown [here](https://github.com/taliaku/SternLab/commit/68abf4f05bd3d44923e365e7cdb6ca2ecb1d19ca) that uses the hit table in CSV form as input. They have a script entitled `blast_utilities.py` [here](https://github.com/taliaku/SternLab/blob/master/blast_utilities.py).

- However, I probably lifted the name / idea from something I vaguely recalled from Titus Brown's Lamprey work, see [here](https://github.com/dib-lab/2013-lamprey/blob/7493710399a05433989e25b3d01b962e6cab3553/notebooks/analyses/Lamprey_E_Protein_Analysis.ipynb). Unfortunately the reference to `peasoup` which imports the `BLAST` methods seems defunct [here](https://github.com/dib-lab/2013-lamprey/tree/7493710399a05433989e25b3d01b962e6cab3553) and so I haven't entirely tracked it back to put a link to the code here yet. (Maybe [Camille Scott's stuff?](https://github.com/dib-lab/2013-lamprey/blob/7493710399a05433989e25b3d01b962e6cab3553/notebooks/analyses/petmar-gtf-overlap.ipynb)?)
