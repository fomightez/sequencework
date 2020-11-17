# hhsuite3-utilities

Utility scripts for working with command line-based [HH-suite3](https://github.com/soedinglab/hh-suite/wiki), which is already installed and working in the Jupyter instances launchable [here](https://github.com/fomightez/hhsuite3-binder).

# The scripts

* hhsuite3_results_to_df.py
> HH-suite3 hhr results file format --> dataframe of data for use in Python

Takes output in the form of `.hhr` results files from command line HH-suite3 programs and makes a dataframe from it for use with Python.

Verified compatible with both Python 2.7 and Python 3.8.

Written to run from command line or pasted/loaded inside a Jupyter notebook cell.  
The main ways to run the script are demonstrated in [a series of notebooks found here](https://github.com/fomightez/hhsuite3-binder). They can be run actively by pressing the `launch binder` button [there](https://github.com/fomightez/hhsuite3-binder). (https://nbviewer.jupyter.org/github/fomightez/sequencework/blob/master/blast-utilities/demo_blast2df_in_python2.ipynb) in Python 2.)

Assumes BLAST is run with the flag below:
```
-outfmt "6 qseqid sseqid stitle pident qcovs length mismatch gapopen qstart qend sstart send qframe sframe frames evalue bitscore qseq sseq"
```
Example `hhblits` command to make the results needed by `hhsuite3_results_to_df.py`:
```
hhblits ...
```

See [HH-suite3 dcoumentation](https://github.com/soedinglab/hh-suite/wiki) for more on generating results files.


Example calls to run the `hhsuite3_results_to_df.py` script from command line:
```
python hhsuite3_results_to_df.py results_file.hhr
```

(Alternatively, upload the script to a Jupyter environment and use `%run hhsuite3_results_to_df.py results_file.hhr` in a Python-backed notebook to run the example.)


**The pickle format ending in extension `.pkl` was chosen as the default output, despite not being human-readable, in order to store data efficiently since searches of entire gemomes has the potential to generate a lot of data files.** If you need to convert to a text-readable form, you can do the following with Jupyter or IPython where the pickled data has been saved to save it as tab separated text:

    import pandas as pd
    df = pd.read_pickle("hhr_derived_df.pkl")
    df.to_csv('hhr_results.tsv', sep='\t',index = False) 



#### For running in a Jupyter notebook:

To use this script after pasting or loading into a cell in a Jupyter notebook, in the next cell define the URL and then call the main function similar to below:
```
df = hhsuite3_results_to_df("results_file.hhr")
```
See [here](????) and for notebooks demonstrating use within a Jupyter notebook; click `luanch binder` to launch the notebooks from there.


Related
-------

- ?
