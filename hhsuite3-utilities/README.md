# hhsuite3-utilities

Utility scripts for working with command line-based [HH-suite3](https://github.com/soedinglab/hh-suite/wiki), which is already installed and working in the Jupyter instances launchable [here](https://github.com/fomightez/hhsuite3-binder).

# The scripts

* hhsuite3_results_to_df.py
> HH-suite3 hhr results file format --> dataframe of data for use in Python

Takes output in the form of `.hhr` results files from command line HH-suite3 programs and makes a dataframe from it for use with Python.

Verified compatible with both Python 2.7 and Python 3.8.

Written to run from command line or pasted/loaded inside a Jupyter notebook cell.  
The main ways to run the script are demonstrated in [a series of notebooks found here](https://github.com/fomightez/hhsuite3-binder). They can be run actively by pressing the `launch binder` button [there](https://github.com/fomightez/hhsuite3-binder). (https://nbviewer.jupyter.org/github/fomightez/sequencework/blob/master/blast-utilities/demo_blast2df_in_python2.ipynb) in Python 2.)


Example `hhblits` command to make an HH-suite3 results file in `.hhr` format needed as input by `hhsuite3_results_to_df.py`:
```
hhblits ...
```

See [HH-suite3 dcoumentation](https://github.com/soedinglab/hh-suite/wiki) for more on generating results files on the command line.

Alternatively, you can use the HHpred webserver to make these result files using your favorite web browser. For example, go to [HHpred site here](https://toolkit.tuebingen.mpg.de/tools/hhpred) and paste in a protein sequence. Feel free to adjust the search options if you'd like before hitting the 'Submit' button in the bottom right. After the submitted job finishes, from the 'Results' tab, you can select 'Download HHR' to retrieve a HH-suite3 results file in `.hhr` format to your local machine. You can upload that to a running session launched from [here](https://github.com/fomightez/hhsuite3-binder) and edit some of the examples in the series of notebooks there to use the hhsuite3_results_to_df.py script on your data present in your downloaded `.hhr` file.


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
See [here](https://github.com/fomightez/hhsuite3-binder) for notebooks demonstrating use within a Jupyter notebook; click `launch binder` to launch a session that will allow you to us the notebooks from there.


Related
-------

- ?
