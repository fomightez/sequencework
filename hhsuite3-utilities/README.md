# hhsuite3-utilities

Utility scripts for working with command line-based [HH-suite3](https://github.com/soedinglab/hh-suite/wiki), which is already installed and working in the Jupyter instances launchable [here](https://github.com/fomightez/hhsuite3-binder).

# The scripts

* hhsuite3_results_to_df.py
> HH-suite3 hhr results file format --> dataframe of data for use in Python

Takes output in the form of `.hhr` results files from command line HH-suite3 programs and makes a dataframe from it for use with Python.

Intended to be compatible with both Python 2.7 and Python 3.8. (Not yet verified though. **UPDATE WHEN VERIFIED**)

Written to run from command line or pasted/loaded inside a Jupyter notebook cell.  
The main ways to run the script are demonstrated in [a series of notebooks found here](https://github.com/fomightez/hhsuite3-binder). They can be run actively by pressing the `launch binder` button [there](https://github.com/fomightez/hhsuite3-binder). (https://nbviewer.jupyter.org/github/fomightez/sequencework/blob/master/blast-utilities/demo_hhsuite3_results_to_df_in_python2.ipynb) in Python 2. **<===TBD**)


Example `hhblits` command to make an HH-suite3 results file in `.hhr` format needed as input by `hhsuite3_results_to_df.py`:
```
hhblits ...
```

See [HH-suite3 documentation](https://github.com/soedinglab/hh-suite/wiki) for more on generating results files on the command line.

Alternatively, you can use the HHpred webserver to make these result files using your favorite web browser. For example, go to [HHpred site here](https://toolkit.tuebingen.mpg.de/tools/hhpred) and paste in a protein sequence. Feel free to adjust the search options if you'd like before hitting the 'Submit' button in the bottom right. After the submitted job finishes, from the 'Results' tab, you can select 'Download HHR' to retrieve a HH-suite3 results file in `.hhr` format to your local machine. You can upload that to a running session launched from [here](https://github.com/fomightez/hhsuite3-binder) and edit some of the examples in the series of notebooks there to use the hhsuite3_results_to_df.py script on your data present in your downloaded `.hhr` file.


Example calls to run the `hhsuite3_results_to_df.py` script from command line:
```
python hhsuite3_results_to_df.py results_file.hhr
```

(Alternatively, upload the script to a Jupyter environment and use `%run hhsuite3_results_to_df.py results_file.hhr` in a Python-backed notebook to run the example.)

**Note:** This script keeps each hit as a separate entry in the resulting dataframe. The biopython hhsuite parser also demonstrated in the notebook 'Using Python to examine HH-suite3 results.ipynb' & LINK GOES HERE,  that demonstrates this script, makes a text table that combines 'high-scoring segment pair' or 'high-scoring pair' (hsp) that correspond to the same accession identifier. You could always use Panda's `groupby` in a downstream step to achieve a similar accounting with the dataframe made by `hhsuite3_results_to_df.py`.

**The pickle format ending in extension `.pkl` was chosen as the default output, despite not being human-readable, in order to store data efficiently since searches of entire gemomes has the potential to generate a lot of data files.** If you need to convert to a text-readable form, you can do the following with Jupyter or IPython where the pickled data has been saved to save it as tab separated text:

    import pandas as pd
    df = pd.read_pickle("hhr_derived_df.pkl")
    df.to_csv('hhr_results.tsv', sep='\t',index = False) 


**CAVEAT**: For now, this script fails with a 'key error' if the hit identifier and the indentifier on corresponding sequence line of the alignment block don't match, for example [hhsearch_q9bsu1_uniclust_w_ss_pfamA_30.hhr](https://github.com/biopython/biopython/blob/master/Tests/HHsuite/hhsearch_q9bsu1_uniclust_w_ss_pfamA_30.hhr). (Contrast that with the details and alignment blocks [here](https://github.com/biopython/biopython/blob/master/Tests/HHsuite/2uvo_hhblits.hhr), where for example `2uvo_A` is after the `>` for hit 'No 1' and is after the 'T' below 'T Consensus'.) I had thought there was a relationship when developing the script and used that as a way to extract the hit identifier ('hid') and relate lines of the alignment. Since I am unsure if this is because some of the biopython stuff is hh-suite format version 2 because it refers to that and not version 3, I am reluctant to spend much time on this now. The easiet way until I decide to change it would be to build in a pre-processing step upstream of the use `hhsuite3_results_to_df.py`( and then possibly build that same step into `hhsuite3_results_to_df.py` later) that would detect the situation and make it more like the convention seen elsewhere.  Alternatively, until I come up with handling that situation, I can just parse such results file with biopython as covered in the bottom of the same notebook where I introduce the `hhsuite3_results_to_df.py` script.


#### For running in a Jupyter notebook:

To use this script after pasting or loading into a cell in a Jupyter notebook, in the next cell define the URL and then call the main function similar to below:
```
df = hhsuite3_results_to_df("results_file.hhr")
```
See [here](https://github.com/fomightez/hhsuite3-binder) for notebooks demonstrating use within a Jupyter notebook; click `launch binder` to launch a session that will allow you to use the notebooks from there.



* make_report_on_equivalents_of_interacting_residues.py
> HH-suite3 hhr results file format using a sequence from a structure --> report on the equivalent of residues in the homolog that interact in the homologous structure. 

Takes output in the form of `.hhr` results files from command line HH-suite3 programs, along with information on the structure of the sequence used in the  and makes a dataframe from it for use with Python. Requires using in conjuction with the notebook, 'Report if residues interacting with a specific chain have equivalent residues in an hhsuite-generated alignment.ipynb', that is meant to be used when you are studing a complex that contains at least two or more proteins because it is meant to generate a report on the residues of a protein chain of interest that interacts with another in a structure and tell you what parts of the chain of interest have equivalents in a remote homolog.

Intended to be compatible with both Python 2.7 and Python 3.8. (Not yet verified though. **UPDATE WHEN VERIFIED**)

Written to run from a Jupyter notebook cell where extra information needed has been defined earlier in the calling notebook.  
The way to run this script is demonstrated in [a series of notebooks found here](https://github.com/fomightez/hhsuite3-binder). 

Example:
```
%run -i make_report_on_equivalents_of_interacting_residues.py
```


* core_reduction_dataframe_steps_for_ECIIIL_genome.py
> input data derived from .hhr files > table about size reduction for orthologs

Takes output in the form of `.hhr` results files from command line HH-suite3 programs and table about size reduction for orthologs.  
Meant to be used in particular notebooks where already defined the other variables and input. See the first few lines of the script for pointers to those notebooks.

This one is hard coded for dealing with 'ECIIIL genome' for locus column; however, that coloumn could be dropped or the handling altered for other, future versions of this utility.



Related
-------

- [ys2ec_snakefiles - snakefiles to process yeast to Ec homologs](https://github.com/fomightez/ys2ec_snakefiles)https://github.com/fomightez/ys2ec_snakefiles
