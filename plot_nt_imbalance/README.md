Scripts for plotting nucleotide imbalance from DNA sequences
============================================================

# The scripts

* nucleotide_difference_imbalance_plot_stylized_like_Figure_8_of_Morrill_et_al_2016.py
> FASTA sequence for a region of a chromosome/contig/scaffold --> plot of GvsC and AvsT imbalance

See [Figure 8, panel B](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4882425/figure/F8/) of [Morrill et al 2016(PMID: 27026700)](https://www.ncbi.nlm.nih.gov/pubmed/27026700).

There is a demonstration of this script available in my set of sequence analysis related demonstrations notebooks [here](https://github.com/fomightez/cl_sq_demo-binder). To run it actively,  launch a binder session by clicking on the `launch binder` badges [here](https://github.com/fomightez/cl_sq_demo-binder) and select the link to 'Demo of script to plot nt imbalance for sequence span' from the list of notebooks to run it actively.  The particular notebook can be viewed statically, nicely displayed [here](https://nbviewer.jupyter.org/github/fomightez/cl_sq_demo-binder/blob/master/notebooks/Demo%20of%20script%20to%20plot%20nt%20imbalance%20for%20sequence%20span.ipynb).

* two_nucleotides_in_proximity_difference_imbalance_plot.py
> FASTA sequence for a region of a chromosome/contig/scaffold --> plot of imbalance of the sum of two bases vs the sum of the other two

This is like `nucleotide_difference_imbalance_plot_stylized_like_Figure_8_of_Morrill_et_al_2016.py`; however, it allows specifying a combination of any two basepairs to see of the amount of those in proximity is out of balance vs the sum of the other two basepairs. For example, `GandC` vs. `AandT`. 

Using the script is similar to `nucleotide_difference_imbalance_plot_stylized_like_Figure_8_of_Morrill_et_al_2016.py` and so look over the demo notebook about that and then to use `two_nucleotides_in_proximity_difference_imbalance_plot.py` add an argument of two nucleotides to group, such as `GC` when calling it on the command line. Or specify `dibase_text1`, such as `dibase_text1="GC"`, when calling the main function of the script when in a Jupyter notebook cell or in IPython.  
Example invocations of the script or main function:

- from command line:

```bash
python two_nucleotides_in_proximity_difference_imbalance_plot sequence.fa 20000 GC
```

- from importing main function when in Jupyter or IPython:

```python
%matplotlib inline
from two_nucleotides_in_proximity_difference_imbalance_plot import two_nucleotides_in_proximity_difference_imbalance_plot
two_nucleotides_in_proximity_difference_imbalance_plot("sequence.fa", 20000, dibase_text1="GC", return_plot=True);
```

Related items created by myself
-------------------------------

?

Related items by others
-----------------------

?
