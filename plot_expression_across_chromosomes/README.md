# plot_expression_across_chromosomes.py

> just a draft of documentation for now.  
> Contact me if you need help until I better document this.

A python script to olot ratio of expression of experimental condition vs. wild-type (or 
baseline state) for genes in sequential order across chromosomes in the genome. It requires two files: 1. a 
file of a genome annotation format in order to parse the locations of genes 
and (approximate) length of chromosomes; 2. a file of expression data to plot.
Currently, the script necessitates that the data is in a tab-delimited format. 
At this point to accomodate replicates, it is assumed you have in your data file for each gene the 
resulting (combined) level metric for the replicates, such as the mean TPM for
your "wild-type" samples and your mean TPM for your experimental samples. In 
the future there will be a related script for "raw" results produced by Salmon 
or HTSeq. The hope being you just need to point the script at the raw data 
files and it will automagically handle the combining and produce a plot 
showing the expression of genes across the chromosomes.


The plotting approach and other aspects borrow from Brent Pedersen's awesome
`manhattan-plot.py` script [here](https://github.com/brentp/bio-playground/blob/master/plots/manhattan-plot.py). (Many thanks to him for sharing that!)
This script should produce a plot similar to a combination of Brent Pedersen's
`manhattan-plot`and Figure 5B from Thorburn et al. 2013 (PMID: 23468524), (see [here](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3639041/)).

There are several optional flags that can be supplied at the time of calling
the script to control options. These are shown if you invoke with `-help` 
flag, i.e., type `python plot_expression_across_chromosomes.py --help` or simply call the script with no additional arguments. Additionally, 
inside the script there are several `USER ADJUSTABLE VALUES` that can be 
edited for easy customization.
Built to be general enough to be easily modified. For example, should be 
fairly straightforward to change to plotting ratio of sequencing coverage 
across chromosomes in a genome.

To avoid compressing the typically important range, by default the y-axis is limited to a range that emphasizes the alterations characteristic of aneuploidy. In that case, points beyond those limits are plotted at the edge of the plot in a manner that distinguishes them from the other data points, i.e., they are represented as open triangles at the edge pointing in an up-or-down direction. This approach is styled on how DESeq2 plotMA handles out of bounds points, see section 1.5.1 of [the DESeq2 vignette](https://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#ma-plot). Depnding on your needs, run with the optional `--no_limits` or `--no_log` flags to see the full plot.

