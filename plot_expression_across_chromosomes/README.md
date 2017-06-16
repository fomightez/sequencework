# plot_expression_across_chromosomes.py

> plots a ratio of expression values across chromosomes or scaffolds of a genome to highlight regions of deviation.

A python script to plot ratio of expression of experimental condition vs. wild-type (or 
baseline state) for genes in sequential order across chromosomes in the genome. It requires two files: 1. a 
file of a genome annotation format in order to parse the locations of genes 
and (approximate) length of chromosomes; 2. a file of expression data to plot.
Currently, the script necessitates that the data is in a tab-delimited format. 
At this point to accomodate replicates, it is assumed you have in your data file for each gene the 
resulting (combined) level metric for the replicates, such as the mean TPM for
your "wild-type" samples and your mean TPM for your experimental samples. In 
the future there will be a related script for "raw" results produced by Salmon 
or HTSeq-count. The hope being you just need to point the script at the raw data 
files and it will automagically handle the combining and produce a plot 
showing the expression of genes across the chromosomes.

Importantly, the plot script provided here is meant to be pipeline-agnostic. This means you should be able to take output from analysis of almost any RNA-Seq and generate such a plot. Every effort was made to make it not depend on what upstream or downstream analyses you used. See more about this under details. 

You can generate chromosome profiles for individual chromosomes or only a few chromosomes or scaffolds by using the `--chrs` flag to limit the analysis to specific chromosomes or scaffolds.

## USAGE SUMMARY

```
usage: plot_expression_across_chromosomes.py [-h] [-cols COLUMNS] [-l]
                                             [-chr CHRS] [-nl] [-nlim] [-s]
                                             [-ed EXP_DESIG] [-bd BASE_DESIG]
                                             [-ndh] [-ac ADVANCE_COLOR]
                                             ANNOTATION_FILE DATA_FILE

plot_expression_across_chromosomes.py plots a ratio of expression values
across chromosomes or scaffolds of a genome to highlight regions of deviation.
Besides the options listed here, there are several `USER ADJUSTABLE VALUES`
inside the script that can be edited for easy customization. A similar plot is
called a Manhattan plot and this implementation borrows the plotting approach
and some of the features from Brent Pedersen's awesome `manhattan-plot.py`
script. **** Script by Wayne Decatur (fomightez @ github) ***

positional arguments:
  ANNOTATION_FILE       Name of file containing the genome annotation.
                        REQUIRED. This is needed to determine the order of
                        individual data points along the chromosome and how to
                        display the data across chromosomes or scaffolds.
  DATA_FILE             Name of file containing the summarized data to plot,
                        such as mean TPM or RPKM, etc. in tab-delimited form.
                        REQUIRED. See my script
                        `plot_expression_across_chromosomes_from_raw.py` if
                        you want supply the individual `raw` data files with
                        the level metric for each sample &/or replicate.

optional arguments:
  -h, --help            show this help message and exit
  -cols COLUMNS, --columns COLUMNS
                        columns for gene, wild-type (baseline state)
                        expression value, experimental condition expression
                        value, in that order. This flag is used to specify the
                        data in the summary file to be plotted. Separate the
                        column identifiers by commas, without spaces. Default
                        is `1,2,3`, where `1` indicates the first column,
                        i.e., how you'd refer to the columns in natural
                        language (no zero-indexing).
  -l, --lines           add this flag to plot the expression level ratio value
                        as lines extending from the x-axis rather than points
                        in space. (The resulting aesthetic may resemble a city
                        skyline for which the `manhattan plot` is named.)
  -chr CHRS, --chrs CHRS
                        use this flag to limit plotting of the data to
                        particular chromosomes or scaffolds you specify
                        immediately following this flag. Separate the
                        chromosome or scaffold identifiers by commas, without
                        spaces. Default when this optional flag is not called
                        is to plot that data for all chromosomes or scaffolds.
  -nl, --no_log         add this flag to keep the expression level ratio to be
                        plotted in the common base 10 instead of converting to
                        log2.
  -nlim, --no_limits    add this flag to not impose a limit of above and below
                        4 in plot window when converting to log2. The cutoff
                        can also be adjusted under `user-adjustable settings`
                        in the script. Issuing this flag has no effect if all
                        values are within +/- the cutoff interval or
                        `--no_log` is used.
  -s, --smooth          add this flag to display a smoothing curve fit to the
                        data points (LOWESS) on a per chromosome basis. This
                        option can enhance visualization of deviations
                        characteristic of aneuploidy and copy number variation
                        across the genome, both within and between
                        chromosomes.
  -ed EXP_DESIG, --exp_desig EXP_DESIG
                        Allows changing the text used in y-axis label to
                        reference experimental sample. Following `--exp_desig`
                        type what you'd like to read there instead of
                        `experimental`.
  -bd BASE_DESIG, --base_desig BASE_DESIG
                        Allows changing the text used in y-axis label to
                        reference wild-type or baseline sample. Following
                        `--base_desig` type what you'd like to read there
                        instead of `wild-type`.
  -ndh, --no_data_header
                        add this flag if there is no data header or no first
                        line of column names in the data file. Otherwise, it
                        is assumed there is and any item read as the first
                        gene identifier from the first line won't be
                        highlighted as missing from annotation. IMPORTANTLY,
                        this only affects feedback provided as script is run.
                        If the first line resembles data, i.e., numbers in
                        specified columns, it will be automagically parsed as
                        if data. Remove the header or column labels line from
                        your summary data file on the off-chance this causes
                        issues in your resulting plot.
  -ac ADVANCE_COLOR, --advance_color ADVANCE_COLOR
                        **FOR ADVANCED USE.*** Allows for advancing the color
                        selection iterator the specified number of times. The
                        idea is it allows the ability to control the color of
                        the chromosome when specifying a chromosome or
                        scaffolds to plot so you could make the color match
                        the one used when all chromsome plotted if needed.
                        Supply the number to advance after the flag on the
                        command line. For example, `-ac 4`.
```


## EXAMPLE COMMANDS AND PLOTS PRODUCED

**Command:**

  	python plot_expression_across_chromosomes.py part_of_sacena_for_test.gtf test_data.tsv
	
Explanation: minimum example plotting points along chromosomes for genome

**Command:**

  	python plot_expression_across_chromosomes.py part_of_sacena_for_test.gtf test_data.tsv --columns 1,2,4
	
Explanation: specify columns in data to use other than default.

**Command:**

  	python plot_expression_across_chromosomes.py part_of_sacena_for_test.gtf test_data.tsv --columns 1,2,4 --lines
	
Explanation: specify columns other than default and plot with vertical lines from zero.

	python plot_expression_across_chromosomes.py part_of_sacena_for_test.gtf test_data.tsv --columns 1,2,4 --no_log
	
specify columns other than default and don't convert ratio to log2 value.

	python plot_expression_across_chromosomes.py part_of_sacena_for_test.gtf test_data.tsv --columns 1,2,4 --smooth

specify columns other than default and add a line indicating a smooth curve fit to the data points

	python plot_expression_across_chromosomes.py part_of_sacena_for_test.gtf test_data.tsv --columns 1,2,4 --chrs IV,VIII,XII 
specify columns other than default and only plot for chromsomes IV, VIII, and XII. IMPORTANTLY, note no spaces between chromosome identifiers.



## DETAILS

The implemented plotting approach and other aspects borrow from Brent Pedersen's awesome
`manhattan-plot.py` script [here](https://github.com/brentp/bio-playground/blob/master/plots/manhattan-plot.py). (Many thanks to him for sharing that!)
This script should produce a plot similar to a combination of Brent Pedersen's
`manhattan-plot` that can be seen [here](https://github.com/brentp/bio-playground/tree/master/plots) and the plot in Figure 5B from Thorburn et al. 2013 (PMID: 23468524) that can be seen [here](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3639041/).

There are several optional flags that can be supplied at the time of calling
the script to control options. These are shown if you invoke with `help` 
flag, i.e., type `python plot_expression_across_chromosomes.py --help` or simply call the script with no additional arguments. Additionally, 
inside the script there are several `USER ADJUSTABLE VALUES` that can be 
edited for easy customization.
Built to be general enough to be easily modified. For example, should be 
fairly straightforward to change to plotting ratio of sequencing coverage 
across chromosomes in a genome.

To avoid compressing the typically important range, by default the y-axis is limited to a range that emphasizes the alterations characteristic of aneuploidy or segmental duplication. In that case, points beyond those limits are plotted at the edge of the plot in a manner that distinguishes them from the other data points, i.e., they are represented as open triangles at the edge pointing in an up-or-down direction. This approach is styled on how DESeq2 plotMA handles out of bounds points, see section 1.5.1 of [the DESeq2 vignette](https://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#ma-plot). Depending on your needs, run with the optional `--no_limits` or `--no_log` flags to see the full plot.

Importantly, the plot script provided here is meant to be pipeline-agnostic. By way of contrast, DESeq2 can be used to plot fold changes in genomic space; however, it requires reading in data with `summarizeOverlaps`, as described [here](https://www.bioconductor.org/help/workflows/rnaseqGene/#plotting-fold-changes-in-genomic-space). Reading in data with `summarizeOverlaps` is just one of many of the routes to supplying count matrices into DESeq2. Additionally, you need to define gene models from GTF files by also constructing a TxDb object in DESeq2 prior to even being able to use `summarizeOverlaps` to generate count matrices. The plot_expression_across_chromosomes.py script provide here obviates many of those preparatory steps and makes it possible to take any route to generating such a plot. You can still use DESeq2 if you choose but with this script you aren't tied to using `summarizeOverlaps` and have the option to use any method, say tximport, to bring in counts. Similarly, NOISeq produces Manhattan plots for chromosomes as shown in [Figure 3]((http://bioinfo.cipf.es/noiseq/lib/exe/fetch.php?media=noiseqbio_techreport.pdf), but then you are tied to NOISeq or have to do the analysis with at leat two different packages to get such a plot. Again, this script is meant to eliminate that but still allow producing such plots.

The highest value for `end` found in the annotation file (typically the `gtf` or `gff` format file) for a feature residing on a chromosome is used to determine the approximate length of each chromosome. You can add a "dummy" entry for a non-existent gene to your annotation file if you find the features in the file don't adequately reflect length of the chromosome.

You can generate chromosome profiles for individual chromosomes or only a few chromosomes or scaffolds by using the `--chrs` flag to limit the analysis to specific chromosomes or scaffolds.

The non-parametric strategy used to fit the curve to the data points is LOWESS (Locally Weighted Scatterplot Smoothing). It highlights deviations in the scatterplots better than one could estimate just by eye given so many points.

Related
-------

I have made a script that will run the plot across the entire yeast genome as well as making plots for each chromosome individually, except mitochondrial one, automagically. It could probably by easily modified for any genome and set of data where `plot_expression_across_chromosomes.py` works. There is a separate script that combines the images from several samples into a report summary for each. The two related scripts are:
* [sheperds_chr_thru_plot_expression_across_chromosomes.py](https://github.com/fomightez/mini-pipelines)
* [generate_reports_for_genome_and_all_chromosomes_various_samples.py](https://github.com/fomightez/mini-pipelines)  

They can both be found [here](https://github.com/fomightez/mini-pipelines) in my [mini-pipelines repository](https://github.com/fomightez/mini-pipelines).

Additionally, I have made a script for generating simulated data files where you can designate what chromosomes or scaffolds differ from the baseline and by what fold. That script, called `mock_expression_ratio_generator.py`, can be found [here](https://github.com/fomightez/simulated_data/tree/master/gene_expression).
