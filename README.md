sequencework
============

Mainly python scripts related to nucleic or protein sequence work

I sorely need to put an index here with links to better guide to the appropriate folders. <== **TO DO**  
(For now look at the title of the folders to try and discern if it is something of interest.)

Descriptions of the scripts are found within `README.md` files in the sub folders.

Several have demonstrations in sessions served by MyBinder.org [from my command line-sequence associated repo](https://github.com/fomightez/cl_sq_demo-binder); however, probably best to follow guide listed with individual scripts so that you quickly find the right location. If you already know where you are going, you can launch a session via this button:

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/fomightez/cl_sq_demo-binder/master?filepath=index.ipynb)

Related gist by me
------------------

- [snippets for dealing with FASTA - `useful_FASTA_handling.py`](https://gist.github.com/fomightez/8cd6d9ba88f975b64e43eba562894dec) includes splitting

Related 'Binderized' Utilities
----------------------------

Collection of links to launchable Jupyter environments where various sequence analysis tools work **WITHOUT ANY NEED FOR ADDITIONAL EFFORT/INSTALLS**. Many of my recent scripts are built with use in these environments in mind:

(Many of these include/feature Biopython, too, such as  but I haven't made a one all encompassing one yet for that since I use it a lot as an underlying library.)

- [patmatch-binder](https://github.com/fomightez/patmatch-binder) - launchable Jupyter sessions for running command line-based PatMatch in Jupyter environment provided via Binder (*Perl and Python-based*).

- [blast-binder](https://github.com/fomightez/blast-binder) - launchable Jupyter sessions for running command line-based BLAST+ in Jupyter environment provided via Binder.

- [Demonstration Jupyter Notebooks for My imporved version ofAdam Bessa's Fasta2Structure - Fasta2Structure-cli](https://github.com/fomightez/Fasta2Structure-cli) - To make it more convenient to use, I've modified the Fasta2Structure script to allow more ways to run it to produce improved_Fasta2Structure (a.k.a Fasta2Structure-cli).  It will run on the command line if you supply arguments specifying files as input or fallback to running on the command line if Tkinter cannot connect to a graphical display.  improved_Fasta2Structure (a.k.a Fasta2Structure-cli) - User-Friendly Tool for Converting Multiple Aligned FASTA Files to STRUCTURE Format, that is even more user-friendly because it doesn't need a user to select files in a GUI (Tkinter-based) and can thus run well anywhere, such as on a computer cluster or in Jupyter running remotely or in conjunction with software to make pipelines like Snakemake & NextFlow. For those reasons, the improved script is more user-friendly for those familiar with computation and allows scaling up.

- [InterMine-binder](https://github.com/fomightez/intermine-binder) - 
Intermine Web Services available in a Jupyter environment running via the Binder service. (See the [ guide to getting started with using Intermine sites and Jupyter using MyBinder-served Jupyter notebooks](https://github.com/fomightez/guide_to_intermine-binder).)

- [mcscan-binder](https://github.com/fomightez/mcscan-binder) - 
MCscan software available in a launchable Jupyter environment running via the Binder service (*Python 2-based*), with an example workflow and some other use examples.

- [mcscan-blast-binder](https://github.com/fomightez/mcscan-blast-binder) - MCscan and BLAST+ command line software available in a launchable Jupyter environment running via the Binder service (*Python 2-based*).

- [synchro-binder](https://github.com/fomightez/synchro-binder) - SynChro software available in a launchable Jupyter environment running via the Binder service with Quick start and some other illustrations of its use.

- [cl_sq_demo-binder](https://github.com/fomightez/cl_sq_demo-binder) - launchable, working Jupyter-based environment that has a collection of demonstrations of useful resources on command line (or useable in Jupyter sessions) for manipulating sequence files. (Note: THIS WAS STARTED AFTER SEVERAL OTHER DEMO NOTEBOOKS (many meant to be static) MADE FOR SEQUENCE SCRIPTs, and hopefully slowly those will be added to here as well to be available in active form.)

- [clausen_ribonucleotides binder](https://github.com/fomightez/clausen_ribonucleotides) - Analyze ribonucleotide incorporation data from Clausen et al. 2015 data using script `plot_5prime_end_from_bedgraph.py`.


- [circos-binder](https://github.com/fomightez/circos-binder) - Circos software available in a launchable Jupyter environment running via the Binder service with tutorials illustrating use (TBD)(*Perl and Python-based*).

Related resources by others
---------------------------

- [genomepy](https://github.com/vanheeringen-lab/genomepy)
>"Install and use genomes & gene annotations the easy way!  
genomepy is designed to provide a simple and straightforward way to download and use genomic data. This includes (1) searching available data, (2) showing the available metadata, (3) automatically downloading, preprocessing and matching data and (4) generating optional aligner indexes. All with sensible, yet controllable defaults. Currently, genomepy supports UCSC, Ensembl and NCBI." - Includes an S. cerevisiae example.
- [rna-tools (previously rna-pdb-tools): a toolbox to analyze sequences, structures and simulations of RNA](https://github.com/mmagnus/rna-tools/blob/master/index-of-tools.md)
- ['seqrequester', a tool for summarizing, extracting, generating and modifying DNA sequences.](https://github.com/marbl/seqrequester). Perl-based.
- [SeqKit - a cross-platform and ultrafast toolkit for FASTA/Q file manipulation](https://github.com/shenwei356/seqkit) by shenwei356 has some excellent utilities for handling FASTA or FASTQ files. See [here](https://github.com/shenwei356/seqkit#subcommands) for subcommands listing. [Documentation](https://bioinf.shenwei.me/seqkit/).
- [SeqFu](https://telatin.github.io/seqfu2/) - an easy to use toolkit for FASTA and FASTQ manipulation and inspection on the commandline. Available from Bioconda as ‘seqfu’.
    >"A general-purpose program to manipulate and parse information from FASTA/FASTQ files, supporting gzipped input files. Includes functions to interleave and de-interleave FASTQ files, to rename sequences and to count and print statistics on sequence lengths. SeqFu is available for Linux and MacOS. - A compiled program delivering high performance analyses - Supports FASTA/FASTQ files, also Gzip compressed    - A growing collection of handy utilities, also for quick inspection of the datasets."
        - [Example uses Biopython to make a Pandas dataframe from FASTA sequences](https://www.biostars.org/p/9505933/#9505933)
- https://bsky.app/profile/robert.bio/post/3lbxcdesrwc2t
  >"I made a little swiss army knife tool for quickly inspecting FASTA/FASTQ sequences, including GC content, codon analysis, reverse complement, and file format conversions:  
 https://42basepairs.com/tools/sequence-analysis"  
  >""I also have calculators for exploring:  
  ➡️ FASTQ base quality: https://42basepairs.com/tools/fastq-base-quality  
  ➡️ SAM flags: https://42basepairs.com/tools/sam-flag "  

  
See also
-------

[My simulated data repo](https://github.com/fomightez/simulated_data) has some useful scripts and resources for generating simulated (mock / fake) sequence data, gene expression data, or gene lists.

[My structurework repo](https://github.com/fomightez/structurework/) - for utilities and code dealing with molecular structures

[My proteomicswork repo](https://github.com/fomightez/proteomicswork/) - for utilities and code dealing with proteomics analysis




