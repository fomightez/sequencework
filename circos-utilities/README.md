# circos-utilities

Utility scripts for working with Circos.

# The scripts

* UCSC_chrom_sizes_2_circos_karyotype.py
> UCSC chrom.sizes files --> karyotype.tab file for use in Circos

Takes a URL for a UCSC `chrom.sizes` file and makes a `karyotype.tab` file from it for use with Circos.

Verified compatible with both Python 2.7 and Python 3.6.

Written to run from command line or pasted/loaded/imported inside a Jupyter notebook cell.  
The main ways to run the script are demonstrated in the notebook [`demo UCSC_chrom_sizes_2_circos_karyotype script.ipynb`](https://github.com/fomightez/sequencework/blob/master/circos-utilities/demo%20UCSC_chrom_sizes_2_circos_karyotype%20script.ipynb) that is included in this repository. (That notebook can be viewed in a nicer rendering [here](https://nbviewer.jupyter.org/github/fomightez/sequencework/blob/master/circos-utilities/demo%20UCSC_chrom_sizes_2_circos_karyotype%20script.ipynb).)

To determine the URL to feed the script, google `YOUR_ORGANISM genome UCSC chrom.sizes`,  where you replace `YOUR_ORGANISM` with your organism name and then adapt the path you see in the best match to be something similar to 
`http://hgdownload.cse.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.chrom.sizes` -or-
`http://hgdownload.cse.ucsc.edu/goldenPath/canFam2/bigZips/canFam2.chrom.sizes`.  
You can get an idea of what is available by exploring the top section [here](http://hgdownload.cse.ucsc.edu/downloads.html); clicking on the arrows creates drop-down lists reveal many genomes for each category.

Importantly, this script is intended for organisms without cytogenetic bands, such as dog, cow, yeast, etc..  
(For organisms with cytogenetic band data: Acquiring the cytogenetic bands information is described [here](http://circos.ca/tutorials/lessons/ideograms/karyotypes/), about halfway down 
the page where it says, "obtain the karyotype structure from...". 
Unfortunately, it seems the output to which one is directed to by those instructions is not
directly useful in Circos(?). Fortunately, though as described [here](http://circos.ca/documentation/tutorials/quick_start/hello_world/), "Circos ships with several predefined karyotype files for common sequence 
assemblies: human, mouse, rat, and drosophila. These files are located in 
data/karyotype within the Circos distribution." And also included there is a script for converting the cytogenetic band data to karyotype, see [here](http://circos.ca/documentation/tutorials/quick_start/hello_world/) and [here](https://groups.google.com/d/msg/circos-data-visualization/B55NlByQ6jY/nKWGSPsXCwAJ).)

Example call to run script from command line:
```
python UCSC_chrom_sizes_2_circos_karyotype.py http://hgdownload.cse.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.chrom.sizes
```
(Alternatively, upload the script to a Jupyter environment and use `%run UCSC_chrom_sizes_2_circos_karyotype.py http://hgdownload.cse.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.chrom.sizes` in a Python-backed notebook to run the example.)

Example input from http://hgdownload.cse.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.chrom.sizes :
```
chrIV   1531933
chrXV   1091291
chrVII  1090940
chrXII  1078177
chrXVI  948066
chrXIII 924431
chrII   813184
chrXIV  784333
chrX    745751
chrXI   666816
chrV    576874
chrVIII 562643
chrIX   439888
chrIII  316620
chrVI   270161
chrI    230218
chrM    85779
```

Example output sent to file (tab-separated):
```
chr -   Sc-chrIV    chrIV   0   1531933 black
chr -   Sc-chrXV    chrXV   0   1091291 black
chr -   Sc-chrVII   chrVII  0   1090940 black
chr -   Sc-chrXII   chrXII  0   1078177 black
chr -   Sc-chrXVI   chrXVI  0   948066  black
chr -   Sc-chrXIII  chrXIII 0   924431  black
chr -   Sc-chrII    chrII   0   813184  black
chr -   Sc-chrXIV   chrXIV  0   784333  black
chr -   Sc-chrX chrX    0   745751  black
chr -   Sc-chrXI    chrXI   0   666816  black
chr -   Sc-chrV chrV    0   576874  black
chr -   Sc-chrVIII  chrVIII 0   562643  black
chr -   Sc-chrIX    chrIX   0   439888  black
chr -   Sc-chrIII   chrIII  0   316620  black
chr -   Sc-chrVI    chrVI   0   270161  black
chr -   Sc-chrI chrI    0   230218  black
chr -   Sc-chrM chrM    0   85779   black
```


#### For running in a Jupyter notebook:

To use this script after pasting or loading into a cell in a Jupyter notebook, in the next cell define the URL and then call the main function similar to below:
```
url = "http://hgdownload.cse.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.chrom.sizes"
species_code = "Ys"
UCSC_chrom_sizes_2_circos_karyotype(url, species_code)
```
-or-

```
UCSC_chrom_sizes_2_circos_karyotype(url)
```
Without supplying a second argument, a species code will be extracted automatically and used.

Note that `url` is actually not needed if you are using the yeast one because that specific one is hardcoded in script as default.
In fact, because I hardcoded in defaults, just `main()` will indeed work for yeast after script pasted in or loaded into a cell.
See [here](https://nbviewer.jupyter.org/github/fomightez/sequencework/blob/master/circos-utilities/demo%20UCSC_chrom_sizes_2_circos_karyotype%20script.ipynb) for a notebook demonstrating use within a Jupyter notebook.


Related
-------

- [circos-binder](https://github.com/fomightez/circos-binder) - for running Circos in your browser without need for downloads, installations, or maintenance.

- ### Related 

- [gos: (epi)genomic visualization in python](https://gosling-lang.github.io/gos/) looks to circos-like images, see that page in the link for representative examples in the image at the top. [The Example Gallery](https://gosling-lang.github.io/gos/gallery/index.html) has a link to a Circos-style example under the 'Others' heading; it includes code [here](https://gosling-lang.github.io/gos/gallery/circos.html).
