Adjust Annotation Files
=======================

**My collection of files to deal with adjusting sequence annotation files.**

---
**Description of each script**

- fix_lsu_rRNA_annotation_in_gff_resulting_from_mfannot.py

>annotation file (gff3) from fungal mitochondria without annotation for 21S rRNA --> annotation file with 21S rRNA

MFannot seems to do poorly with annotating the large ribosomal subunit, even if all the strains are cerevisiae it is hit or miss for some reason, `fix_lsu_rRNA_annotation_in_gff_resulting_from_mfannot.py` fixes that for cereivisae and other budding yeast(?).
Besides the annoation file in gff3 format produced by MFAnnot followed by conversion to GFF3 format (via [`mfannot2gff3.pl`](https://github.com/yjx1217/LRSDAY/blob/master/scripts/mfannot2gff3.pl) or [`mfannot2gff.pl`](https://github.com/kbseah/mitonotate/blob/master/mfannot2gff.pl)) it needs the corresponding sequence file that MFannot used so it can locate the best match to the coding region of S. cereivisiae 21S rRNA. It also will check for the presence of the Omega intron in the course of locating the 21S rRNA candidate.

There isn't a demo notebook for this script yet, but it is described being used as a step in a workflow/pipeline/ad-hoc series of processing in a [demo notebook for my `measure_intergenic_regions_in_mito_annotations.py` script in this repo](https://github.com/fomightez/cl_sq_demo-binder). Launch a binder session from there and select 'Demo of script to get intergenic gap sizes from annotation file' to run it actively.  The particular notebook can be viewed statically, nicely displayed [here](https://nbviewer.jupyter.org/github/fomightez/cl_sq_demo-binder/blob/master/notebooks/Demo%20of%20script%20to%20get%20intergenic%20gap%20sizes%20from%20annotation%20file.ipynb).

There are some related scripts in the ['omega-presence' subfolder in my sequencework repo](https://github.com/fomightez/sequencework/tree/master/omega-presence).

- makes_length_annotation_file.py

> annotation file --> differently-arranged annotation file  
`makes_length_annotation_file.py` takes a file of `biomart_length.txt` obtained using [Biomart tool to access Ensemble data](http://useast.ensembl.org/info/data/biomart/index.html) and returns a file with gene_id
and then gene length for use with NOISeq. See the User's Guide for NOISeq
[here](http://www.bioconductor.org/packages/release/bioc/html/NOISeq.html) concerning this.  
The provided file has the higher number always second, but I made the script so that it will work no matter the order the gene end or gene start is provided. It will be more robust this way.  
For ease I'd suggest the `biomart_length.txt`, or your equivalent, be in the same directory with the script. Users can change the name of that file to suit their needs, but must change a line under `USER ADJUSTABLE VALUES` in the script to correspond.  
A file of the output will be saved in the same directory in which the script is run.  

#####example of input and output for `makes_length_annotation_file.py`:

sample of original input in `biomart_length.txt`:
```
Ensembl Gene ID	Transcript Start (bp)	Transcript End (bp)
YHR055C	214533	214718
YPR161C	864449	866422
YOL138C	61325	65350
YDR395W	1263324	1266158
YGR129W	750400	751047
YPR165W	875368	875997
YPR098C	728947	729528
YPL015C	525810	526883
YCL050C	37836	38801
YAL069W	335	649
YMR193W	650036	650812
YLR031W	204225	204785
YIL014C-A	325212	325526
YGR053C	594986	595837
YCR101C	302482	303030
YER087C-B	332582	332830
YOR280C	844992	845792
YGR097W	678695	682135
YHR215W	552099	553502
YKL025C	390240	392279
YKR011C	461635	462696
```

output file produced after:
```
YHR055C	186
YPR161C	1974
YOL138C	4026
YDR395W	2835
YGR129W	648
YPR165W	630
YPR098C	582
YPL015C	1074
YCL050C	966
YAL069W	315
YMR193W	777
YLR031W	561
YIL014C-A	315
YGR053C	852
YCR101C	549
YER087C-B	249
YOR280C	801
YGR097W	3441
YHR215W	1404
YKL025C	2040
YKR011C	1062


```
 
`
`
`
`

 ----------------------------------------------------------------------
 ----------------------------------------------------------------------
 ----------------------------------------------------------------------
CLOSELY RELATED sub-directory
--------------------

I realized `Adjust` didn't suit some of the things I was doing with annotation files and so see the `annotation-utilities` subdirectory [here[(https://github.com/fomightez/sequencework/tree/master/annotation-utilities).


Related scripts
---------------

Several other of my code repositories hold code related to the manipulation of annotations or lists of annotation. Here are some:

- [collection of scripts to Adjust lists in Sequencework repo](https://github.com/fomightez/sequencework/blob/master/Adjust_lists/)

- [Collection of scripts to Extract Details or Annotations in sequencework repo](https://github.com/fomightez/sequencework/tree/master/Extract_Details_or_Annotation)

* `UCSC_chrom_sizes_2_circos_karyotype.py` takes a list of chromosome sizes from UCSC `chrom.sizes` files and reworks it so it can be used as a karyotype file in Circos. It can be found in [my collection of circos-related utility scripts (Python)](https://github.com/fomightez/sequencework/tree/master/circos-utilities).
