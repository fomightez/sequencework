Adjust Annotation Files
=======================

**My collection of files to deal with adjusting sequence annotation files.**

---
**Description of each script**


- cufflinksIDs_list_to_systematic_names.py

> a file listing gene ids --> a file listing systematic gene ids
`cufflinksIDs_list_to_systematic_names.py` uses a `gtf` file that is output from Cufflinks to convert a list of `XLOC_`-style gene ids in a file to systematic gene ids.  
# The list should be gene ids each on a separate line of the file. 
As written, it needs a `gtf` annotation file called `merged.gtf` but this can be changed within the script under `USER ADJUSTABLE VALUES `.  
This script was written for use with yeast gene annotation file but should work with any `gtf` file from Cufflinks assuming it has entries for `oId`, `gene_id`, and (possibly) `nearest_ref`.   
A file of the output will be saved in the same directory in which the provided gene list file occurs.  

#####example of input and output for makes_length_annotation_file.py:

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
 
 
Related scripts
---------------

Several other of my code repositories hold code related to the manipulation of gene lists. Here are some:

- [My YeastMine repository](https://github.com/fomightez/yeastmine) has several useful scripts for comverting yeast gene lists, such as `geneID_list_to_systematic_names.py`.

- [My text-mining/text manipulation code repository](https://github.com/fomightez/text_mining) has several useful scripts involving lists, such as, the `find_overlap_in_lists.py` script and related `find_overlap_in_lists_with_Venn.py`.


