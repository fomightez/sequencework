Adjust Annotation Files
=======================

**My collection of files to deal with adjusting sequence annotation files.**

---

- makes_length_annotation_file.py

> annotation file --> differently-arranged annotation file  
`makes_length_annotation_file.py` takes a file of `biomart_length.txt` obtained from from [here](http://useast.ensembl.org/info/data/biomart/index.html) and returns a file with gene_id
and then gene length for use with NOISeq. See the User's Guide for NOISeq
[here](http://www.bioconductor.org/packages/release/bioc/html/NOISeq.html) concerning this.  
The provided file has the higher number always second, but I made the script so that it will work no matter the order the gene end or gene start is provided. It will be more robust this way.  
For ease I'd suggest the `biomart_length.txt`, or your equivalent, be in the same directory with the script. Users can change the name of that file to suit their needs, but must change a line under `USER ADJUSTABLE VALUES` in the script to correspond.  
A file of the output will be saved in the same directory in which the script is run.  

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
