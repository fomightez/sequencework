Adjust Lists
=======================

**My collection of files to deal with adjusting lists related to sequence information.**

---
**Description of each script**


- `cufflinksIDs_list_to_systematic_names.py`

> a file listing gene ids --> a file listing systematic gene ids  
`cufflinksIDs_list_to_systematic_names.py` uses a `gtf` file that is output from Cufflinks to convert a list of `XLOC_`-style gene ids in a file to systematic gene ids.  
The list should be gene ids each on a separate line of the file. 
As written, it needs a `gtf` annotation file called `merged.gtf` but this can be changed within the script under `USER ADJUSTABLE VALUES `.  
This script was written for use with yeast gene annotation file but should work with any `gtf` file from Cufflinks assuming it has entries for `oId`, `gene_id`, and (possibly) `nearest_ref`.   
A file of the output will be saved in the same directory in which the provided gene list file occurs.  

**Usage**  

```
usage: cufflinksIDs_list_to_systematic_names.py [-h] FILE  

cufflinksIDs_list_to_systematic_names.py uses a `gtf` file that is output from
Cufflinks to convert a list of `XLOC_` gene ids in a file to systematic yeast
gene ids. The provided list should be gene ids (`XLOC_` form) each on a
separate line of the file.  
*** Script by Wayne Decatur (fomightez @ github) ***  

positional arguments:
  FILE        Names of file containing `XLOC` list to fix. REQUIRED.  
  
optional arguments:
  -h, --help  show this help message and exit
```

##### example of input and output for `cufflinksIDs_list_to_systematic_names.py`:

first part of annotation file in `merged.gtf`:
```
I	Cufflinks	exon	335	649	.	+	.	gene_id "XLOC_000001"; transcript_id "TCONS_00000001"; exon_number "1"; gene_name "YAL069W"; oId "YAL069W"; nearest_ref "YAL069W"; class_code "="; tss_id "TSS1"; p_id "P2";
I	Cufflinks	exon	538	792	.	+	.	gene_id "XLOC_000001"; transcript_id "TCONS_00000002"; exon_number "1"; gene_name "YAL068W-A"; oId "YAL068W-A"; nearest_ref "YAL068W-A"; class_code "="; tss_id "TSS2"; p_id "P1";
I	Cufflinks	exon	2480	2707	.	+	.	gene_id "XLOC_000002"; transcript_id "TCONS_00000003"; exon_number "1"; gene_name "YAL067W-A"; oId "YAL067W-A"; nearest_ref "YAL067W-A"; class_code "="; tss_id "TSS3"; p_id "P3";
I	Cufflinks	exon	10091	10399	.	+	.	gene_id "XLOC_000003"; transcript_id "TCONS_00000004"; exon_number "1"; gene_name "YAL066W"; oId "YAL066W"; nearest_ref "YAL066W"; class_code "="; tss_id "TSS4"; p_id "P4";
I	Cufflinks	exon	12046	12426	.	+	.	gene_id "XLOC_000004"; transcript_id "TCONS_00000005"; exon_number "1"; gene_name "YAL064W-B"; oId "YAL064W-B"; nearest_ref "YAL064W-B"; class_code "="; tss_id "TSS5"; p_id "P5";
I	Cufflinks	exon	21566	21850	.	+	.	gene_id "XLOC_000005"; transcript_id "TCONS_00000006"; exon_number "1"; gene_name "YAL064W"; oId "YAL064W"; nearest_ref "YAL064W"; class_code "="; tss_id "TSS6"; p_id "P6";
I	Cufflinks	exon	31567	32940	.	+	.	gene_id "XLOC_000006"; transcript_id "TCONS_00000007"; exon_number "1"; gene_name "GDH3"; oId "YAL062W"; nearest_ref "YAL062W"; class_code "="; tss_id "TSS7"; p_id "P7";
I	Cufflinks	exon	33448	34701	.	+	.	gene_id "XLOC_000007"; transcript_id "TCONS_00000008"; exon_number "1"; gene_name "BDH2"; oId "YAL061W"; nearest_ref "YAL061W"; class_code "="; tss_id "TSS8"; p_id "P8";
I	Cufflinks	exon	35155	36303	.	+	.	gene_id "XLOC_000008"; transcript_id "TCONS_00000009"; exon_number "1"; gene_name "BDH1"; oId "YAL060W"; nearest_ref "YAL060W"; class_code "="; tss_id "TSS9"; p_id "P9";
I	Cufflinks	exon	36509	37147	.	+	.	gene_id "XLOC_000009"; transcript_id "TCONS_00000010"; exon_number "1"; gene_name "ECM1"; oId "YAL059W"; nearest_ref "YAL059W"; class_code "="; tss_id "TSS10"; p_id "P10";
I	Cufflinks	exon	37464	38972	.	+	.	gene_id "XLOC_000010"; transcript_id "TCONS_00000011"; exon_number "1"; gene_name "CNE1"; oId "YAL058W"; nearest_ref "YAL058W"; class_code "="; tss_id "TSS11"; p_id "P11";
I	Cufflinks	exon	39259	41901	.	+	.	gene_id "XLOC_000011"; transcript_id "TCONS_00000012"; exon_number "1"; gene_name "GPB2"; oId "YAL056W"; nearest_ref "YAL056W"; class_code "="; tss_id "TSS12"; p_id "P12";
I	Cufflinks	exon	42177	42719	.	+	.	gene_id "XLOC_000012"; transcript_id "TCONS_00000013"; exon_number "1"; gene_name "PEX22"; oId "YAL055W"; nearest_ref "YAL055W"; class_code "="; tss_id "TSS13"; p_id "P13";
I	Cufflinks	exon	45899	48250	.	+	.	gene_id "XLOC_000013"; transcript_id "TCONS_00000014"; exon_number "1"; gene_name "FLC2"; oId "YAL053W"; nearest_ref "YAL053W"; class_code "="; tss_id "TSS14"; p_id "P14";
I	Cufflinks	exon	48564	51707	.	+	.	gene_id "XLOC_000014"; transcript_id "TCONS_00000015"; exon_number "1"; gene_name "OAF1"; oId "YAL051W"; nearest_ref "YAL051W"; class_code "="; tss_id "TSS15"; p_id "P15";
I	Cufflinks	exon	54584	54913	.	+	.	gene_id "XLOC_000015"; transcript_id "TCONS_00000016"; exon_number "1"; gene_name "YAL047W-A"; oId "YAL047W-A"; nearest_ref "YAL047W-A"; class_code "="; tss_id "TSS16"; p_id "P16";
I	Cufflinks	exon	57518	57850	.	+	.	gene_id "XLOC_000016"; transcript_id "TCONS_00000017"; exon_number "1"; gene_name "YAL044W-A"; oId "YAL044W-A"; nearest_ref "YAL044W-A"; class_code "="; tss_id "TSS17"; p_id "P17";
I	Cufflinks	exon	61316	62563	.	+	.	gene_id "XLOC_000017"; transcript_id "TCONS_00000018"; exon_number "1"; gene_name "ERV46"; oId "YAL042W"; nearest_ref "YAL042W"; class_code "="; tss_id "TSS18"; p_id "P18";
I	Cufflinks	exon	62840	65404	.	+	.	gene_id "XLOC_000018"; transcript_id "TCONS_00000019"; exon_number "1"; gene_name "CDC24"; oId "YAL041W"; nearest_ref "YAL041W"; class_code "="; tss_id "TSS19"; p_id "P19";
I	Cufflinks	exon	71786	73288	.	+	.	gene_id "XLOC_000019"; transcript_id "TCONS_00000020"; exon_number "1"; gene_name "CDC19"; oId "YAL038W"; nearest_ref "YAL038W"; class_code "="; tss_id "TSS20"; p_id "P20";
I	Cufflinks	exon	74020	74823	.	+	.	gene_id "XLOC_000020"; transcript_id "TCONS_00000021"; exon_number "1"; gene_name "YAL037W"; oId "YAL037W"; nearest_ref "YAL037W"; class_code "="; tss_id "TSS21"; p_id "P21";
I	Cufflinks	exon	76427	79435	.	+	.	gene_id "XLOC_000021"; transcript_id "TCONS_00000022"; exon_number "1"; gene_name "FUN12"; oId "YAL035W"; nearest_ref "YAL035W"; class_code "="; tss_id "TSS22"; p_id "P22";
I	Cufflinks	exon	79718	80587	.	+	.	gene_id "XLOC_000022"; transcript_id "TCONS_00000023"; exon_number "1"; gene_name "MTW1"; oId "YAL034W-A"; nearest_ref "YAL034W-A"; class_code "="; tss_id "TSS23"; p_id "P23";
I	Cufflinks	exon	82706	83227	.	+	.	gene_id "XLOC_000023"; transcript_id "TCONS_00000024"; exon_number "1"; gene_name "POP5"; oId "YAL033W"; nearest_ref "YAL033W"; class_code "="; tss_id "TSS24"; p_id "P24";
I	Cufflinks	exon	84669	84977	.	+	.	gene_id "XLOC_000024"; transcript_id "TCONS_00000025"; exon_number "1"; gene_name "YAL031W-A"; oId "YAL031W-A"; nearest_ref "YAL031W-A"; class_code "="; tss_id "TSS25"; p_id "P25";
```

sample gene list (content in file to be specified when script called):
```
XLOC_000058
XLOC_000318
XLOC_000427
XLOC_000594
XLOC_000652
XLOC_000831
XLOC_001067
XLOC_001248
XLOC_002503
XLOC_002571
XLOC_004181
XLOC_006140
XLOC_006597
```


**command:**

    cufflinksIDs_list_to_systematic_names.py gene_list.txt  

**output after run:**  
(text in a file, called `gene_list_converted_to_sys_id.txt `, with the contents below)
```
YAR053W
YBR250W
YBR029C
YCL055W
YCR081W
YDL154W
YDR290W
YDL159C-B
YGL168W
YGL053W
YKR089C
YOL083C-A
YPR100W

```
 
 
------------------------
 - `systematic_names_to_standard_names_using_cufflinks_gtf.py`

> a file listing yeast gene ids in systematic form --> a file listing standard (common) gene names, where they exist    
`systematic_names_to_standard_names_using_cufflinks_gtf.py` uses a `gtf` file that is output from Cufflinks to convert a list of yeast systematic gene ids in a file to standard (common) gene names, where they exist.  
The list should be systematic gene ids each on a separate line of the file. 
As written, it needs a `gtf` annotation file called `merged.gtf`, but this can be changed within the script under `USER ADJUSTABLE VALUES `.  
This script was written for use with yeast gene annotation file but should work with any `gtf` file from Cufflinks assuming it has entries for `oId`, `gene_id`, `gene_name`, and (possibly) `nearest_ref`.   
A file of the output will be saved in the same directory in which the provided gene list file occurs.  
Since this is a standard conversion relying on SGD-specific data that just happens to also be already in the gtf file, I plan to make a version similar to this that doesn't require the gtf file but instead access data at YeastMine to allow the conversion.

 ```
 usage: systematic_names_to_standard_names_using_cufflinks_gtf [-h] FILE  
 
systematic_names_to_standard_names_using_cufflinks_gtf.py uses a `gtf` file
that is output from Cufflinks to convert a list of systematic gene ids in a
file to standard (common) gene names, where they exist. The list should be
gene ids each on a separate line of the file. 
**** Script by Wayne Decatur (fomightez @ github)***  

positional arguments:
  FILE        Names of file containing `systematic ids` list to convert.
              REQUIRED.  
              
optional arguments:
  -h, --help  show this help message and exit
 ```
 
 
##### example of input and output for `cufflinksIDs_list_to_systematic_names.py`:
 
see above about needed file `merged.gtf`.

sample gene list (content in file to be specified when script called):
```
YAR053W
YBR250W
YBR029C
YCL055W
YCR081W
YDL154W
YDR290W
YDL159C-B
YGL168W
YGL053W
YKR089C
YOL083C-A
YPR100W
```


**command:**

    systematic_names_to_standard_names_using_cufflinks_gtf.py gene_list.txt  

**output after run:**  
(text in a file, called `gene_list_converted_to_std_id.txt `, with the contents below)
```
YAR053W
SPO23
CDS1
KAR4
SRB8
MSH5
YDR290W
YDL159C-B
HUR1
PRM8
TGL4
YOL083C-A
MRPL51

```
 
 ------------------------
 
 
 
Related scripts
---------------

Several other of my code repositories hold code related to the manipulation of gene lists. Here are some:

- [My YeastMine repository](https://github.com/fomightez/yeastmine) has several useful scripts for comverting yeast gene lists, such as `geneID_list_to_systematic_names.py`.

- [My text-mining/text manipulation code repository](https://github.com/fomightez/text_mining) has several useful scripts involving lists, such as, the `find_overlap_in_lists.py` script and related `find_overlap_in_lists_with_Venn.py`.


