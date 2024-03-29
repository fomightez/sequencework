Adjust Lists
=======================

**My collection of files to deal with adjusting lists related to sequence information.**

---
**Description of each script**


- `tx2gene_generator.py`

> a file listing Salmon results --> a file listing transcript to gene names for tximport   
`tx2gene_generator.py` takes a `quant.sf` file from Salmon output where there is a 1:1 transcript: gene relationship and makes the `tx2gene.csv` (comma-separated values) file needed for using tximport to bring Salmon data into DESeq2. By default, the files read and saved will be those, respectively, but they can be specified as arguments when executing the script. See [here](https://www.bioconductor.org/help/workflows/rnaseqGene/) and [here](https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html) for information about tximport use with DESeq2.  
In order to keep things simple with the argument parsing system, one has to specify both input and output file if looking to not use the standard default for output.  
Even though it is stated as being made to work with Salmon's `quant.sf`-type output, it will actually work with any gene/transcript list where there is a single header/column names line first and the gene/transcript designation is the first word or column (tab-separated) for pretty much anything, but comma-separated input. It could be easily edited to work to take a comma-separated list by editing the `split` command. The caveat here is this is made to work where there is a 1:1 relationship with transcripts to genes as there is generally assumed in practice with early branching eukaryotes, like for baker's/budding yeast.
On a even more generalized note, this script is a good for reading in one list file and repuprosing it to another form and saving output and uses command line arguments as file names, and so it would be a good place to start for modifying for other uses. In particular, the argparse handling of the file names and use of `with open` is better/more pythonic at present.

 ```
usage: tx2gene_generator.py [-h] [INPUT_FILE] [OUTPUT_FILE]  

tx2gene_generator.py takes a `quant.sf` file from Salmon output and makes the
`tx2gene.csv` file needed for using tximport to bring Salmon data into DESeq2.
By default, the files read and saved will be those, respectively, but they can
be specified as arguments when executing the script. **** Script by Wayne
Decatur (fomightez @ github) ***  

positional arguments:
  INPUT_FILE   **OPTIONAL**Name of the file generated by Salmon when run with
               your transcriptome of interest. Usually, this is 'quant.sf' &
               if no input file name is provided then this will be used by
               default.
  OUTPUT_FILE  **OPTIONAL**Name of file to save results. If BOTH input and
               output file are not provided, 'tx2gene.csv', will be used.  
               
optional arguments:
  -h, --help   show this help message and exit
 ```
 
 
##### example of input and output for `tx2gene_generator.py  `:
 
**input:**  
Salmon results, `quant.sf` file


**typical commands:**

  python tx2gene_generator.py  
  
    -or-
    
  python tx2gene_generator.py quant.sf tx2geneSPECIAL.csv

**output after run:**  
(text in a file, called `tx2gene.csv `)

 
 ------------------------


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
  FILE        Name of file containing `XLOC` list to fix. REQUIRED.  
  
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

  python cufflinksIDs_list_to_systematic_names.py gene_list.txt  

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
Commonly you would use this script in conjunction with `cufflinksIDs_list_to_systematic_names.py` that is also in this repository and requires an associated gtf file.

> A related script that is independent of any gtf file and instead employs YeastMine data in order to convert systematic names to yeast standard (common) gene names is available, see `genes_in_list_with_SGD_Systematic_Name_to_standard_name.py` in the [YeastMine repository](https://github.com/fomightez/yeastmine). A script to aid in comparing results between the two, called `compare_results_systematic_to_std.py`, is found in the `Evaluation` subdirectory in the [YeastMine repository](https://github.com/fomightez/yeastmine).

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

  python  systematic_names_to_standard_names_using_cufflinks_gtf.py gene_list.txt  

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
 
 
 
 
 
My Related scripts
---------------

Several other of my code repositories hold code related to the manipulation of gene or genetic information lists. Here are some:

- [My YeastMine repository](https://github.com/fomightez/yeastmine) has several useful scripts for converting yeast gene lists, such as `geneID_list_to_systematic_names.py`.

- [My text-mining/text manipulation code repository](https://github.com/fomightez/text_mining) has several useful scripts involving lists, such as, the `find_overlap_in_lists.py` script (and related `find_overlap_in_lists_with_Venn.py`) and `extract_data_on_line_using_word_list.py`.

* `UCSC_chrom_sizes_2_circos_karyotype.py` takes a list of chromosome sizes from UCSC `chrom.sizes` files and reworks it so it can be used as a karyotype file in Circos. It can be found in [my collection of circos-related utility scripts (Python)](https://github.com/fomightez/sequencework/tree/master/circos-utilities).


Related code by others
----------------------

From [Ming Tang early 2019](https://twitter.com/tangming2005/status/1087776160463888384) using `awk` to create a tx2gene file/`tx2gene.csv` from an ensembl gtf:

>"create a tx2gene mapping file from ensemble gtf retaining the version number of genes and transcripts.
awk -F "\t" '$3 == "transcript" { print $9 }' myensembl.gtf| tr -s ";" " "   | cut -d " " -f2,4|  sed 's/\"//g' | awk '{print $1"."$2}' > genes.txt"
"awk -F "\t" '$3 == "transcript" { print $9 }' myensembl.gtf| tr -s ";" " "   | cut -d " " -f6,8|  sed 's/\"//g' | awk '{print $1"."$2}' > transcripts.txt

paste transcripts.txt genes.txt > tx2genes.txt"


Related ideas by others
-----------------------


From https://twitter.com/DrAnneCarpenter/status/1193521158697734145    November 2019
>"I have a list of 167 genes and I need to choose ~30 of them for an experiment. I want the most diverse set, i.e., not a bunch of genes in the same pathways.
Ideas for how to do this expediently?"

Small sampling of sme answers that intrguided me:

  >"What I often do:
  1. Pathways & master regulator(MR) analysis.
  2. Rank the pathways.
  3. Extract MR in each of those pathways.
  MR can be also via Motif analysis to obtain TFBS. One can also obtain network hubs & rank them based on disease prioritization/drug response signatures."
  >"If the tools or packages are not narrowed down. I will suggest to take a look at @MaayanLab website for their amazing resources to use.http://labs.icahn.mssm.edu/maayanlab/resources/  And in regulon or MR space I am biased with VIPER https://bioconductor.org/packages/release/bioc/vignettes/viper/inst/doc/viper.pdf + scRNASeq GRNs if cell type info is needed"  
  >"I would define a set scoring function that would reward the individual importance of each gene but penalize the similarity (similarity could be defined based on shared pathways, GO semantic similarity, etc) and then use the greedy algorithm with that scoring function.1/2"
  >"We used a similiar idea here: https://arxiv.org/abs/1812.02497 ) and here: https://dx.doi.org/10.1109/TCBB.2019.2935437  but the selection is not on genes.  2/2"  
  >"Cc @jmschreiber91 use Jacob's submodular optimization tool that allows you to do exactly that"  
  >"Here is the tool https://github.com/jmschrei/apricot . It's a general method for selecting optimal representative subsets"  
  >"Take 30 random. If your orininal set, as you seem to suspect, is very redundaut then you are probably ok . Once you have your 30 you use the remaining time allocated to this task to do some due diligence on them and perhaps replace some obvious redundancies."  
  >"I’d pick  @XDarzacq’s idea for the quick&dirty first pass. However, ideal method would be something akin to @ilyakorsunsky ’s suggestion; find the 30 genes least correlated to others in the set, THEN swap out any that had little/no prior knowledge."  

