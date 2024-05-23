# Sequence Retrieval

Repo for my computational resources dealing with retrieving sequences from remote resources.

## Running these scripts
You need to do two things to run these scripts:

- These rely on the Biopython module to run. PythonAnywhere.com has this module [installed by default](https://www.pythonanywhere.com/batteries_included/). To install on [Sourcelair.com](https://www.sourcelair.com/) and many other places, type `pip install biopython` at the terminal. Often ` easy_install biopython` will work at a terminal line as well. See [here](http://biopython.org/DIST/docs/install/Installation.html) for further information about installing the module if the previous suggestion didn't work for your situation. (In regards to PythonAnywhere... as of early February 2015, the version of the Biopython module listed for 2.7 was significantly behind and running `GetmRNAorCDSforProtein.py` actually elicited a warning of along the lines of `ESummary needing DTDs` unless the module updated. To update the module in your own PythonAnywhere account enter on the command line, `pip install --user --upgrade biopython`.) 

- Put your e-mail address in the User_Email variable under USER ADJUSTABLE VALUES. This is because NCBI wants a way to contact you if you are violating their policies.  See [Frequency, Timing and Registration of E-utility URL Requests](http://www.ncbi.nlm.nih.gov/books/NBK25497/).

---

**Description of each script**

- `GetmRNAforProtein.py`      

> protein --> mRNA  
`GetmRNAforProtein.py` takes a list of protein sequence records in FASTA format and gets the mRNA sequence corresponding to each.

**example of input and output for `GetmRNAforProtein.py`:**

original input:
```
>gi|354497310|ref|XP_003510764.1| PREDICTED: glycoprotein hormone alpha-2-like [Cricetulus griseus]
MPMAPRVLLLCLLVLAVIESHCWEAAIPGCHLHSFNVTVRSGRHGTCQGSHVAQACVGHCESSAFPSRHS
VLVASGYRHNITSVSQCCTISSLKKVKVWLSCMGNQRGELEIFTARACQCDMCRLSRY

>gi|379645662|gb|AFD04550.1| gonadotropin common alpha subunit [Cynoglossus semilaevis]
MEGKATAASTMGSVKSATLSLLLLTFSLYVADSYHSKDLQKLGCESCTLGKNDLFSLYGPVYQCQGCCFS
RAFPTPLTTLETMESRKNITSEATCCVARSSYEVVVAGIVVRNHTDCHCSTCKYHKI
```

Command:

    python GetmRNAforProtein.py sequence_list.fa  


Output:
```
>gi|625218702|ref|XM_003510716.2| PREDICTED: Cricetulus griseus glycoprotein hormone alpha 2 (Gpha2), mRNA
CCTAACCCTAGTCTGCCACCAAGCACAGCTTAGTGCCCCCAATGGTGGGACACCCCTTAAACAAAATAAT
CTTGTTCTTGGACTGGGGATGTACCTCCATTGGTTGAGTGCCTCCCTAGCATGCACAGGCCCTGGATTTG
ACCCCAGAACAACATAAACTGGATGTGGGGCTCATGCCTTTTATTCTAGGACTTGGAAGGTAGAGGCAGA
AGGATCAGAAGTTCAAGGTTTCTTTTACTATTACTTTTCTTTTACTATTACTTTTTAGGACAATGTGAAT
CTCTTCTTTTCCCTATTACTTTGCAATCACTAAAAGGCAGATAGGCAACAGGACTGAGAGTAGTCATTTT
TAGACAGAGAATAACTGTTCAAAGAAGAGTCCTGCCCAAGCTCACAGAGTGGGAAGAGGACCCCGCTGCA
CACGTCACACTCCTGAGCTTGCTTCACTCTCTCTAGTACTCTTCCCTTCTGGGGAAGGGAAATGTCAGCA
GATCACTCTCTGAAGAAACCTCAGATTGTGATCCCCAGGTGGGACTGTGCTGAGGCTTCAGGCTAGGAGG
AGCCTGAAAGATGATATCTAGGAAGGCGTGGCTCTAAAACCCAGGCAAATAAAAGCCTGCGGAGCCAGCA
CCAGTTTGGAGACCAGCAGGAGGCACTGGGAGTCTACAGCCGTCTGTTTCTGCTGGACTGGTGAGCTGCT
GAAGCAGGAGACCAGGTGGAGCTGGGAAGGGCAGGAGGGTTGTGGTTAGGGCCTTTTGGCCAAAGCAGCA
TGTGGAGGGGGGAAGAACTTTGTACCCCAGATGGTGTGGGACATTCCAGGGAAAATCTGCACAGAAGGAT
GACAAGTCCCATCCAGAAATGGCTACCCTTTCAGAGTAGCCCAGGCAGAGTGGCACATCCCTACCAAAGA
CAGGATAGACAGCACCACCCAGGCACACCTCACTAGACAGCAGACCCAGCAGCCTGGCCTGCCAGTAGAG
AGTCATGTATTCAAAGACCTCCCCTTCCACCCATTCCTTTTCAGATGCCCATGGCACCACGTGTCCTGCT
CCTCTGCCTGCTGGTTCTGGCAGTCATTGAAAGCCATTGTTGGGAGGCAGCCATCCCGGGCTGCCACTTG
CACTCCTTCAATGTGACAGTGCGAAGTGGCCGCCATGGCACCTGCCAGGGCTCCCATGTAGCGCAGGCCT
GTGTAGGACATTGTGAGTCTAGTGCCTTCCCTTCCCGGCACTCTGTGTTGGTGGCCAGTGGCTATCGGCA
CAACATCACCTCTGTCTCTCAGTGCTGTACCATTAGCAGTCTAAAAAAGGTAAAGGTGTGGCTGAGCTGC
ATGGGGAACCAGCGGGGGGAACTTGAGATCTTTACCGCCAGAGCCTGCCAGTGTGATATGTGCCGTCTCT
CCCGCTACTAGTCCCTGCCTGAAGGGCTCAGGCCCAGGTCCTGCCACTGATATGCTATAGGATCTCTCAA
ATGAGGGGGCTACTTCAGTTCTGACCCTCTTTAGAGTTAGTGAAGATGCCTAGCATT

>gi|379645661|gb|JQ364953.1| Cynoglossus semilaevis gonadotropin common alpha subunit mRNA, complete cds
ACATGGGGAGTGCCAGGAGTTCTCTACAGACGCACCATGGAAGGGAAGGCAACCGCTGCGTCCACAATGG
GCTCGGTGAAATCAGCAACTCTGTCTCTTCTTCTGTTGACCTTTTCTCTTTATGTAGCTGACTCTTACCA
CAGCAAAGACCTACAGAAATTGGGCTGCGAGAGTTGCACTCTGGGAAAGAATGATTTATTCTCACTATAT
GGTCCAGTCTACCAGTGCCAGGGCTGCTGCTTCTCACGAGCGTTCCCCACTCCTCTAACGACATTGGAAA
CGATGGAAAGTCGAAAGAACATCACTTCAGAGGCGACGTGCTGCGTGGCCAGGTCCAGCTATGAGGTAGT
GGTGGCTGGTATTGTGGTGAGAAACCACACAGACTGCCACTGTAGTACCTGTAAATACCACAAGATATGA
CAGAACAGGGGACCACGCTGCAGAGCTCAGCTTCACGGCACATCATTATTCATTTAGTGAACTGTGCCAA
AGATAGTTTTTCTTTTTCAAAAATGTGTTTTAAAACCGATCACACTTTTGCAGTGATCCTTAGTCTGTGA
CGTGTAATTAGTCCACACACTGTCCTTTGTAGATGTACTCTGTGTAATCGATGATGTAATGGAAAGCAAT
TAAATGAAAACATGCATCTTTATAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
```

---

- `GetmRNAorCDSforProtein.py`

> protein --> mRNA or CDS  
`GetmRNAorCDSforProtein.py` takes a list of protein sequence records in FASTA format and gets the mRNA or CDS sequence corresponding to each.

**example of input and output for `GetmRNAorCDSforProtein.py`:**

original input:
```
>gi|354497310|ref|XP_003510764.1| PREDICTED: glycoprotein hormone alpha-2-like [Cricetulus griseus]
MPMAPRVLLLCLLVLAVIESHCWEAAIPGCHLHSFNVTVRSGRHGTCQGSHVAQACVGHCESSAFPSRHS
VLVASGYRHNITSVSQCCTISSLKKVKVWLSCMGNQRGELEIFTARACQCDMCRLSRY

>gi|357618926|gb|EHJ71710.1| glycoprotein hormone alpha 2 [Danaus plexippus]
MFLRNFILLLTLSHLLVAQSYKKPGCHRQGHTRSISIPDCVEFKITTNACRGYCESYSLPSIMLGFKRHP
VTSLGQCCNIMESEDIPVKVLCLDGERNLVFKSAVTCACYHCQKE

>gi|379645662|gb|AFD04550.1| gonadotropin common alpha subunit [Cynoglossus semilaevis]
MEGKATAASTMGSVKSATLSLLLLTFSLYVADSYHSKDLQKLGCESCTLGKNDLFSLYGPVYQCQGCCFS
RAFPTPLTTLETMESRKNITSEATCCVARSSYEVVVAGIVVRNHTDCHCSTCKYHKI
```

Command:

    python GetmRNAorCDSforProtein.py sequence_list.fa  
 
Output:
```
>gi|625218702|ref|XM_003510716.2| PREDICTED: Cricetulus griseus glycoprotein hormone alpha 2 (Gpha2), mRNA
CCTAACCCTAGTCTGCCACCAAGCACAGCTTAGTGCCCCCAATGGTGGGACACCCCTTAAACAAAATAAT
CTTGTTCTTGGACTGGGGATGTACCTCCATTGGTTGAGTGCCTCCCTAGCATGCACAGGCCCTGGATTTG
ACCCCAGAACAACATAAACTGGATGTGGGGCTCATGCCTTTTATTCTAGGACTTGGAAGGTAGAGGCAGA
AGGATCAGAAGTTCAAGGTTTCTTTTACTATTACTTTTCTTTTACTATTACTTTTTAGGACAATGTGAAT
CTCTTCTTTTCCCTATTACTTTGCAATCACTAAAAGGCAGATAGGCAACAGGACTGAGAGTAGTCATTTT
TAGACAGAGAATAACTGTTCAAAGAAGAGTCCTGCCCAAGCTCACAGAGTGGGAAGAGGACCCCGCTGCA
CACGTCACACTCCTGAGCTTGCTTCACTCTCTCTAGTACTCTTCCCTTCTGGGGAAGGGAAATGTCAGCA
GATCACTCTCTGAAGAAACCTCAGATTGTGATCCCCAGGTGGGACTGTGCTGAGGCTTCAGGCTAGGAGG
AGCCTGAAAGATGATATCTAGGAAGGCGTGGCTCTAAAACCCAGGCAAATAAAAGCCTGCGGAGCCAGCA
CCAGTTTGGAGACCAGCAGGAGGCACTGGGAGTCTACAGCCGTCTGTTTCTGCTGGACTGGTGAGCTGCT
GAAGCAGGAGACCAGGTGGAGCTGGGAAGGGCAGGAGGGTTGTGGTTAGGGCCTTTTGGCCAAAGCAGCA
TGTGGAGGGGGGAAGAACTTTGTACCCCAGATGGTGTGGGACATTCCAGGGAAAATCTGCACAGAAGGAT
GACAAGTCCCATCCAGAAATGGCTACCCTTTCAGAGTAGCCCAGGCAGAGTGGCACATCCCTACCAAAGA
CAGGATAGACAGCACCACCCAGGCACACCTCACTAGACAGCAGACCCAGCAGCCTGGCCTGCCAGTAGAG
AGTCATGTATTCAAAGACCTCCCCTTCCACCCATTCCTTTTCAGATGCCCATGGCACCACGTGTCCTGCT
CCTCTGCCTGCTGGTTCTGGCAGTCATTGAAAGCCATTGTTGGGAGGCAGCCATCCCGGGCTGCCACTTG
CACTCCTTCAATGTGACAGTGCGAAGTGGCCGCCATGGCACCTGCCAGGGCTCCCATGTAGCGCAGGCCT
GTGTAGGACATTGTGAGTCTAGTGCCTTCCCTTCCCGGCACTCTGTGTTGGTGGCCAGTGGCTATCGGCA
CAACATCACCTCTGTCTCTCAGTGCTGTACCATTAGCAGTCTAAAAAAGGTAAAGGTGTGGCTGAGCTGC
ATGGGGAACCAGCGGGGGGAACTTGAGATCTTTACCGCCAGAGCCTGCCAGTGTGATATGTGCCGTCTCT
CCCGCTACTAGTCCCTGCCTGAAGGGCTCAGGCCCAGGTCCTGCCACTGATATGCTATAGGATCTCTCAA
ATGAGGGGGCTACTTCAGTTCTGACCCTCTTTAGAGTTAGTGAAGATGCCTAGCATT

>gi|AGBW01005318.1:3193-3286(+) glycoprotein hormone alpha 2 [Danaus plexippus] mRNA cds , corresponds to cds for protein sequence GI_number 357618926
ATGTTCCTTAGAAATTTTATTTTACTTTTAACCCTTAGTCATTTACTCGTTGCACAAAGT
TACAAAAAACCTGGTTGCCACAGACAAGGTCATACAAGAAGTATTAGTATCCCGGATTGT
GTTGAATTCAAAATAACAACAAACGCCTGCCGTGGTTATTGCGAATCATATTCACTACCC
AGCATTATGCTTGGCTTCAAAAGACATCCCGTTACTTCGCTGGGACAATGCTGCAATATT
ATGGAATCGGAAGACATCCCAGTGAAAGTGCTATGCTTGGACGGAGAGAGAAATTTAGTA
TTTAAATCCGCTGTTACCTGTGCCTGCTACCATTGTCAAAAAGAATAA

>gi|379645661|gb|JQ364953.1| Cynoglossus semilaevis gonadotropin common alpha subunit mRNA, complete cds
ACATGGGGAGTGCCAGGAGTTCTCTACAGACGCACCATGGAAGGGAAGGCAACCGCTGCGTCCACAATGG
GCTCGGTGAAATCAGCAACTCTGTCTCTTCTTCTGTTGACCTTTTCTCTTTATGTAGCTGACTCTTACCA
CAGCAAAGACCTACAGAAATTGGGCTGCGAGAGTTGCACTCTGGGAAAGAATGATTTATTCTCACTATAT
GGTCCAGTCTACCAGTGCCAGGGCTGCTGCTTCTCACGAGCGTTCCCCACTCCTCTAACGACATTGGAAA
CGATGGAAAGTCGAAAGAACATCACTTCAGAGGCGACGTGCTGCGTGGCCAGGTCCAGCTATGAGGTAGT
GGTGGCTGGTATTGTGGTGAGAAACCACACAGACTGCCACTGTAGTACCTGTAAATACCACAAGATATGA
CAGAACAGGGGACCACGCTGCAGAGCTCAGCTTCACGGCACATCATTATTCATTTAGTGAACTGTGCCAA
AGATAGTTTTTCTTTTTCAAAAATGTGTTTTAAAACCGATCACACTTTTGCAGTGATCCTTAGTCTGTGA
CGTGTAATTAGTCCACACACTGTCCTTTGTAGATGTACTCTGTGTAATCGATGATGTAATGGAAAGCAAT
TAAATGAAAACATGCATCTTTATAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
```

## Impending issues and to do
- These scripts may have issues working as written once the GI sequence identifiers get phased out as described in section 1.4.1 of the [June 15 2015 Distribution Release Notes for NCBI-GenBank Flat File Release 208.0](ftp://ftp.ncbi.nih.gov/genbank/gbrel.txt) that can be found at ftp://ftp.ncbi.nih.gov/genbank/gbrel.txt. The GI accession form will be phased out by mid 2016 according to those release notes.
- I should work on switching them over to the Accession.Version system. 

# Related scripts by me
For retrieving S. cerevisiae S288C sequences via YeastMine, see the following:

* [`get_protein_seq_as_FASTA.py` in my 'YeastMine' code repository](https://github.com/fomightez/yeastmine)
* [`get_gene_genomic_seq_as_FASTA.py` in my 'YeastMine' code repository](https://github.com/fomightez/yeastmine)
* [`get_chromosomal_coordinates_as_FASTA.py` in my 'YeastMine' code repositor](https://github.com/fomightez/yeastmine)

For retrieving/extracting sequences from 'local' sources (meaning where I supply most of sequences themselves), see:

* ['Sequencework/Extract_from_FASTA' code repository](https://github.com/fomightez/sequencework/tree/master/Extract_from_FASTA/)
* ['Sequencework/FindSequence' code repository](https://github.com/fomightez/sequencework/tree/master/FindSequence/)
* ['Sequencework/alignment-utilities' code repository](https://github.com/fomightez/sequencework/tree/master/alignment-utilities/)


I have coded several for extracting various sequences or subsequences from FASTA files using biopython. For example, see `drafting_extracting_mito_chr_function.py`. 

# Related scripts/APIs/Tools by others

- Haibao Tang's jcvi/MCscan software, see [here for links and info](https://github.com/fomightez/mcscan-binder), contains a script that
collect sequennces matching features annotated in gff3 files:

    >" Parses the selected features out of GFF, with subfeatures concatenated.
        For example, to get the CDS sequences, do this:
        $ gff.py load athaliana.gff athaliana.fa --parents mRNA --children CDS"

    It can be called with `!python2 -m jcvi.formats.gff load` in a Jupyter envioronment with jcvi/MCscan installed (easily launched from [here](https://github.com/fomightez/mcscan-binder)). Use that command alone to get 'help' information/usage from which that quote above i taken.
    
- [genomepy](https://github.com/vanheeringen-lab/genomepy)
>"Install and use genomes & gene annotations the easy way!  
genomepy is designed to provide a simple and straightforward way to download and use genomic data. This includes (1) searching available data, (2) showing the available metadata, (3) automatically downloading, preprocessing and matching data and (4) generating optional aligner indexes. All with sensible, yet controllable defaults. Currently, genomepy supports UCSC, Ensembl and NCBI." - Includes an S. cerevisiae example.

- I came across someone else's (Leighton Pritchard) take on getting coding sequences for proteins with the protein sequence in early 2018:
Maybe useful stuff or better? See [here](https://twitter.com/widdowquinn/status/963148471334158336)

    >"Need to work backwards from proteins to their coding sequences? Tired of trawling through the database to get the CDS? Me, too. So I wrote this: https://widdowquinn.github.io/ncfp/     I hope you find it useful. Bugs, issues etc. here, please: https://github.com/widdowquinn/ncfp/issues …"


- [pyani](https://github.com/widdowquinn/pyani) has a script called `genbank_get_genomes_by_taxon.py` that  enables download of genomes from NCBI, specified by taxon ID.

- [genomepy](https://github.com/vanheeringen-lab/genomepy) See [here](https://twitter.com/svheeringen/status/1260085341198913538).

        >"Want a way to easily download genomes from Ensembl, UCSC, NCBI or any specific URL? Automatic installation of annotation? Easy access of genome files and accompanying files via Python API?"
        >"Get version 0.8.1 via conda or pip! Use in your pipelines! Recommend it to your colleagues if you like it! Let us know what could be improved: https://github.com/vanheeringen-lab/genomepy/issues! Lots of work behind the scenes, the majority of the codebase is now covered by tests."
        >"Visible changes include:
        White heavy check markImproved download (and sanity checking) of annotations.
        White heavy check markSearch by NCBI taxonomy ID.
        White heavy check markCleaner search output, which now includes assembly accessions.
        White heavy check markMore complete metadata in the genome README."  
  There's also a snamkemake wrapper for genomepy, see [here](https://twitter.com/svheeringen/status/1300664254840950784) 
    >"Download genomes and annotation as part of your snakemake workflow with the genomepy wrapper! https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/genomepy.html "

- [`ffq` (Fetch FastQ)](https://github.com/pachterlab/ffq) is a command line tool for finding sequencing data from public databases. "ffq receives an accession and returns the metadata for that accession as well as the metadata for all downstream accessions following the connections between GEO, SRA, EMBL-EBI, DDBJ, and Biosample" - SOURCE: https://github.com/pachterlab/ffq

- [iSeq](https://github.com/BioOmics/iSeq): An integrated tool to fetch public sequencing data and metadata. It is a BASH script installable by `conda`.
