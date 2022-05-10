# Alignment utilities

Repo for my own computational resources dealing with aligning and alignments (both nucleic and protein), plus give me a place to reference other handy resources.

### Description of the utlity scripts:

* add_actual_position_numbering_to_sequence_in_top_line_of_clustal_output.py
> alignment text --> alignment text with top sequence in alignment blocks annotated with numbering according to actual contiguous sequence (i.e., excluding gaps)

info here

* MSA_to_corresponding_residue_numbers.py
> alignment --> information about individual aligning positions from the MSA

There is a [demo notebook for this script as part of my series of structurework-command line demos](https://github.com/fomightez/cl_demo-binder). You can launch the series from [here](https://github.com/fomightez/cl_demo-binder) and then selecting from the index to go to the 'Determine residues that match to a reference from MSA and use to construct fit commands' page. The direct link to a nicely-rendered, static version of that page is [here](https://nbviewer.org/github/fomightez/cl_demo-binder/blob/master/notebooks/Determine%20residues%20that%20match%20to%20a%20reference%20from%20MSA%20and%20use%20to%20construct%20fit%20commands.ipynb). Usage of the script is illustrated there because this script is ideal for bridging sequence analysis to molecular structure analysis and I had an idea for a good illustration of that bridging concept. That demo also has a functional, fuller accounting of documentation for this script as well.

* BIOALIGNAMER (* NOT HERE YET BUT BELONGS HERE*)
>alignment text with cryptic identifiers for organisms --> alignment text with common names as ids

I'd need to see if it needs updating/fixing (for example, I think it is purely 2.7) but I have one of the more fleshed out Python scripts I wrote that corrects names in alignments to add the common name automagically. Created for facilating pipelines for submitting tasks to the GARLI webserver.
I intended to put up a web-served version of it [here](http://fomightez.pythonanywhere.com/) but I recall I hadn't figured out how to get the toggles to work at the bottom at the time for controlling the command line parameters. (I know I got some text entry forms working later, for example see [here](http://fomightez.pythonanywhere.com/ammonium_screen/) and [here](http://fomightez.pythonanywhere.com/spartan_fixer/) and [here](http://fomightez.pythonanywhere.com/triplet_freq/), and so I probably could figure it out now.)  
Something to explore related to this is if [ETE Toolkit:Dealing with the NCBI Taxonomy database](http://etetoolkit.org/docs/latest/tutorial/tutorial_ncbitaxonomy.html) is any easier to deal with than Entrez for getting the common names.  
Related to my BIOALIGNNAMER but more limited scope: [FungusNameCheck.R](https://github.com/Rowena-h/FungusNameCheck) by Rowena Hill (@Rowena-h) - "A script to check fungi names for the most up to date nomenclature according to [Species Fungorum](http://www.speciesfungorum.org/)."


* pretty_msa_maker_from_clustal_nucleic.py (* NOT HERE YET BUT SHOULD BE HERE*)
> alignment in text form --> alignment fully annotated with features and numbers in vector graphics form

info here

* extract_regions_from_clustal_alignment.py
> alignment of large region ---> alignment of a portion, plus saves sequence of ungapped corresponding block as well

Intended to be used when aligning the entire chromosome along full length, end-to-end to insure corresponding sub-regions best aligned (assuming no rearrangements). Under 'Whole genome alignment', about 4/5 through set of slides, [here](http://www.utdallas.edu/~pradiptaray/teaching/1_mult_seq_align.pdf) refers to this method as '• Identify “collinear” (orthologous) regions or blocks and perform piecewise alignment'.)

info here.


* reverse_complement_of_clustal_alignment.py
> alignment ---> alignment with same sequences in reverse complement form

Intended to be used after `extract_regions_from_clustal_alignment.py` for sequences on the other (a.k.a. reverse) strand because some might be on other strand when extract from entire chromosome and sometimes still it would make more sense to view the extracted sequence in the other direction. For example, if typically used to looking at the regions in one direction, etc.

info here.

* score_columns_in_clustal_msa.py
>alignment ---> two scores assessing conservation of an MSA overall in a quick-n-dirty manner

Intended to be used when comparing separate alignments of related but distinct genetic elements or protein paralogs/orthologs/gene family members/gene cluster members in order to objectively assess which of the related elements are more highly conserved overall.

info here

* calculate_cons_for_clustal_nucleic.py
> nucleic alignment without consensus line ---> alignment with consensus symbols line

Meant to add a consensus symbol line for cases where alignment produced doesn't have it (e.g., Clustal output from Stretcher) or it has been lost. Consensus symbols line it produces is similar to the Clustal style consensus symbols MUSCLE adds below each reside of the alignment. 
info here. 
(Keep in mind [MView](https://www.ebi.ac.uk/Tools/msa/mview/) good if just need visual representation, albeit not as concise as this script makes or MUSCLE outputs.)



* calculate_cons_for_clustal_protein.py
> protein alignment without consensus line ---> alignment with consensus symbols line

Meant to add a consensus symbol line for cases where alignment produced doesn't have it (e.g., Clustal output from Stretcher or extracted I-TASSER thread template listing or exported from aligned templates at SWISS-MODEL) or it has been lost. Consensus symbols line it produces is similar to the Clustal style consensus symbols MUSCLE adds below each reside of the alignment.  
(Keep in mind [MView](https://www.ebi.ac.uk/Tools/msa/mview/) good if just need visual representation, albeit not as concise, suitable for saving, or machine-parseable as this script makes or MUSCLE outputs.)

The script `calculate_cons_for_clustal_protein.py` is demonstrated in MyBinder sessions you can launch from [here](https://github.com/fomightez/cl_sq_demo-binder).  To see the script in action, go [here](https://github.com/fomightez/cl_sq_demo-binder) and click on `launch binder`. Once the session spins up, select 'Use biopython to make valid CLUSTAL formatted MSAs, check sequence of manually edited alignment, and add consensus line' from the list of available notbeooks. The section with use of `calculate_cons_for_clustal_protein.py` is towards the bottom of the Jupyter notebook under the section 'Add a consensus symbol line to an MSA'.  You can see the Jupyter notebook statically and pertinent section [here](https://nbviewer.jupyter.org/github/fomightez/cl_sq_demo-binder/blob/master/notebooks/Use%20biopython%20to%20make%20valid%20CLUSTAL%20formatted%20MSAs%2C%20check%20sequence%20of%20manually%20edited%20alignment%2C%20and%20add%20consensus%20line.ipynb#Add-a-consensus-symbol-line-to-an-MSA).

* check_seq_in_MSAclustal_consistent_with_FASTA.py
> multiple sequence alignment (clustal) ----> whether or not the sequences in the alignment match sequence provided in FASTA format

This is used to check nothing lost when hand-editing an alignment.   
Demonstrated [here](https://nbviewer.org/github/fomightez/cl_sq_demo-binder/blob/master/notebooks/Use%20biopython%20to%20make%20valid%20CLUSTAL%20formatted%20MSAs%2C%20check%20sequence%20of%20manually%20edited%20alignment%2C%20and%20add%20consensus%20line.ipynb#Check-sequence-of-manually-edited-alignment).

* score_sequences_in_clustal_msa.py 
> multiple sequence alignment (clustal) ---> dataframe of the id and assessment score of each sequence 

This ranks the individual sequences of a multipe sequence alignment for how well well they compare to the others in the alignment. Specifically, it takes a multiple seuence alignment (MSA) in clustal format and generates a quick-n-dirty 'score' of percent match to the 'majority' for each individual sequence in the alignment. Returns a dataframe with the id and score for each sequence in the alignment. Mainly meant to be a rough guide to replace (or confirm) visual inspection to identify what is supposedly 'aligned' doesn't really match up. This situation can arise when looking at small sections from much larger alignments, say for a region from a chromosome or genome alignment. In some regions, the sequence available to align doesn't really match up, but an attempt to align will be made by the aligning software. This script should be useful for recognizing weak scoring sequences from those that appear more reasonable, i.e., represent homologs.

See related script, `score_sequences_in_clustal_msa_favoring_top_line.py`.

RELATED:
see about Torsten Seemann's `snp-dists` below.


* score_sequences_in_clustal_msa_favoring_top_line.py 
> multiple sequence alignment (clustal) ---> dataframe of the id and assessment score of each sequence where emphasis given to matches to top line of alignment

Takes a multiple sequence alignment (MSA) in clustal format and  generates a quick-n-dirty 'score' of percent match to the top line sequence for each individual sequence in the alignment. Mainly meant to be a rough guide to replace (or confirm) visual inspection. Returns a dataframe with the id and score for each sequence in the alignment. It is not exactly relative the top line because gaps in the top line are not counted as conserved and thus even the top line won't be 100% if in the alignment the top line has gaps; this was a compromise to get those sequences that are all gaps to score as `0.0`, which makes more sense then them accumulating score where the gaps happen to match with gaps in top line. Hence `favoring` in the name of the script and not `relative`. This is only meant to be run when `score_sequences_in_clustal_msa.py` fails and you have additional knowledge that the sequence on the top line of the multiple sequence alignment is a good guide, i.e., it is a reference seqeunce. The tell-tale sign that `score_sequences_in_clustal_msa.py` failed is when all sequences have scores of 0.0. `score_sequences_in_clustal_msa.py`scores purely on a match to the majority; however, if the majority of sequences included in an alignment are unrelated to several sequences, the subset of 'related' sequences will be scored poorly even though they really are the subject of interest. This script allows placing the subject of interest on the first line of the mutliple sequence alignment to indicate it is to be considered as 'the reference'. (This can be done, for example, by using `input` as the order in Clustal Omega by selecting that under `More options` at https://www.ebi.ac.uk/Tools/msa/clustalo/ when submitting sequences for alignment.)

Mainly meant to be a rough guide to replace (or confirm) visual inspection to identify what is supposedly 'aligned' doesn't really match up. This situation can arise when looking at small sections from much larger alignments, say for a region from a chromosome or genome alignment. In some regions, the sequence available to align doesn't really match up, but an attempt to align will be made by the aligning software. This script should be useful for recognizing weak scoring sequences from those that appear more reasonable, i.e., represent homologs.

See related script `score_sequences_in_clustal_msa.py`.

RELATED:  
see about Torsten Seemann's `snp-dists` below.


* categorize_residues_based_on_conservation_relative_consensus_line.py
> multiple sequence alignment (clustal with consensus symbols) ---> listings of the positions categorized by conservation for a sequence in the alignment

Meant to be run once you have a multiple sequence alignment and want to use the conservation details for further work, such as authoring commands for molecular visualization of a related structure. 

Takes an multiple sequence alignment (in CLUSTAL format) that has a consensus symbols line, produced from, say, [MUSCLE](https://www.ebi.ac.uk/Tools/msa/muscle/) , and for a specific sequence in the alignment categorizes the residues that are identical, strongly, similar, or 
weakly similar in the alignment.  Plus, the script also categorizes unconserved, while at it. Importantly, residue positions in the results are in common terms where first residue is number one.

There is a [demo notebook for this script as part of my series of structurework-command line demos](https://github.com/fomightez/cl_demo-binder). You can launch the series from [here](https://github.com/fomightez/cl_demo-binder) and then selecting from the index to go to the 'Categorize conservation in a MSA and use that to generate molvis commands' page. The direct link to a nicely-rendered, static version of that page is [here](https://nbviewer.jupyter.org/github/fomightez/cl_demo-binder/blob/master/notebooks/Categorize%20conservation%20in%20a%20MSA%20and%20use%20that%20to%20generate%20molvis%20commands.ipynb). Usage of the script is illustrated there because this script is ideal for bridging sequence analysis to molecular structure analysis and I had an idea for a good illustration of that bridging concept. That demo also has a functional, fuller accounting of documentation for this script as well.

* mview_to_CLUSTAL.py 
> MView output of an alignment with few 'bells-n-whistles' ---> alignment in CLUSTAL format

Meant to be used when wanting to adjust width of sequence characters per line.

Takes a multiple sequence alignment output by MView with options set so mainly just the lines of sequence and not html or consenuss information and reformats it to CLUSTAL.

There is a [demo notebook for this script as part of my series of sequencework-command line demos](https://github.com/fomightez/cl_sq_demo-binder) where this script is demonstrated as part of a notebook combining a typical workflow. To get an active form of the specific demo notebook, you can launch the series from [here](https://github.com/fomightez/cl_sq_demo-binder) and then select from the index to go to the 'Use biopython to make valid CLUSTAL formatted MSAs, check sequence of manually edited alignment, and add consensus line' notebook. Alternatively, the direct link to a nicely-rendered, static version of that page is [here](https://nbviewer.jupyter.org/github/fomightez/cl_sq_demo-binder/blob/master/notebooks/Use%20biopython%20to%20make%20valid%20CLUSTAL%20formatted%20MSAs%2C%20check%20sequence%20of%20manually%20edited%20alignment%2C%20and%20add%20consensus%20line.ipynb).

* roughly_score_relationships_to_subject_seq_pairwise_premsa.py
> sequences in FASTA format (single multi-FASTA file)  ---> quick assessment of similarity of first sequence to each of the others

Meant to be run prior to a multiple sequence alignment task in order to tell quickly what is similar to one another and maybe leave out redundant sequences if you are looking to reduce the submitted alignment job.

Produces a ranked list where the highest scoring one is closest to the first sequence in the provided sequence.

The core of this is based on [this script](https://github.com/berrisfordjohn/adding_stats_to_mmcif/blob/master/adding_stats_to_mmcif/pairwise_align.py) by [John Berrisford](https://github.com/berrisfordjohn).

There is a [demo notebook for this script in this repo](https://github.com/fomightez/sequencework/blob/master/alignment-utilities/demo%20roughly_score_relationships_to_subject_seq_pairwise_premsa.ipynb) that can be viewed nicely displayed [here](https://nbviewer.jupyter.org/github/fomightez/sequencework/blob/master/alignment-utilities/demo%20roughly_score_relationships_to_subject_seq_pairwise_premsa.ipynb).

CAVEAT:  
The original intention was to align the entire length of provided sequences but for moderate- or large-sized sequences (>5 kb) this is not possible unless you have a substantial amount of computing power due to demanding memory resources from Biopython's aligmnet algorithmn (see [here](https://github.com/biopython/biopython/pull/1655)). Thus, it was modified to sample what seems to work when using the script in Jupyter sessions launched from MyBinder.org via [here](https://github.com/fomightez/blast-binder). Therefore, the reported score may not represent the entire range of the provided sequences but that is reported back if that is the case.

RELATED:  
see about Torsten Seemann's `snp-dists` below.


* score_differences_between_sequences_by_pairwise_alignment.py
> sequences in FASTA format (single multi-FASTA file)  ---> matrix of difference from pairwise alignment comparisons

See about the impetus for this script [here](https://git.io/fj9ES ). Karin Lagesen [asked](https://twitter.com/karinlag/status/1153676961690071040) about a script that compared seqeuences to each other and I realized I had a similar script already (`roughly_score_relationships_to_subject_seq_pairwise_premsa.py`) if the sequences were *NOT* previously aligned. It takes the sequences and aligns them to each other sequence in turn and calculates the number of differences in the alignment. (Plus, I realized I had a need for such a script now.) The matrix output is based on Torsten Seemann's [`snp-dists`](https://github.com/tseemann/snp-dists) which calculates all the differenes for sequences *already aligned*.  
By default it also produces a heatmap plot of the differences.  
The output should be informative for deciding how to arrange submissions to alignment tools where the order of the input is preserved in the default as it seems with some aligning tools. 

Has some nice features for use in a notebook, too.  
There is a [demo notebook for this script as part of my series of sequencework-command line demos](https://github.com/fomightez/cl_sq_demo-binder). You can launch the series from [here](https://github.com/fomightez/cl_sq_demo-binder) and then selecting from the index to go to the 'Demo of script to calculate differences between sequences in multiFASTA file' page. The direct link to a nicely-rendered, static version of that page is [here](https://nbviewer.jupyter.org/github/fomightez/cl_sq_demo-binder/blob/master/notebooks/Demo%20of%20script%20to%20calculate%20differences%20between%20sequences%20in%20multiFASTA%20file.ipynb).

RELATED:  
see about Torsten Seemann's `snp-dists` below.


- Conversion of MSA in FASTA format to clustal.ipynb
>Conversion of a multiple sequence alignment (MSA) in one format to another based on section 'File Format Conversion' at https://biopython.org/wiki/AlignIO .

notebook in static form: [Conversion of MSA in FASTA format to clustal using Biopython](https://nbviewer.jupyter.org/github/fomightez/cl_demo-binder/blob/master/notebooks/Conversion%20of%20FASTA%20alignment%20format%20to%20clustal%20using%20Biopython.ipynb)

This notbeook is available and can be run by going [here](https://github.com/fomightez/cl_demo-binder), clicking on `launch binder`, and then selecting 'Convert MSA in FASTA format to clustal format using Bioython' from the index listing available notebooks.




Related utilities in my other repositories
------------------------------------------

- `collapse_large_unknown_blocks_in_DNA_sequence.py`

  In order to better judge some sequences extracted from alignments against a chromosome/genome, I've found that reducing apparently arbitrarily sized repeats of unknown nucleotides, represented as uninterrupted repeats of `N`s, such as "NNNNNNNNNNNNNNN" for example, can help in assessing the sequences extracted from a much larger alignment or help prepare them for aligning again, especially when additional knowledge, such as matches in flanking sequence seem to suggest the spacing doesn't match the number of unknown nucleotides shown in the represtation. In other words, often in assemblies generated from Illumina pipelines 'educated guesses' or abtirary numbers (such as 50) `N`s in a row will be introduced for smallish 'gaps' (i.e., ones where software has concluded an intact scaffold is still indicated) in the sequence assembly and reducing these to a smaller size can sometimes make alignments of sub-elements in different rows more obvious. It is important to be clear this was done and keep both sets of data. I have made a script that does this reduction in the size of the blocks; it is called `collapse_large_unknown_blocks_in_DNA_sequence.py` and it can be found in my ['Sequencework/ConvertSeq' code repository](https://github.com/fomightez/sequencework/tree/master/ConvertSeq).

- `get_seq_following_seq_from_multiFASTA.py`

  Often when you use an alignment of entire chromosome to extract a collinear region to identify orthologs sometimes the ends aren't exactly perfect. You may wish to look at the sequence just beyond the edge of collected sequence, and so you can use the end of that extracted seqeunce to search that pattern agains the original source to get the downstream sequence. I made a script to do this post-alignment processing. It is called `get_seq_following_seq_from_multiFASTA.py` and it can be found in my ['Sequencework/Extract_from_FASTA' code repository](https://github.com/fomightez/sequencework/tree/master/Extract_from_FASTA).

- `delete_seq_following_pattern_within_multiFASTA.py`

  Often when you use an alignment of entire chromosome to extract a collinear region to identify orthologs sometimes the ends aren't exactly perfect. You may wish to remove sequences after an element or pattern in order to not have long trailing overhangs beyond the aliged region. I made a script to do this post-alignment processing. It is called `delete_seq_following_pattern_within_multiFASTA.py` and it can be found in my ['Sequencework/AdjustFASTA_or_FASTQ' code repository](https://github.com/fomightez/sequencework/tree/master/AdjustFASTA_or_FASTQ).
  
 - `report_diff_between_two_seq_strings.py`

    This script directly compares two Python strings using an approach meant for aligning biological strings. This is like a simple, direct version of my script `score_differences_between_sequences_by_pairwise_alignment.py`, which resides in this sub-repo.  
  This script is located [here](https://github.com/fomightez/sequencework/tree/master/Compare_biological_seq_strings) in my [compare_biological_seq_strings utilities sub-repo](https://github.com/fomightez/sequencework/tree/master/Compare_biological_seq_strings).
  
 - utility to reformat sequence alignment in result of TM-align structural alignment from all one long sequence to blocks of paired sequences

    Notebook referenced [here in my TM-align-utilities sub-repo](https://github.com/fomightez/structurework/tree/master/tm-align-utilities) demonstrates code to take the alignment TM-align returns as part of its results page and reformats it into blocks of aligned sequences. The reformatted version fits better in other documents & reports and is easier to view & compare the N- and C-terminal sequences at the same time.  
    The residue pairs come from a structural alignment and so this is truly a structure-related utility; however, since it involves involves an alignment featuring sequences I thought I'd leave a link to it here.


Utilities by others
------------------

- [alnvu](https://github.com/nhoffman/alnvu) - Python package to reformat and condense multiple sequence alignments to highlight variability
>"alnvu makes a multiple alignment of biological sequences more easily readable by condensing it and highlighting variability.
Produces formatted multiple alignments in plain text, html, and pdf."  
Related to my scripts `score_columns_in_clustal_msa.py`, `calculate_cons_for_clustal_nucleic.py`, `calculate_cons_for_clustal_protein.py`, `mview_to_CLUSTAL.py`, `categorize_residues_based_on_conservation_relative_consensus_line.py`, and `score_sequences_in_clustal_msa.py`. 

- From a [Biostars reply by Mensur Dlakic to 'Tool to visualize alignment comparisons / Alignment diff'](https://www.biostars.org/p/9505294/#9505319) where OP wanted to quickly visualize the differences between multiple resulting MSAs of the same sequences from different aligners:
>"There is a tool called hhalign from the hh-suite that takes multiple alignments and aligns them into a master alignment."

- [Adding unaligned sequences into an existing alignment using MAFFT and LAST](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3516148/)

I have several listed [here](http://proteopedia.org/wiki/index.php/User:Wayne_Decatur/Sequence_analysis_tools) such as those that interconvert/reformat alignments in various forms. Just search `align` on that page to step through examining them. In particular note:

- Under 'Installable software for fine-tuning sequence alignments' there, I discuss using SEQOTRON on my work Mac and  `cons` alignment consensus program and many others at EMBOSS explorer website to put in conservation. Along the way it looks like it adds numbering similar to Mview where it is just numbering for total number on line INCLUDING GAPS. (Not to self, see my `edited t-coffee alignment for true XXXXX XXXXX.md` for example.)

- A MSA viewer for right in the terminal based on Rich/Textual:   
  https://twitter.com/althonos/status/1524024017325371393   May 2022
  >"55 lines of #Python code to build a Multiple Sequence Alignment viewer for the terminal with the #Rich library (from @textualizeio) #bioinformatics"
  >"If you want to test it, here's a repository, I guess I'll release it to PyPI later on once I've added a few tests and documentation: https://github.com/althonos/rich-msa"

- [trimAl](http://trimal.cgenomics.org/) - a tool for the automated removal of spurious sequences or poorly aligned regions from a multiple sequence alignment. [Corresponding publication for trimAl](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2712344/). Was recommended for " I would still remove the columns that have only gaps (trimAl can do that)."- [Source](https://www.biostars.org/p/9510167/#9510169)

- [The Newick Utilities](https://github.com/tjunier/newick_utils/wiki) 
>"Here's a useful package for working with Newick style #tree files: Newick utilities!
https://github.com/tjunier/newick_utils/wiki
I've used this extensively to prepare files for ASTRAL-III. 
Will incorporate new functions to automate nw utils in next PIrANHA release." (Source: [here](https://twitter.com/justincbagley/status/1163496265100824576)

- [snp-dists](https://github.com/tseemann/snp-dists) - Pairwise SNP distance matrix from a FASTA sequence alignment. 
>"Convert a FASTA alignment to SNP distance matrix."    <--- related to my scripts `roughly_score_relationships_to_subject_seq_pairwise_premsa.py`,`score_sequences_in_clustal_msa.py`, and `score_sequences_in_clustal_msa_favoring_top_line.py`..

-  [HamHeat: Hamming distance calculation from multiple sequence data for heatmap visualization](https://github.com/alexeyrakov/HamHeat)   <--- related to my scripts `score_differences_between_sequences_by_pairwise_alignment.py` and `roughly_score_relationships_to_subject_seq_pairwise_premsa.py`

    >"The script (HamHeat) was applied to the real Salmonella allelic data of our recently published paper Rakov et al., 2019. Figure S1 (additional file 5) from this publication shows the heatmap for 70 virulence factors alleles for 500 Salmonella genomes. For this figure, the HamHeat results for each of the 70 alleles were combined in one file to be used as the input file for the [Morpheus matrix visualization software](https://software.broadinstitute.org/morpheus/)."

- [GToTree](https://github.com/AstrobioMike/GToTree/wiki/What-is-GToTree%3F)
>"GToTree is a program that aims to give more researchers the capability to generate phylogenomic trees to help guide their work. At its heart it just takes in genomes and outputs an alignment and phylogenomic tree based on the specified single-copy gene set."

- [for selecting only differing residues based on an alignment for identifying differences to feed Chimera](https://gist.github.com/idyoung/266427969f99f869bf9aaaf684a1497e)

- [CIAlign](https://github.com/KatyBrown/CIAlign) - "CIAlign allows you to remove sources of noise from an MSA and produce publication-ready visualisations. CIAlign can be used to remove poorly aligned sequence ends, highly divergent sequences, non-majority insertions, short sequences and gaps from a multiple sequence alignment. It also allows the user to plot 'mini alignments' - a novel visualisation method to give an overview of a large MSA in a manageable image size", see [here](https://twitter.com/thekatybrown1/status/1503684156739534849). (Not limited to visualizing nucleic; publication says works with protein alignments as well, although I didn't yet find an example shown.)[Associated publication](https://peerj.com/articles/12983/). [Now availabele through bioconda](https://twitter.com/thekatybrown1/status/1521119123178246144). Also see there about, "And we've added an option to set a threshold proportion of sequences for the remove insertions function."

Somewhat related:  

- [pyani](https://github.com/widdowquinn/pyani) for calculating average nucleotide identity (ANI) and related measures for whole genome comparisons, and rendering relevant graphical summary output.

- [FastANI](https://github.com/ParBLiSS/FastANI) - FastANI is developed for fast alignment-free computation of whole-genome Average Nucleotide Identity (ANI)  
>"FastANI rapidly compares fasta sequences and will output global identity. It is meant only for relatively similar genomes, and will not create any output if identity is below 75-80%." - from Mensur Dlakic's Biostars answer [here](https://www.biostars.org/p/9509271/#9509767)

- There are other ways to score alignments. For now I have been using simplistic assessments based on differences. This can be done much fancier based on probability models, with different scores for sections in between the gaps as well. See [How sequence alignment scores correspond to probability models. Martin C Frith. Bioinformatics, btz576, https://doi.org/10.1093/bioinformatics/btz576](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btz576/5536873) for recent coverage of this, although I don't see an implementation linked to there on my cursory examination.

- [Cool R script to colorfully annotate adapter sequences in a sequence - could be useful to adapt for displaying PCR primers, etc](https://twitter.com/clintcodesbio/status/1339947174239612929); R script [here](https://gitlab.com/gringer/bioinfscripts/-/blob/master/read_annotator.r)
