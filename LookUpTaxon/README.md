Taxon Lookup Scripts
===================

- LookUpTaxonFA.py   

>LookUpTaxonFA.py takes as input a file of sequences in FASTA format from NCBI and gives a report of the taxons involved. It laso makes an output file with each taxon added to the description line of each FASTA record.

- LookUpTaxonFAComD.py   

>LookUpTaxonFAComD.py takes as input a file of sequences in a modified FASTA format where my BioalignerNamer script has added the common name of the organisms and gives a report of the taxons involved. It also makes an output file with each taxon added to the description line of each FASTA record.

example of input and output for LookUpTaxonFAComD.py:
before:
```
>chicken |gi|45383616|ref|NP_989588.1| follitropin subunit beta [Gallus gallus]
---------------------------------------M-----------KTLN-----
-CYV-L----------------------------------L--FCWKA-ICCYS-CELTN
I-TIAVERE-----ECELCITVNATWCSGYCFTRD-PVYKYPP--V------------SS
VQQICTFKEV--VY----ETVKIPGCGDHPE---SFYSYPVATECHCETCDTDSTDCTVR
GLGP--SYCSFSHNGSNQ------------------------------------------

```

after:
```
>chicken |gi|45383616|ref|NP_989588.1| follitropin subunit beta [Gallus gallus]|Aves
---------------------------------------M-----------KTLN-----
-CYV-L----------------------------------L--FCWKA-ICCYS-CELTN
I-TIAVERE-----ECELCITVNATWCSGYCFTRD-PVYKYPP--V------------SS
VQQICTFKEV--VY----ETVKIPGCGDHPE---SFYSYPVATECHCETCDTDSTDCTVR
GLGP--SYCSFSHNGSNQ------------------------------------------
-
```
Related report also produced not shown here.
