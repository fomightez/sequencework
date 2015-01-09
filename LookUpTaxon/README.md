Taxon Lookup Scripts
===================

- LookUpTaxonFA.py   

>LookUpTaxonFA.py takes as input a file of sequence records in FASTA format from NCBI and gives a report of the taxons involved. It laso makes an output file with each taxon added to the description line of each FASTA record.  The FASTA record can be modified, for example include alignment data, as the script only deals with the descriptor line and just duplicates the other lines.

#####example of input and output for LookUpTaxonFA.py:

original input:
```
>gi|45383616|ref|NP_989588.1| follitropin subunit beta [Gallus gallus]
MKTLNCYVLLFCWKAICCYSCELTNITIAVEREECELCITVNATWCSGYCFTRDPVYKYPPVSSVQQICT
FKEVVYETVKIPGCGDHPESFYSYPVATECHCETCDTDSTDCTVRGLGPSYCSFSHNGSNQ
```

ouput after:
```
>gi|45383616|ref|NP_989588.1| follitropin subunit beta [Gallus gallus]|Aves
MKTLNCYVLLFCWKAICCYSCELTNITIAVEREECELCITVNATWCSGYCFTRDPVYKYPPVSSVQQICT
FKEVVYETVKIPGCGDHPESFYSYPVATECHCETCDTDSTDCTVRGLGPSYCSFSHNGSNQ


```
#####Related report also produced not shown here.
`
`
`
`

 ----------------------------------------------------------------------
 ----------------------------------------------------------------------
 ----------------------------------------------------------------------
 `
 `
 `
 `

- LookUpTaxonFAComD.py   

>LookUpTaxonFAComD.py takes as input a file of sequence records in a modified FASTA format where my BioalignerNamer script has added the common name of the organisms and gives a report of the taxons involved. It also makes an output file with each taxon added to the description line of each FASTA record.  The FASTA record can be modified, for example include alignment data, as shown below, as the script only deals with the descriptor line and just duplicates the other lines.

#####example of input and output for LookUpTaxonFAComD.py:

original input:
```
>chicken |gi|45383616|ref|NP_989588.1| follitropin subunit beta [Gallus gallus]
---------------------------------------M-----------KTLN-----
-CYV-L----------------------------------L--FCWKA-ICCYS-CELTN
I-TIAVERE-----ECELCITVNATWCSGYCFTRD-PVYKYPP--V------------SS
VQQICTFKEV--VY----ETVKIPGCGDHPE---SFYSYPVATECHCETCDTDSTDCTVR
GLGP--SYCSFSHNGSNQ------------------------------------------

```

output after:
```
>chicken |gi|45383616|ref|NP_989588.1| follitropin subunit beta [Gallus gallus]|Aves
---------------------------------------M-----------KTLN-----
-CYV-L----------------------------------L--FCWKA-ICCYS-CELTN
I-TIAVERE-----ECELCITVNATWCSGYCFTRD-PVYKYPP--V------------SS
VQQICTFKEV--VY----ETVKIPGCGDHPE---SFYSYPVATECHCETCDTDSTDCTVR
GLGP--SYCSFSHNGSNQ------------------------------------------
-
```

######Related report also produced not shown here.
