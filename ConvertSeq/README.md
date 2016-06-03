Sequence Converters
===================

- ConvertFASTAdnaSEQtoRNA.py      DNA --> mRNA

>ConvertFASTAdnaSEQtoRNA.py takes as input a file of nucleotide sequences in FASTA format from NCBI and changes all 'T's to 'U's in the sequence, effectively changing it to mRNA sequence. Importantly, it doesn't alter the [description line](http://blast.ncbi.nlm.nih.gov/blastcgihelp.shtml) of the FASTA entries, i.e., the line beginning '>' before the sequence. Note there is no check if these are protein or nucleic sequences. Thus it will change T to U in protein if you try to give it protein sequences.

#####example of input and output for ConvertFASTAdnaSEQtoRNA.py:

original input:
```
>gi|392920960|ref|NM_073597.2| Caenorhabditis elegans Protein FLR-2 (flr-2) mRNA, complete cds
ATGGGCTCCAAAGCACGAGCACGACGACGTTTAAGTTGTTTTTTAAGCGTTTTTGTTGTGACATGCTTAT
TACAGTACTGCACAGCAGGTGTTACTAAGAATAATAGTTGCAAAAAAGTTGGAGTGGAGGAACTTATAGA
TGAAGAAGGCTGTGATTTGATGATAATTCGAATCAATCGATGCAGTGGGCATTGCTTCTCATTTACATTT
CCTAATCCCTTAACGAAAAAATATTCAGTGCATGCGAAGTGCTGCCGGATGGTTGAATGGGAAATGCTTG
AAACAGAATTAAAATGTTCCAAAGGAAACCGAAATCTTCGAATACCATCTGCAACACAATGTGAATGTTT
TGATTGTCTTGTTCGATAG

>gi|XM_003122551.1:0-390(+) PREDICTED: glycoprotein hormone alpha-2-like [Sus scrofa] mRNA cds , corresponds to cds for protein sequence GI_number 311247365
ATGCCCATGGCCTCCCCCCAAACCCTGCTCCTCTGCCTGCTGGTCCTGGCAATCCCTGAA
GGCCAGGGTCAGCAGGCAGCCATCCCAGGCTGCCACTTGCACCCCTTCAACGTGACCGTG
CGAAGTGACCGCCAAGGCACCTGCCAGGGCTCCCATGTGGCACAGGCCTGTGTGGGCCAC
TGTGAGTCCAGTGCCTTCCCATCCCGGTACTCCGTGCTGGTGGCCAGCGGCTATCGACAC
AACATCACCTCCGTCTCTCAGTGCTGCACCATCAGCAGCCTGAGGAAGGTGAAGGTGCAG
CTGCACTGTGGGGGGGACCGGAGGGAGGAGCTGGAGATCTTCACGGCCAGGGCCTGCCAG
TGCGACATGTGTCGCCTCTCACGCTACTAG
```

output after:
```
>gi|392920960|ref|NM_073597.2| Caenorhabditis elegans Protein FLR-2 (flr-2) mRNA, complete cds
AUGGGCUCCAAAGCACGAGCACGACGACGUUUAAGUUGUUUUUUAAGCGUUUUUGUUGUGACAUGCUUAU
UACAGUACUGCACAGCAGGUGUUACUAAGAAUAAUAGUUGCAAAAAAGUUGGAGUGGAGGAACUUAUAGA
UGAAGAAGGCUGUGAUUUGAUGAUAAUUCGAAUCAAUCGAUGCAGUGGGCAUUGCUUCUCAUUUACAUUU
CCUAAUCCCUUAACGAAAAAAUAUUCAGUGCAUGCGAAGUGCUGCCGGAUGGUUGAAUGGGAAAUGCUUG
AAACAGAAUUAAAAUGUUCCAAAGGAAACCGAAAUCUUCGAAUACCAUCUGCAACACAAUGUGAAUGUUU
UGAUUGUCUUGUUCGAUAG

>gi|XM_003122551.1:0-390(+) PREDICTED: glycoprotein hormone alpha-2-like [Sus scrofa] mRNA cds , corresponds to cds for protein sequence GI_number 311247365
AUGCCCAUGGCCUCCCCCCAAACCCUGCUCCUCUGCCUGCUGGUCCUGGCAAUCCCUGAA
GGCCAGGGUCAGCAGGCAGCCAUCCCAGGCUGCCACUUGCACCCCUUCAACGUGACCGUG
CGAAGUGACCGCCAAGGCACCUGCCAGGGCUCCCAUGUGGCACAGGCCUGUGUGGGCCAC
UGUGAGUCCAGUGCCUUCCCAUCCCGGUACUCCGUGCUGGUGGCCAGCGGCUAUCGACAC
AACAUCACCUCCGUCUCUCAGUGCUGCACCAUCAGCAGCCUGAGGAAGGUGAAGGUGCAG
CUGCACUGUGGGGGGGACCGGAGGGAGGAGCUGGAGAUCUUCACGGCCAGGGCCUGCCAG
UGCGACAUGUGUCGCCUCUCACGCUACUAG


```
 
`
`
`
`

 ----------------------------------------------------------------------
 ----------------------------------------------------------------------
 ----------------------------------------------------------------------
 `
 `
 
In the `RetrieveSeq` folder is a script that takes a list of protein sequence records in FASTA format and gets the mRNA sequence correspoding to each. See `GetmRNAforProtein.py`.
 ----------------------------------------------------------------------
 ----------------------------------------------------------------------
 ----------------------------------------------------------------------
 `
 `
 `
 `

In the `RetrieveSeq` folder is a script that takes a list of protein sequence records in FASTA format and gets the mRNA sequence or at least the coding sequence corresponding to each one. See `GetmRNAorCDSforProtein.py`.


 ----------------------------------------------------------------------
 ----------------------------------------------------------------------
 ----------------------------------------------------------------------
 `
 `
 `
 `
 
In the `UGENE_help` folder in a different repository is some guidance and one script for converting between Unipro's UGENE sequence anysis software and Serial CLoner. See [here](https://github.com/fomightez/UGENE_help) for `UGENE help`.
