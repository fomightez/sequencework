Sequence Converters
===================

- ConvertFASTAdnaSEQtoRNA.py      DNA --> mRNA

>ConvertFASTAdnaSEQtoRNA.py takes as input a file of nucleotide sequences in FASTA format from NCBI and changes all 'T's to 'U's in the sequence, effectively changing it to mRNA sequence. Importantly, it doesn't alter the description line of the FASTA entries, i.e., the line beginning '>' before the sequence. Note there is no check if these are protein or nucleic sequences. Thus it will change T to U in protein if you try to give it protein sequences.

#####example of input and output for ConvertFASTAdnaSEQtoRNA.py:

original input:
```
> 
```

output after:
```
> 


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
 `
 `

- Another scipt?
