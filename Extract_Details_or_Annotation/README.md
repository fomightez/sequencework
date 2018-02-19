# Collection of scripts to Extract Details or Annotations

# The scripts

* mine_mito_features_from_transcriptome.py
> transcriptome in FASTA format --> dataframe of mitochdondrial chromsome genes/features

Takes a transcriptome file of entries in FASTA format and mines the 
mitochondrial details from definition lines like below:


```
>Q0065 cdna chromosome:R64-1-1:Mito:13818:21935:1 gene:Q0065 gene_biotype:protein_coding transcript_biotype:protein_coding gene_symbol:AI4 description:Endonuclease I-SceII; encoded by a mobile group I intron within the mitochondrial COX1 gene; intron is normally spliced by the BI4p maturase but AI4p can mutate to acquire the same maturase activity [Source:SGD;Acc:S000007264]
>Q0143 cdna chromosome:R64-1-1:Mito:51277:51429:1 gene:Q0143 gene_biotype:protein_coding transcript_biotype:protein_coding description:Dubious open reading frame; unlikely to encode a functional protein, based on available experimental and comparative sequence data [Source:SGD;Acc:S000007277]
```

Returns a dataframe with the details for the mito features.
Specifically, in this case systematic_id, start, end, midpoint, gene_symbol
(standard name) if there is one noted

* mine_SUT_and merge_into_mito_transcripts_df.py 
> SUTs data and transcriptome data --->  dataframe of mitochdondrial chromsome genes/features

Intended to be used in conjuction with `mine_mito_features_from_transcriptome.py`.  
Returns a dataframe with the details for the mito features, now with yeast mitochondrial SUTS added.

# Related

Put some links here
