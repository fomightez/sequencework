- [degenotate - Annotate degeneracy of sites in coding regions of a genome](https://github.com/harvardinformatics/degenotate) (Learned of it from [this Biostar's question seeking help to calculate degreee of degeneracy](https://www.biostars.org/p/9594758/#9594839).)
  >"degenotate takes as input either a genome FASTA file and a corresponding annotation file (GFF or GTF) OR file or directory of files that contain coding sequences in FASTA format and outputs a bed-like file that contains the degeneracy score (0-, 2-, 3-, or 4-fold) of every coding site." (Zero means "any mutation will change the amino acid" and 4 meands, "four nucleotides at the position code for the same AA, so all 3 possible mutations are synonymous"  

   Vaguely related to my 'score' utilities featured in [sequencework/CompareFASTA_or_FASTQ](https://github.com/fomightez/sequencework/tree/master/CompareFASTA_or_FASTQ), but bringing in the concept of coding potential for coding regions because it scores degeneracy and then for non-degenerate sites it sumarizes mutations and then goes towards comparing when start going to 'MK site counts and tests'.


  
