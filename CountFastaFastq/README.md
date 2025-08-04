Counting FASTA / FASTQ entries in a multi-fasta fastq file
============

If local and I am just examining where I don't need to progammatically generate a report, I just use REGEX based on start of the description line pattern.  
Such as for FASTQ:

```text
^@.*
```

Or for FASTA:

```text
^>.*
```

Progammatic generating report
-----------

See my Python script `count_fastq_entries.py` here. That also includes notes for bash and awk approaches.
