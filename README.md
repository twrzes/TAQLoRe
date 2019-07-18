# TAQLoRe
TAQLoRe - Transcript Annotation and Quantification using Long Reads

## Documentation

Full documentation is available at https://taqlore.readthedocs.io.

## Description

**TAQLoRE** is a Snakemake-based pipeline to improve existing annotations and to quantify transcripts coming from long read amplicon-based cDNA sequencing technologies (Oxford Nanopore Technologies, PacBio). It was tested on Linux (CentOS 6).

Briefly, it uses LAST to align all reads to the transcriptome, then it discovers new exons by looking at insertions in alignments, it creates meta-gene with all known and novel exons, aligns all reads to it and generates a TMM-normalised read counts, together with expression heatmaps and PCA plots. It also identifies new splice sites by looking at perfectly aligned reads to the genome, and correcting all splice sites to the closest most abundant canonical ones.

## Authors

### Developers:

- Wilfried Haerty (Earlham Institute)
- Tomasz Wrzesinski (Earlham Institute)

### Contributors:

- Elizabeth Tunbridge (University of Oxford)
- Michael Clark (University of Melbourne)
- Nicola Hall (University of Oxford)
- Syed Hussain (University of Oxford)
- Hami Lee (University of Oxford)
