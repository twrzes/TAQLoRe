###################################################################
TAQLoRE - Transcript Annotation and Quantification using Long Reads
###################################################################

TAQLoRE is a Snakemake-based pipeline to improve existing annotations and to quantify transcripts coming from long read amplicon-based cDNA sequencing technologies (Oxford Nanopore Technologies, PacBio). It was tested on Linux (CentOS 6) but it should work on Mac as well. Briefly, it uses LAST to align all reads to the transcriptome, then it discovers new exons by looking at insertions in alignments, it creates meta-gene with all known and novel exons, aligns all reads to it and generates a TMM-normalised read counts, together with expression heatmaps and PCA plots. It also identifies new splice sites by looking at perfectly aligned reads to the genome, and correcting all splice sites to the closest most abundant canonical ones. For more information, refer to :doc:`description`.

************
Installation
************

- Source code: `GitHub <https://github.com/twrzes/TAQLoRe>`_
- Issue tracker: `Issue tracker <https://github.com/twrzes/TAQLoRe/issues>`_

*********
Citations
*********

Our pipeline is based on following software:

- Snakemake: `Köster J, Rahmann S. “Snakemake - A scalable bioinformatics workflow engine”. Bioinformatics. 2018 Oct 15;34(20):3600. <https://academic.oup.com/bioinformatics/article/28/19/2520/290322>`_
- LAST: `Kiełbasa SM, Wan R, Sato K, Horton P, Frith MC. "Adaptive seeds tame genomic sequence comparison". Genome Res. 2011 Mar;21(3):487-93. <https://genome.cshlp.org/content/21/3/487.long>`_
- GMAP: `Wu TD, Watanabe CK. GMAP: a genomic mapping and alignment program for mRNA and EST sequences. Bioinformatics. 2005 May 1;21(9):1859-75. <https://academic.oup.com/bioinformatics/article/21/9/1859/409207>`_
- Bedtools: `Quinlan AR, Hall IM. BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics. 2010 Mar 15;26(6):841-2. <https://academic.oup.com/bioinformatics/article/26/6/841/244688>`_
- Pybedtools: `Dale RK, Pedersen BS, Quinlan AR. Pybedtools: a flexible Python library for manipulating genomic datasets and annotations. Bioinformatics. 2011 Dec 15;27(24):3423-4. <https://academic.oup.com/bioinformatics/article/27/24/3423/304825>`_

*************
Things to add
*************

- Splice-site-based pipeline (part4 and part5).
- Usage of splice-site-based approach (part4 and part5).
- Description of output files.
- Description of scripts.
- Description of example datasat.

.. toctree::
   :maxdepth: 2
   :caption: Documentation index

   description
   installation
   usage_exon_based
   usage_splice_site_based
   output_files_exon_based
   output_files_splice_site_based
