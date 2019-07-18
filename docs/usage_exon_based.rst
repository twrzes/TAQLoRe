##########################################
Running the pipeline - exon-based approach
##########################################

As you may have read in :doc:`description`, the pipeline as a whole is comprised of five parts, three of which are used for exon-based analyses, and remaining two for splice-site-based method.

*********************
Snakemake config file
*********************

The first step of the pipeline is editing the Snakemake JSON config file ``config_TAQLoRe.json``. All the parameters and paths to input files are stored there. There are five main keys ('sections') in the configuration file.

workdir
=======

In this section there is only one key - ```workdir``, denoting the path to the working directory.

samples
=======

This section of config file contains two keys:

- ``fasta_file_dir``: Directory where FASTA files are stored.

.. warning::
  All the FASTA files must have consistent naming in the format of {run_prefix}.{barcode}.fa,
  where {run_prefix} is a prefix of the run (i.e. date of sequencing), and {barcode} is a barcode
  number. This must also be consistent with barcode-to-sample mapping file, as described in :ref:`mappings`.

- ``downsampled_reads_num``: Number of reads to downsample to in third part of the pipeline - :ref:`part3`.

input_files
===========

This section contains paths to all necessary input files:

- ``transcriptome_fasta``: Path to the FASTA file with all transcript sequences.

.. warning::
  All transcript names (i.e. headers of each FASTA sequence) **MUST NOT** contain any spaces, as they cannot be parsed correctly.

  Moreover, FASTA file should be single-line FASTA file (i.e. each sequence must contain two lines - first being FASTA header, and the second one being a sequence).

- ``last_transcriptome_index``: Prefix of LAST index for a transcriptome.
- ``last_transcriptome_index_dummy_file``: Path to the dummy file to be created after generating LAST index for a transcriptome (required to maintain the order of steps in the pipeline).
- ``genome_fasta``: Path to the file with genome FASTA sequence.

.. warning::
  All chromosome names must be consistent between this file and both GTF file and sections of Snakemake config file. The file **MUST NOT** contain any spaces and be single-line FASTA file.

- ``chrom_sizes``: Path to the file containing chromosome sizes for a genome.
- ``last_genome_index``: Prefix of the LAST index for a genome.
- ``last_genome_index_dummy_file``: Path to the dummy file being created after generating LAST index for a genome (required to maintain the order of steps in the pipeline).
- ``gtf_file``: Path to the GTF file with gene annotations.

.. warning::
  All chromosome names must be consistent between this file and both genome FASTA file and sections of Snakemake config file.

- ``all_exons_positions_ENSEMBL``: Path to the file containing gene annotations from ENSEMBL BioMart. You can find how to create this file in this section: :ref:`bioMart_input`.
- ``barcode_to_sample_file``: Path to the file with barcode-to-sample mappings. The structure of the file can be found in this section: :ref:`mappings`.
- ``barcode_to_sample_file_downsampled``: Path to the file with barcode-to-sample mappings for downsampling. The structure of the file can be found in this section: :ref:`mappings`. More about downsampling can be found in :ref:`part3`.

.. _bioMart_input:

How to create an annotation file with ENSEMBL bioMart
-----------------------------------------------------

To create a file needed in the pipeline, you need to use ENSEMBL BioMart. The reason behind it is that some information required to run the pipeline may not be included in the GTF file (e.g. UTRs).

The file can be created with following steps:

1. Go to http://www.ensembl.org/biomart/martview
2. From **ENSEMBL Genes** database choose the dataset of interest (e.g. **Human Genes**).
3. From **Filters** section on the left-hand side of the website choose the gene of interest (e.g. select **Input external references ID list [Max 500 advised]**, and put ENSEMBL gene ID into the field).
4. From **Attributes** section on the left-hand side of the website select **Structures**, deselect all attributes (from **GENE** subsection), and select following attributes (**IN THIS PARTICULAR ORDER**) (subsections of attributes are written in square brackets):

- Gene stable ID [GENE]
- Gene start (bp) [GENE]
- Gene end (bp) [GENE]
- Transcript stable ID [GENE]
- Transcript start (bp) [GENE]
- Transcript end (bp) [GENE]
- Transcription start site (TSS) [GENE]
- Exon stable ID [EXON]
- Exon region start (bp) [EXON]
- Exon region end (bp) [EXON]
- Exon rank in transcript [EXON]
- cDNA coding start [EXON]
- cDNA coding end [EXON]
- Genomic coding start [EXON]
- Genomic coding end [EXON]

5. Click **Results** button on top-left side of the website.
6. From **Export  all results to ** section, choose **File** and **TSV** format, then click **Go** and save the file in the desired location.

.. _mappings:

How to create a barcode-to-sample mapping file
----------------------------------------------

The file with metadata is a tab-delimited file (without header) containing three columns:

+--------------+-----------+--------------------+
| Run prefix   | Barcode   |    Sample name     |
+==============+===========+====================+
| 2017_01_13   | barcode01 | Jan_5238_cingulate |
+--------------+-----------+--------------------+
| 2017_06_15   | barcode12 | Jun_5346_striatum  |
+--------------+-----------+--------------------+

where:

- ``Run prefix`` - is the run prefix for each file. This can be e.g. a date of the sequencing or any string that denotes different sequencing batches.
- ``Barcode`` - is the barcode of each file.
- ``Sample name`` - is the underscore-separated string with a following structure: ``{run_name}_{sample_id}.{sample_sub_id}```. In the example above ``{run_name}`` denotes two sequencing runs (one in January, second one in June), ``{sample_id}`` denotes different individuals, and ``{sample_sub_id}`` denotes different brain regions.

gene_info
=========

This section of config file contains the information about analysed gene.

- ``gene_name``: Gene name for the gene of interest.
- ``chromosome_name``: Chromosome name for the gene of interest.

.. warning::
  All chromosome names must be consistent between this section and both genome FASTA file and GTF file.

- ``gene_start``: Start position of the gene (0-based coordinates).
- ``gene_end``: End position of the gene (0-based coordinates).
- ``strand``: Strand of the gene ('+' or '-').

parameters
==========

This section of config file contains all the parameters being used by scripts.

- ``min_prop``: Minimum proportion of read covering a transcript.
- ``min_prop_align``: Minimum proportion of transcript being covered by aligned read.
- ``min_insert``: Minimum length of potential novel exon (insertion in the alignment).
- ``min_exon_distance``: Minimum distance from annotated exon.
- ``distance_between_exons``: Minimum distance between exons (both known and novel).
- ``exon_coverage_threshold``: Minimum coverage of exon for exon to be included in transcripts.
- ``min_exon_length``: Minimum length of the exon.
- ``min_reads_threshold``: Minimum number of reads covering the exon for the exon to be included in transcripts.
- ``min_num_individuals_threshold``: Minimum number of different individuals (``{sample_id}``) having at least ``min_reads_threshold`` reads per exon.
- ``min_num_libraries_threshold``: Minimum number of different tissues (``{sample_sub_id}``) having at least ``min_reads_threshold`` reads per exon.
- ``sum_threshold``: Minimum sum of reads for a transcript to be included in the annotation.
- ``reads_in_sample_threshold``: Minimum number of reads per sample for a transcript to be included in the annotation.
- ``sample_threshold``: Minimum number of samples having at least ``reads_in_sample_threshold`` reads for a transcript to be included in the annotation.

********************
Running the pipeline
********************

Local computer/server
=====================

To run each part of the pipeline on a computer/server, you can run it by typing:

.. prompt:: bash $

  cd /path/to/TAQLoRe
  conda activate taqlore
  snakemake -j {number_of_cores} -s TAQLore_part_1
  snakemake -j {number_of_cores} -s TAQLore_part_2

etc.

where ``{number_of_cores}`` is the number of cores to use for a Snakemake run.

.. note::
  LAST alignements are very memory-intensive - for 100k reads LAST uses ~48G of memory (human genome/transcriptome). Therefore, the best way to run the pipeline is to use HPC (see below).

.. note::
  The way how snakemake is run in the example above requires a user to pre-install the whole environment located in `envs/taqlore.yaml`. To pre-install this environment please refer to :doc:`installation`. If a user wants to create conda environment from scratch, they can use `snakemake -j {number_of_cores} -s TAQLore_part_1 --use-conda` command which will create a local copy of the whole environment in the working directory.

High Performance Computing
==========================

In order to run Snakemake pipeline on a computational cluster (preferred way), two additional files must be created

.. _cluster_config:

Cluster configuration file
--------------------------

For additional information, refer to `Snakemake documentation <https://snakemake.readthedocs.io/en/stable/>`_.

This file contains all parameters for each jobs to be used by the scheduler, such as time, memory, number of CPUs, partition name, etc.

The example JSON file to run using SLURM scheduler (Earlham Institute's infrastructure) looks like this:

.. code-block:: json

  {
    "__default__" :
      {
          "nodes" : 1,
          "time" : "7-00:00:00",
          "n" : 1,
          "ntasks-per-node": 1,
          "cpu" : 1,
          "partition" : "medium",
          "memory" : "64G",
          "job_name" : "{rule}.{wildcards}",
          "out" : "slurm.%N.%j.{rule}.{wildcards}.out",
          "err" : "slurm.%N.%j.{rule}.{wildcards}.err"
      },
      "last_index_transcriptome" :
      {
          "time" : "1-00:00:00",
          "cpu" : 8,
          "memory" : "64G"
      },
      "last_train_gap_mismatch_transcriptome" :
      {
          "time" : "1-00:00:00",
          "cpu" : 8,
          "memory" : "64G"
      },
      "last_align_transcriptome" :
      {
          "time" : "7-00:00:00",
          "cpu" : 8,
          "memory" : "64G"
      },
      "last_index_genome" :
      {
          "time" : "1-00:00:00",
          "cpu" : 8,
          "memory" : "64G"
      },
      "novel_exons_alignment_to_genome_parsing_last_maf" :
      {
          "time" : "00:45:00",
          "cpu" : 8,
          "partition" : "short",
          "memory" : "32G"
      },
      "novel_exons_genomic_coordinates":
      {
          "time" : "00:45:00",
          "partition" : "short",
          "memory" : "16G"
      },
      "filtering_genomic_positions_gene_boundaries" :
      {
          "time" : "00:05:00",
          "partition" : "short",
          "memory" : "4G"
      },
      "novel_exons_per_library" :
      {
          "time" : "00:45:00",
          "partition" : "short",
          "memory" : "4G"
      },
      "novel_exons_summary" :
      {
          "time" : "00:45:00",
          "partition" : "short",
          "memory" : "4G"
      },
      "generate_BedGraph_sum" :
      {
          "time" : "00:45:00",
          "partition" : "short",
          "memory" : "4G"
      },
      "novel_exons_file_1nt_coordinates" :
      {
          "time" : "00:01:00",
          "partition" : "short",
          "memory" : "2G"
      },
      "coordinates_genomic_meta_gene_exons" :
      {
          "time" : "00:05:00",
          "partition" : "short",
          "memory" : "16G"
      }
  }

.. note::
  ``__default__`` section specifies parameters of default job, while parameters under rule names (copied from Snakemake files, e.g. ``last_index_transcriptome``) denote deviation(s) from default rule (e.g. different time, memory, number of CPUs, etc.).

.. _batch_script_to_submit:

Batch script to submit
----------------------

In order to run each part of the pipeline using HPC, a batch script to submit must be created. An example of the batch script (SLURM scheduler, Earlham Institute's infrastructure) can be seen below:

.. code-block:: bash

  #!/bin/bash
  #SBATCH -p medium
  #SBATCH -N 1
  #SBATCH -n 1
  #SBATCH -c 1
  #SBATCH --mem 4G
  #SBATCH -t 7-00:00:00
  #SBATCH -o slurm.%N.%j.out
  #SBATCH -e slurm.%N.%j.err
  #SBATCH --mail-type=ALL
  #SBATCH --mail-user=some.user@some.insitute.ac.uk

  module load conda

  conda activate taqlore

  srun snakemake -s TAQLoRe_part1 --latency-wait 60 -j {number_of_jobs} --cluster-config /path/to/cluster/config.json --cluster "sbatch -p {cluster.partition} -N {cluster.nodes} -n {cluster.n} --ntasks-per-node={cluster.ntasks-per-node} -c {cluster.cpu} -t {cluster.time} --mem {cluster.memory} -J {cluster.job_name} -o slurm.%N.%j.out -e slurm.%N.%j.err --mail-type=FAIL --mail-user=some.user@some.insitute.ac.uk"

where:

- ``--cluster-config`` denotes path to the :ref:`cluster_config`.
- ``{number_of_jobs}`` denotes number of jobs submitted to the cluster at once.

The script can be submitted with command (SLURM scheduler):

.. prompt:: bash $

  sbatch batch_script_to_submit.sh

where ```batch_script_to_submit.sh`` is the name of the file which contents are shown above (:ref:`batch_script_to_submit`).

.. _after_part1:

*************************************************
Things to do after running part 1 of the pipeline
*************************************************

The first part of the pipeline ends with creating a meta-gene annotation file in ```{workdir}/results/meta_gene_construction/meta_gene_genomic_exon_coordinates.txt``.

In order to run the second part of the pipeline, the user must annotate UTRs in the ``meta_gene_genomic_exon_coordinates.txt`` file and/or add/remove additional/unnecessary exons. The reason behind it is that there are different sources of UTRs (e.g. ENSEMBL, GENCODE, RefSeq, UCSC, etc.) and the differences between annotations can be significant.

UTRs must be added as the last column to the ``meta_gene_genomic_exon_coordinates.txt`` file, for each exon position in the meta-gene. The last column should have 'UTR' string denoting that a position is an UTR, and any other string (e.g. 'NA') denoting the non-UTR/unknown status of a genomic region in a gene. Example:

- ``meta_gene_genomic_exon_coordinates.txt`` after part1 of the pipeline::

    1	277	1970786	1971062	ENST00000543114	ENSE00001774617	No
    278	416	1971063	1971201	ENST00000543114	ENSE00001774617	Yes
    417	656	2053298	2053537	ENST00000335762;ENST00000399655	ENSE00001539923;ENSE00001539923	No;No
    657	681	2053538	2053562	ENST00000335762;ENST00000399655;ENST00000480911	ENSE00001539923;ENSE00001539923;ENSE00001839973	No;No;No
    682	730	2053563	2053611	ENST00000335762;ENST00000399655;ENST00000480911;ENST00000399595;ENST00000399644;ENST00000399638;ENST00000399597;ENST00000399621;ENST00000399637;ENST00000399591;ENST00000399641;ENST00000347598;ENST00000399606;ENST00000399601;ENST00000344100;ENST00000399629;ENST00000327702;ENST00000399649;ENST00000402845;ENST00000399603;ENST00000399634;ENST00000399617;ENST00000406454	ENSE00001539923;ENSE00001539923;ENSE00001839973;ENSE00001539466;ENSE00001539466;ENSE00001539466;ENSE00001539466;ENSE00001539466;ENSE00001539466;ENSE00001539466;ENSE00001539466;ENSE00001539466;ENSE00001539466;ENSE00001539466;ENSE00001539466;ENSE00001539466;ENSE00001539466;ENSE00001539466;ENSE00001539466;ENSE00001539466;ENSE00001539466;ENSE00001539466;ENSE00001539466	Yes;Yes;Yes;Yes;Yes;Yes;Yes;Yes;Yes;Yes;Yes;Yes;Yes;Yes;Yes;Yes;Yes;Yes;Yes;Yes;Yes;Yes;Yes
    ...
    23196	23495	2690900	2691199 ENST00000335762;ENST00000399655;ENST00000399595;ENST00000399644;ENST00000399638;ENST00000399597;ENST00000399621;ENST00000399637;ENST00000399591;ENST00000399641;ENST00000347598;ENST00000399606;ENST00000399601;ENST00000344100;ENST0000039962 ;ENST00000327702;ENST00000399649;ENST00000402845;ENST00000399603;ENST00000399634;ENST00000399617;ENST00000406454;ENST00000616390	ENSE00002228600;ENSE00001539473;ENSE00001724521;ENSE00001724521;ENSE00001724521;ENSE00001724521;ENSE00001724521;ENSE00001724521;ENSE00001724521;ENSE00001724521;ENSE00001724521;ENSE00001724521;ENSE00001724521;ENSE00001724521;ENSE00001724521;ENSE00001724521;ENSE00001724521;ENSE00001724521;ENSE00001539391;ENSE00001539391;ENSE00001539391;ENSE00001539391;ENSE00003738703	Yes;Yes;Yes;Yes;Yes;Yes;Yes;Yes;Yes;Yes;Yes;Yes;Yes;Yes;Yes;Yes;Yes;Yes;Yes;Yes;Yes;Yes;Yes
    23496	23563	2691200	2691267	ENST00000335762;ENST00000399655;ENST00000399595;ENST00000399644;ENST00000399638;ENST00000399597;ENST00000399621;ENST00000399637;ENST00000399591;ENST00000399641;ENST00000347598;ENST00000399606;ENST00000399601;ENST00000344100;ENST00000399629;ENST00000327702;ENST00000399649;ENST00000402845;ENST00000399603;ENST00000399634;ENST00000399617;ENST00000406454;ENST00000616390	ENSE00002228600;ENSE00001539473;ENSE00001724521;ENSE00001724521;ENSE00001724521;ENSE00001724521;ENSE00001724521;ENSE00001724521;ENSE00001724521;ENSE00001724521;ENSE00001724521;ENSE00001724521;ENSE00001724521;ENSE00001724521;ENSE00001724521;ENSE00001724521;ENSE00001724521;ENSE00001724521;ENSE00001539391;ENSE00001539391;ENSE00001539391;ENSE00001539391;ENSE00003738703	No;No;No;No;No;No;No;No;No;No;No;No;No;No;No;No;No;No;No;No;No;No;No
    23564	23589	2691268	2691293	ENST00000335762;ENST00000399655;ENST00000399595;ENST00000399644;ENST00000399638;ENST00000399597;ENST00000399621;ENST00000399637;ENST00000399591;ENST00000399641;ENST00000347598;ENST00000399606;ENST00000399601;ENST00000344100;ENST00000399629;ENST00000327702;ENST00000399649;ENST00000402845;ENST00000399603;ENST00000399634;ENST00000399617;ENST00000406454	ENSE00002228600;ENSE00001539473;ENSE00001724521;ENSE00001724521;ENSE00001724521;ENSE00001724521;ENSE00001724521;ENSE00001724521;ENSE00001724521;ENSE00001724521;ENSE00001724521;ENSE00001724521;ENSE00001724521;ENSE00001724521;ENSE00001724521;ENSE00001724521;ENSE00001724521;ENSE00001724521;ENSE00001539391;ENSE00001539391;ENSE00001539391;ENSE00001539391	No;No;No;No;No;No;No;No;No;No;No;No;No;No;No;No;No;No;No;No;No;No
    23590	23706	2691294	2691410	ENST00000335762;ENST00000399655;ENST00000399603;ENST00000399634;ENST00000399617;ENST00000406454	ENSE00002228600;ENSE00001539473;ENSE00001539391;ENSE00001539391;ENSE00001539391;ENSE00001539391	No;No;No;No;No;No
    23707	24455	2691411	2692159	ENST00000399655;ENST00000399603;ENST00000399634;ENST00000399617;ENST00000406454	ENSE00001539473;ENSE00001539391;ENSE00001539391;ENSE00001539391;ENSE00001539391	No;No;No;No;No
    24456	30246	2692160	2697950	ENST00000399655	ENSE00001539473	No

- ``meta_gene_genomic_exon_coordinates.txt`` after adding UTR annotations (before running part2 of the pipeline)::

    1	277	1970786	1971062	ENST00000543114	ENSE00001774617	No	UTR
    278	416	1971063	1971201	ENST00000543114	ENSE00001774617	Yes	NA
    417	656	2053298	2053537	ENST00000335762;ENST00000399655	ENSE00001539923;ENSE00001539923	No;No	UTR
    657	681	2053538	2053562	ENST00000335762;ENST00000399655;ENST00000480911	ENSE00001539923;ENSE00001539923;ENSE00001839973	No;No;No	UTR
    682	730	2053563	2053611	ENST00000335762;ENST00000399655;ENST00000480911;ENST00000399595;ENST00000399644;ENST00000399638;ENST00000399597;ENST00000399621;ENST00000399637;ENST00000399591;ENST00000399641;ENST00000347598;ENST00000399606;ENST00000399601;ENST00000344100;ENST00000399629;ENST00000327702;ENST00000399649;ENST00000402845;ENST00000399603;ENST00000399634;ENST00000399617;ENST00000406454	ENSE00001539923;ENSE00001539923;ENSE00001839973;ENSE00001539466;ENSE00001539466;ENSE00001539466;ENSE00001539466;ENSE00001539466;ENSE00001539466;ENSE00001539466;ENSE00001539466;ENSE00001539466;ENSE00001539466;ENSE00001539466;ENSE00001539466;ENSE00001539466;ENSE00001539466;ENSE00001539466;ENSE00001539466;ENSE00001539466;ENSE00001539466;ENSE00001539466;ENSE00001539466	Yes;Yes;Yes;Yes;Yes;Yes;Yes;Yes;Yes;Yes;Yes;Yes;Yes;Yes;Yes;Yes;Yes;Yes;Yes;Yes;Yes;Yes;Yes	NA
    ...
    23196	23495	2690900	2691199	ENST00000335762;ENST00000399655;ENST00000399595;ENST00000399644;ENST00000399638;ENST00000399597;ENST00000399621;ENST00000399637;ENST00000399591;ENST00000399641;ENST00000347598;ENST00000399606;ENST00000399601;ENST00000344100;ENST00000399629;ENST00000327702;ENST00000399649;ENST00000402845;ENST00000399603;ENST00000399634;ENST00000399617;ENST00000406454;ENST00000616390	ENSE00002228600;ENSE00001539473;ENSE00001724521;ENSE00001724521;ENSE00001724521;ENSE00001724521;ENSE00001724521;ENSE00001724521;ENSE00001724521;ENSE00001724521;ENSE00001724521;ENSE00001724521;ENSE00001724521;ENSE00001724521;ENSE00001724521;ENSE00001724521;ENSE00001724521;ENSE00001724521;ENSE00001539391;ENSE00001539391;ENSE00001539391;ENSE00001539391;ENSE00003738703	Yes;Yes;Yes;Yes;Yes;Yes;Yes;Yes;Yes;Yes;Yes;Yes;Yes;Yes;Yes;Yes;Yes;Yes;Yes;Yes;Yes;Yes;Yes	NA
    23496	23563	2691200	2691267	ENST00000335762;ENST00000399655;ENST00000399595;ENST00000399644;ENST00000399638;ENST00000399597;ENST00000399621;ENST00000399637;ENST00000399591;ENST00000399641;ENST00000347598;ENST00000399606;ENST00000399601;ENST00000344100;ENST00000399629;ENST00000327702;ENST00000399649;ENST00000402845;ENST00000399603;ENST00000399634;ENST00000399617;ENST00000406454;ENST00000616390	ENSE00002228600;ENSE00001539473;ENSE00001724521;ENSE00001724521;ENSE00001724521;ENSE00001724521;ENSE00001724521;ENSE00001724521;ENSE00001724521;ENSE00001724521;ENSE00001724521;ENSE00001724521;ENSE00001724521;ENSE00001724521;ENSE00001724521;ENSE00001724521;ENSE00001724521;ENSE00001724521;ENSE00001539391;ENSE00001539391;ENSE00001539391;ENSE00001539391;ENSE00003738703	No;No;No;No;No;No;No;No;No;No;No;No;No;No;No;No;No;No;No;No;No;No;No	UTR
    23564	23589	2691268	2691293	ENST00000335762;ENST00000399655;ENST00000399595;ENST00000399644;ENST00000399638;ENST00000399597;ENST00000399621;ENST00000399637;ENST00000399591;ENST00000399641;ENST00000347598;ENST00000399606;ENST00000399601;ENST00000344100;ENST00000399629;ENST00000327702;ENST00000399649;ENST00000402845;ENST00000399603;ENST00000399634;ENST00000399617;ENST00000406454	ENSE00002228600;ENSE00001539473;ENSE00001724521;ENSE00001724521;ENSE00001724521;ENSE00001724521;ENSE00001724521;ENSE00001724521;ENSE00001724521;ENSE00001724521;ENSE00001724521;ENSE00001724521;ENSE00001724521;ENSE00001724521;ENSE00001724521;ENSE00001724521;ENSE00001724521;ENSE00001724521;ENSE00001539391;ENSE00001539391;ENSE00001539391;ENSE00001539391	No;No;No;No;No;No;No;No;No;No;No;No;No;No;No;No;No;No;No;No;No;No	UTR
    23590	23706	2691294	2691410	ENST00000335762;ENST00000399655;ENST00000399603;ENST00000399634;ENST00000399617;ENST00000406454	ENSE00002228600;ENSE00001539473;ENSE00001539391;ENSE00001539391;ENSE00001539391;ENSE00001539391	No;No;No;No;No;No	UTR
    23707	24455	2691411	2692159	ENST00000399655;ENST00000399603;ENST00000399634;ENST00000399617;ENST00000406454	ENSE00001539473;ENSE00001539391;ENSE00001539391;ENSE00001539391;ENSE00001539391	No;No;No;No;No	UTR
    24456	30246	2692160	2697950	ENST00000399655	ENSE00001539473	No	UTR

After adding UTR annotation to the ``meta_gene_genomic_exon_coordinates.txt`` file, the pipeline can be run as usual (file: ``TAQLoRe_part2``).

.. _part3:

*******************************
Read counts bias - downsampling
*******************************

Third part of the pipeline () can be run in order to remove read counts bias (i.e. different number of reads in analysed samples). This step is necessary if you have huge differences in read numbers (> 5000 reads difference between samples with highest and lowest number of reads).

Before running part3 of the pipeline
====================================

Before running part3 of the pipeline two additional steps are necessary:

- Before running the pipeline ``config_TAQLoRe.json`` needs to be edited (in the section ``"samples"``, sub-section ``downsampled_reads_num``), to reflect the number of reads to choose for all the files. The downsampling will be done on ``/path/to/workdir/results/meta_gene_exon_counts_splicing_patterns/{run_prefix}.{barcode}_splicing_patterns_cds.tmp`` files, where ``{run_prefix}`` and ``{barcode}`` first and second column from :ref:`mappings` file, respectively. Thus, to see the number of reads in each file (sorted by the number of reads) the following command may be invoked to see the number of reads in each sample:

.. prompt:: bash $

  cd /path/to/workdir
  for i in `ls results/meta_gene_exon_counts_splicing_patterns/*_splicing_patterns_cds.tmp`; do j=`cat $i | wc -l`; printf "${i}\t${j}\n"; j=''; done | sort -k2,2n

- A meta-data file needs to be edited, as removing some outliers with the lowest numbers of reads might be necessary to accurately compare the expression between samples. The file has the same structure as the one in :ref:`mappings`. The path to this file should be put in ``config_TAQLoRe.json``, section ``input_files``, sub-section ``barcode_to_sample_file_downsampled``.

After these steps the pipeline can be run as usual (file: ``TAQLore_part3``).
