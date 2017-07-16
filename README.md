SPARpipe
========

A pipeline to extract results from raw SPAR-seq data.

No bells and whistles. Comes strictly 'as-is' - use at your own risk.

Input
-----
* FASTQ files with one read each for the fwd barcode (I1), rev barcode (I2),
  fwd event (R1), and rev event (R2). In a big project multiplexed over several
  Illumina lanes/runs, there will be 'batches' with one of each of the above.
* Fwd and rev barcode bowtie libraries
* A barcode table containing the expected combinations of fwd and rev barcodes.
* Fwd and rev juncion libraries in bowtie format that contain the expected reads.
* An event table created concomitantly with the junction library, detailing which
  junctions to use for which (part of an) event
* An treatment table [### pending implementation ###]

Output
------
* Lots of intermediate files per sample such as FASTQ files of unmapped reads, 
  BAM files of mapped reads
* Per 'batch', files with PSI, reads per event, inclusion/exclusion reads etc.
* Demultiplexing and mapping stats

Dependencies
------------
* bowtie, samtools, R

Workflow
--------
1. Map fwd and rev barcode reads to barcode libraries with 
     'bowtie -v 2 -k 1 -m 1 --best --strata [--nofw|--norc] -S --sam-nohead <BC> <FASTQ> <BCx.sam>'
2. Demultiplex the fwd and rev event read files using the SAM files just created 
     and the script 1_demultiplex.pl. If demultiplexing is done otherwise, make sure 
     that file names are compatible with downstream steps as sample numbers and batch IDs are
     taken from file names.
3. Map demultiplexed FASTQ files to junction libraries using 2_align.pl.
4. Extract read counts and metrics from BAM files using 3_count.pl. In case the
     project is distributed over several 'batches', each one is run separately.
5. Combine batches  ###PENDING###
     