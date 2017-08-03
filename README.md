SPARpipe
========

A pipeline to extract results from raw SPAR-seq data.

This collection of scripts serves to help you get from raw sequencing output to tables with 
effect size per treatment. 

Few bells and whistles, and still under development. Comes strictly 'as-is' - use at your own risk.

Reference
---------
If you use this code, or parts of it, in published work, please cite:

Han H, Braunschweig U, Gonatopoulos-Pournatzis T, Weatheritt RJ, Hirsch CL, Ha KC, Radovani E, Nabeel-Shah S, Sterne-Weiler T, Wang J, O'Hanlon D, Pan Q, Ray D, Zheng H, Vizeacoumar F, Datti A, Magomedova L, Cummins CL, Hughes TR, Greenblatt JF, Wrana JL, Moffat J, Blencowe BJ. Multilayered Control of Alternative Splicing Regulatory Networks by Transcription Factors. Mol Cell. 2017 Feb 2;65(3):539-553.e7. [PubMed](https://www.ncbi.nlm.nih.gov/pubmed/28157508)

For feedback and questions, mail to u.braunschweig@utoronto.ca.

Input
-----
* FASTQ files with one read each for the fwd barcode (I1), rev barcode (I2),
  fwd event (R1), and rev event (R2). In a big project multiplexed over several
  Illumina lanes/runs, there will be 'batches' with one of each of the above.
* Fwd and rev barcode bowtie libraries that contain the expected reads, plus BED files for both 
  (junction_name 0 junction_length).
  These can be generated from FASTA junction files produced by the script `MakeJunctionsFASTA.R` 
  in the accessories folder. As input, it takes a CSV spreadsheet with the sequences of all the 
  elements in each amplicon and a string code indicating how they are connected, as well as a CSV
  file with primer sequences.
* A barcode table containing the expected combinations of fwd and rev barcodes (also produced by
  the accessory script)
* An event table created concomitantly with the junction library (also produced by the accessory
  script), detailing which junctions should be used to calculate the PSI for each (part of an) 
  event
* A treatment table specifying barcode numbers and batch along with treatment ID and replicate

Output
------
* Tables with raw and normalized PSI, dPSI, and SSMD
* Plots to follow normalization, check coverage, and monitor screen performance
* Demultiplexing and mapping stats
* Intermediate files per sample and batch such as FASTQ files of unmapped reads, 
  BAM files of mapped reads, raw PSI and count files

Dependencies
------------
* bowtie, samtools, R

Workflow
--------
0. Map fwd and rev barcode reads to barcode libraries with 
   `bowtie -v 2 -k 1 -m 1 --best --strata [--nofw|--norc] -S --sam-nohead <BC> <FASTQ> <BCx.sam>`
   Whether to use --nofw/norc may depend on the sequencing platform (HiSeq/NextSeq/MySeq)
   but it can probably be left out without much loss of specificity.
1. Demultiplex the fwd and rev event read files using the SAM files just created 
   and the script `1_demultiplex.pl`. If demultiplexing is done otherwise, make sure 
   that file names are compatible with downstream steps as sample numbers and batch IDs are
   taken from file names.
2. Map demultiplexed FASTQ files to junction libraries using `2_align.pl`.
3. Extract read counts and metrics from BAM, and generate files with raw PSI, RPM, read counts, 
   and pseudo-inclusion/exclusion reads using `3_count.pl`. In case the project is distributed 
   over several 'batches', they will be combined. Events and samples are checked against the 
   provided templates.
4. Perform normalization of raw PSI with `4_normalize.R` (currently, only plate normalization 
   and weighthed plate normalization are implemented) to get dPSI and SSMD scores.
   This step can also be performed on raw PSI data from other sources, as long as a suitable
   treatment table and a table with read counts per event are provided.
     
Known issues
------------
- There is no module yet to extract expression changes
