SPARpipe
========

**Tools to help design and extract results from raw SPAR-seq data**

Few bells and whistles, and still under development. Comes strictly 'as-is' - use at your own risk.

SPAR-seq is a strategy to obtain a sequencing-based readout for a selected set of RNA processing
events of interest from an arrayed screen using any kind of perturbation. Portions of endogenous
transcripts are simulatenously amplified in each well of a multi-well plate, barcoded to allow
tracking of the well of origin, pooled, and sequenced. 

This collection of scripts serves to help you get from raw sequencing output to tables with  
tables and plots showing the effect sizes of your treatments on inclusion or removal of parts
of transcripts. Also included are scripts to check primers for multiplex-RT-PCR for 3'-overlaps,
which are the major source of un-assignable reads; as well as for generating the junction library
for alignment of event reads and other required files for analysis.


Reference
---------
If you use this code, or parts of it, in published work, please cite:

Han H, Braunschweig U, Gonatopoulos-Pournatzis T, Weatheritt RJ, Hirsch CL, Ha KC, Radovani E, Nabeel-Shah S, Sterne-Weiler T, Wang J, O'Hanlon D, Pan Q, Ray D, Zheng H, Vizeacoumar F, Datti A, Magomedova L, Cummins CL, Hughes TR, Greenblatt JF, Wrana JL, Moffat J, Blencowe BJ. Multilayered Control of Alternative Splicing Regulatory Networks by Transcription Factors. Mol Cell. 2017 Feb 2;65(3):539-553.e7. [PubMed](https://www.ncbi.nlm.nih.gov/pubmed/28157508)

For feedback and questions, mail to u.braunschweig@utoronto.ca.


Dependencies
------------
* bowtie, samtools, R


Workflow
--------
### Setup
1. Design primers for your set of RNA processing events of interest. Avoid amplicons of
   dramatically different length to minimize amplification bias, and genes with dramatically 
   different expression lest you sequence mostly the few most highly expressed events.
2. Check primers using `CheckPrimers.R`.
3. Create a CSV file detailing the coordinates and sequences of the different parts of each 
   amplicon. Their properties and connections are provided by a string in column *structure*:
   * C1: upstream constitutive exon (part)
   * C2: dowmstream constitutive exon (part)
   * E1, E2, ...: alternative exons. Capital E for segments for which a score shall be produced.
   * [3] : Suffix to elements that are due to an alternative 3' splice site.
   * [5] : Same, for 5' splice site.
   * [i] : Suffix for intron retention event
   * `-` and `:` indicate connections via a splice and adjacent elements, respectively.

   Example: C1-E1:e2[5]:e3[5]-e4[3]:C2 indicates a situation with one cassette exon (E1) that
            has two nested alternative 5' splice sites (e2 and e3), and an alternative 3' splice 
            site at the downstream constitutive exon (e4). Only E1 will be monitored.

   Required columns are *gene*, *event*, *structure*, *chrom*, *strand*, *C1.start*, *C1.end*,
   *C2.start*, *C2.end*, and *X.start* and *X.end* columns for every alternative segment E1-En,
   (1-based coordinates) followed by sequences (5' to 3') *C1.seq*, *C2.seq* and *X.seq* for E1-En. 
4. Run `MakeJunctionsFASTA.R` to generate junction libraries in FASTA format, report the minimum
   edit distance between junctions (which should be as big as possible but > 1), and BED files 
   of the events required downstream.

### Experiment
Order your primers, conduct your experiment, sequence your amplicon pools with paired-end reads and
barcode reads (or see below).

### Analysis
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


Input for analysis
------------------
* Fwd and rev barcode bowtie libraries that contain the expected reads, plus BED files for both 
  (junction_name 0 junction_length).
  These can be generated from FASTA junction files produced by the script `MakeJunctionsFASTA.R` 
  in the accessories folder. As input, it takes a CSV spreadsheet with the sequences of all the 
  elements in each amplicon and a string code indicating how they are connected, as well as a CSV
  file with primer sequences.
* FASTQ files with one read each for the fwd barcode (I1), rev barcode (I2),
  fwd event (R1), and rev event (R2). In a big project multiplexed over several
  Illumina lanes/runs, there will be 'batches' with one of each of the above.
* A barcode table containing the expected combinations of fwd and rev barcodes
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


Known issues
------------
- There is no module yet to extract expression changes. See the paper for a strategy.
