SPARpipe
========

**Tools to help design a SPAR-seq screen and extract results from raw sequencing data**

**_NOTE: If you are looking for the application to detection of SAR!S-CoV-2 (Aynaud, Hernandez, Barucu et al., Nature Communications 2021), please checkout the 'COVID-19' branch_**

Few bells and whistles, and still under development. Comes strictly 'as-is' - use at your own risk.

SPAR-seq is a strategy to obtain a sequencing-based readout for a selected set of RNA processing
events of interest from an arrayed screen using any kind of perturbation. Portions of endogenous
transcripts are simulatenously amplified in each well of a multi-well plate, barcoded to allow
tracking of the well of origin, pooled, and sequenced. 

This collection of scripts serves to help you develop a SPAR-seq screen and extract tables with 
raw and normalized data from sequencing output. It will also produce plots documenting mapping 
and normalization and showing the effect sizes of your treatments on inclusion or removal of parts 
of transcripts. 
Also included are scripts to check primers for multiplex-RT-PCR for 3'-overlaps,
which are the major source of un-assignable reads; as well as for generating the junction library
for alignment of event reads and other required files for analysis.


Dependencies
------------
* SPARpipe is only available for Linux and has not been tested on MacOS.
* bowtie >= 1.2.1, samtools, R, pigz (available in your PATH)


Workflow
--------
### Setup
1. **Design primers** for your set of RNA processing events of interest using your favourite software, taking into consideration:
   * Avoid amplicons of dramatically different lengths to minimize amplification bias. For intron retention events (or some other cases), that can make it necessary to use a 3-primer stragegy with two primers annealing to the flanking exon, and an intron-internal primer to pick up intron-retained transcripts. This has worked well for us; more complex strategies with more primers are possible, but make sure that there is no cross-amplification which would result in amplicon abundances that do not reflect the isoform distribution in the cell.
   * Avoid targeting genes with dramatically different expression lest you sequence mostly the few most highly expressed events.
   * Different universal adaptors (Illumina) are added to the 5'-ends of the forward and reverse primers. These are then used to add multiplexing barcodes in a second round of PCR - see publications below.
   * Check primers for 3'-overlaps using the `CheckPrimers.R` accessory. We have found that overlaps of 5 bp 
   and more may cause primer dimer problems. Not all 5-bp overlaps result in dimers, but all dimers have them. Longer overlaps are worse.
   * Make sure that amplicons are sufficiently different from one another; using the script mentioned in step 3 will check for that. Usually, designing the primers close to the monitored splice junctions maximizes the differential sequence space that can be monitored.

2. **Test each primer set** on its own to make sure it amplifies the correct bands. Then test them all in a pool - when using several dozen primer sets, you should see a smear with few strong bands sticking out.

3. Run `MakeJunctionsFASTA.R` to **generate junction libraries** in FASTA format, report the minimum
   edit distance between junctions (which should be as big as possible but > 1), and BED files
   of the events required downstream by 3_combine.pl. This script also checks if the elements of
   the RNA events that you are interested in can be reliably measured given a read length.
   

   The **eventFile** is a CSV file detailing the coordinates and sequences of the different parts of 
   each amplicon. Required columns are *gene*, *event*, *structure*, *chrom*, *strand*, *C1.start*,
   *C1.end*, *C2.start*, *C2.end*, and *X.start* and *X.end* columns for every alternative segment 
   E1-En, (1-based coordinates) followed by sequences (5' to 3') *C1.seq*, *C2.seq* and *X.seq* for 
   E1-En. 

   The column *structure* contains a string that described the properties of events and how the
   elements making up the event are connected:
   * C1: upstream constitutive exon (part)
   * C2: dowmstream constitutive exon (part)
   * E1, E2, ...: alternative exons. Capital E for segments for which a score shall be produced.
   * [3] : Suffix to elements that are due to an alternative 3' splice site.
   * [5] : Same, for 5' splice site.
   * [c] : Internal exon is constitutive. If a primer anneals across exons, designate the most distal
           as C1 or C2 and the other as an intervening element with suffix [c].
   * [i] : Suffix for intron retention event
   * `-` and `:` indicate connections via a splice and adjacent elements, respectively.

   Example: C1-E1:e2[5]:e3[5]-e4[3]:C2 indicates a situation with one cassette exon (E1) that
           has two nested alternative 5' splice sites (e2 and e3), and an alternative 3' splice 
           site at the downstream constitutive exon (e4). Only E1 will be monitored.

   
   The **primerFile** is a CSV file with columns *event*, *seqF* and *seqR* (5' to 3', including 
   adaptors)

4. Use the forward and reverse junction FASTA files to **generate (separate) bowtie indices**.

5. Generate a **treatment table** for use downstream. It must have these columns:
   * *ID*, a unique number starting with 'T' (e.g. T001) for each individual treatment being tracked. 
         This can include similar treatments done multiple times.
   * *Replicate*
   * *Batch*: If a large experiment requires re-using of Barcodes, designate data coming from one batch
       (e.g., Illumina lane) with a unique string (cannot contain '_').
   * *Barcode*: Unique number starting with 'W' (e.g. W001) idenfifying a unique combination of fwd and
       rev barcodes.
   * *Treatment*: short string (no '_' characters) describing the treatment, e.g. 'DMSO' or 'siSrrm4'.
   * *Type*: One of 'posCtl' (positive control), 'negCtl' (negative control), 'ctl' (other control),
       or 'exp' (experimental). Type designations are used during normalization and for plots.
   * *Plate*: Number of physical multi-well plate a given sample was derived from. Used for normalization.
   

### Experiment
Conduct your experiment, sequence your amplicon pools with paired-end reads and
barcode reads (or see below). We strongly recommend doing a low-scale test sequencing run to allow recognizing and removing primer sets that do not amplify well, cause excessive dimers, or take up too much sequence space. If a large experiment requires re-use of
barcode combinations, that is not a problem as long as each batch of barcodes is kept separate.


### Analysis
0. Map fwd and rev barcode reads to barcode libraries with 
   `bowtie -v 2 -k 1 -m 1 --best --strata [--nofw|--norc] -S --sam-nohead <BC> <FASTQ> <BCx.sam>`
   Whether to use --nofw/norc may depend on the sequencing platform (HiSeq/NextSeq/MiSeq...)
   but it can probably be left out without much loss of specificity. Run on only 1 thread to
   guarantee read order is maintained. Not necessary if you receive sequencing data demultiplexed.
1. Demultiplex the fwd and rev event read files using the SAM files just created 
   and the script `1_demultiplex.pl`. If demultiplexing is done otherwise, make sure 
   that file names are compatible with downstream steps as sample numbers and batch IDs are
   taken from file names. File names must conform to the pattern '[BATCH]_W[NUMBER]_[OPTIONAL]_[Fwd/Rev].fq(.gz)',
   where BATCH is a string that does not contain '_', and NUMBER is a unique number identifying the barcode
   combination.
2. Map demultiplexed FASTQ files to junction libraries using `2_align.pl`. This is done separately
   for fwd and rev reads.
3. Extract read counts and metrics from BAM, and generate files with raw PSI, RPM, read counts, 
   and pseudo-inclusion/exclusion reads using `3_count.pl`. In case the project is distributed 
   over several 'batches' (e.g. re-using barcodes and run in different lanes), they will be combined.
   Events and samples are checked against the provided templates specifying the design. Additionally,
   mapping stats per sample are documented.
4. Perform normalization of raw PSI with `4_normalize.R` to get dPSI and SSMD scores. Currently,
   several versions are supported that each try to mitigate plate-derived batch effects, using
   either the negative controls (desireable, but requires a reasonable number on each plate),
   all wells on the plate (assuming random placement of experimental treatments), or intermediates.
   This step can also be performed on raw PSI data from other sources, as long as a suitable
   treatment table and a table with read counts per event are provided.


Input for analysis
------------------
* Fwd and rev barcode bowtie indexes that contain the expected reads, plus BED files for both 
  (junction_name 0 junction_length).
  We strongly recommend using the `MakeJunctionsFASTA.R` accessory script to produce FASTA files to 
  make bowtie indexes from, see 'Setup' above.
* FASTQ files with one read each for the fwd barcode (i5), rev barcode (i7),
  fwd event (R1) read, and rev event (R2) read. In a big project multiplexed over several
  Illumina lanes/runs, there will be 'batches' with one of each of the above.
  Barcode combinations must be unique within a batch but can be re-used between batches.
  Alternatively, if samples are demultiplexed already, the file names need to conform to SPARpipe standard
  and barcode read files are not necessary (see above).
* A barcode table containing the expected combinations of fwd and rev barcodes (if not de-multiplexed).
* The event table created concomitantly with the junction library (also produced by the accessory
  script `MakeJunctionsFASTA.R`). See above.
* A treatment table specifying barcode numbers and batch along with treatment ID and replicate.
  See above.


Output
------
* Tables with raw and normalized PSI, dPSI, and SSMD
* Demultiplexing and mapping stats
* Plots to follow normalization, check coverage, and monitor screen performance
* Intermediate files per sample and batch such as FASTQ files of unmapped reads, 
  BAM files of mapped reads, raw PSI and count files


Known issues
------------
- There is no module yet to extract expression changes. See the 2017 paper for a strategy.


Citation
--------
If you use this code, or parts of it, in published work, please cite:

Systematic exploration of dynamic splicing networks reveals conserved multistage regulators of neurogenesis. Han H, Best AJ, Braunschweig U, Mikolajewicz N, Li JD, Roth J, Chowdhury F, Mantica F, Nabeel-Shah S, Parada G, Brown KR, O'Hanlon D, Wei J, Yao Y, Zid AA, Comsa LC, Jen M, Wang J, Datti A, Gonatopoulos-Pournatzis T, Weatheritt RJ, Greenblatt JF, Wrana JL, Irimia M, Gingras AC, Moffat J, Blencowe BJ. *Mol Cell* 2022 Nov 1;72(3):510-524.e12. [PubMed](https://pubmed.ncbi.nlm.nih.gov/35914530/)

Version-specific DOI:
[![DOI](https://zenodo.org/badge/97053989.svg)](https://zenodo.org/badge/latestdoi/97053989)


References
----------
A multiplexed, next generation sequencing platform for high-throughput detection of SARS-CoV-2. Aynaud MM, Hernandez JJ, Barutcu S, Braunschweig U, Chan K, Pearson JD, Trcka D, Prosser SL, Kim J, Barrios-Rodiles M, Jen M, Song S, Shen J, Bruce C, Hazlett B, Poutanen S, Attisano L, Bremner R, Blencowe BJ, Mazzulli T, Han H, Pelletier L, Wrana JL. *Nat Commun* 2021 Mar 3;12(1):1405. [PubMed](https://pubmed.ncbi.nlm.nih.gov/33658502/)

Genome-wide CRISPR-Cas9 Interrogation of Splicing Networks Reveals a Mechanism for Recognition of Autism-Misregulated Neuronal Microexons. Gonatopoulos-Pournatzis T, Wu M, Braunschweig U, Roth J, Han H, Best AJ, Raj B, Aregger M, O'Hanlon D, Ellis JD, Calarco JA, Moffat J, Gingras AC, Blencowe BJ. *Mol Cell* 2018 Nov 1;72(3):510-524.e12. [PubMed](https://pubmed.ncbi.nlm.nih.gov/30388412/)

Multilayered Control of Alternative Splicing Regulatory Networks by Transcription Factors. Han H, Braunschweig U, Gonatopoulos-Pournatzis T, Weatheritt RJ, Hirsch CL, Ha KC, Radovani E, Nabeel-Shah S, Sterne-Weiler T, Wang J, O'Hanlon D, Pan Q, Ray D, Zheng H, Vizeacoumar F, Datti A, Magomedova L, Cummins CL, Hughes TR, Greenblatt JF, Wrana JL, Moffat J, Blencowe BJ. *Mol Cell* 2017 Feb 2;65(3):539-553.e7. [PubMed](https://www.ncbi.nlm.nih.gov/pubmed/28157508)


Feedback
--------
For feedback and questions, mail to u.braunschweig@utoronto.ca.



