#!/usr/bin/env perl
### U. Braunschweig, 2017-2020

use strict;
use Getopt::Long;
use Cwd qw(abs_path);
use File::Spec;

my $helpFlag = 0;
my $junc;
my $outDir = ".";
my $cores  = 1;
my $trim5  = 0;
my $trim3  = 0;
my $btopt = "--best -v 3 -k 1 -m 1";
my $bin;

# Initialize
my $path = abs_path($0);
$0 =~ s/^.*\///;
$path =~ s/\/$0$//;

GetOptions("help"        => \$helpFlag,
           "junc=s"      => \$junc,
	   "outDir:s"    => \$outDir, 
	   "bin:s"       => \$bin,
	   "cores:i"     => \$cores,
	   "trim5:i"     => \$trim5,
	   "trim3:i"     => \$trim3,
	   "btopt:s"     => \$btopt
    );

### Input checking 
# Usage message
if ($helpFlag | !defined($junc) | !scalar(@ARGV)) {
  die "                                                                                                                                                                     
*** Aligning of SPAR-seq read data with junction library ***
Input:  Demultiplexed reads as (gzipped) FASTQ files
        Run separately for fwd and rev reads with appropriate junction library
Output: Mapped read BAM files, unmapped reads FASTQ files, mapping stats

Usage: $0 --junc=FILE [options] INPUTFILES
Options:
         INPUTFILES  One or more demultiplexed FASTQ files, optionally gzipped
                     (space separated or specified using wildcards)
         --junc      Bowtie index (base name) for junction library
        [--outDir]   Output directory [default: ./]
        [--cores]    Computing cores [default: $cores]
        [--trim5]    Nt to trim off the 5' end of the read [default: $trim5]
        [--trim3]    Nt to trim off the 3' end of the read [default: $trim3]
        [--btopt]    Bowtie mapping options (use double quotes) [default: $btopt]

";
}

# Create output dir if missing
$outDir =~ s/\/$//;
mkdir $outDir        unless (-e $outDir);
mkdir $outDir."/map" unless (-e $outDir."/map");

# Create log file
my $logFile = $outDir."/map.log";
open (my $logFH, ">>", $logFile)
    or die "cannot open > $logFile: $!";

# Run bowtie on each file
foreach my $input (@ARGV) { 
    ## Check conformity of file names
    my $volume;
    my $dirs;
    my $sfile;
    ($volume,$dirs,$sfile) = File::Spec->splitpath($input);

    $sfile =~ /^[^_]+_W[0-9]+.*_(Fwd|Rev)\.(fq|fastq)(\.gz|)$/ ||
	die "Input file name $sfile is misshapen. Must match pattern ^[^_]+_W[0-9]+.*_(Fwd|Rev)\.(fq|fastq)(\.gz|)$\n";
    
    my $cmdline = "$path/$0 --junc=$junc --outDir=$outDir --cores=$cores --trim5=$trim5 --trim3=$trim3";

    $cmdline .= " $btopt $input";
    print $logFH "$cmdline";
    system "$path/bin/run_bowtie.sh $input $junc ".$outDir."/map $cores $trim5 $trim3 \"$btopt\" $path/bin" and
	die "[2] Error during bowtie mapping of $input\n";

    $sfile =~ s/\.(fastq|fq).*//;
    my $ssize = -s "$outDir/map/$sfile.stats.txt";
    if ($ssize > 0) {
	print $logFH "\tOK\n";
    } else {
	print $logFH "\tERROR\n";
	die "[2] Error during bowtie mapping: Illegal bowtie input?\n"
    }
}

close ($logFH)
    || warn "close failed: $!";

