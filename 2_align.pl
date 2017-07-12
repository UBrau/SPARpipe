#!/usr/bin/env perl
### U. Braunschweig, 04/2017

use strict;
use Getopt::Long;
use Cwd qw(abs_path);

my $helpFlag = 0;
my $junc;
my $outDir = ".";
my $cores  = 2;
my $trim5  = 5;
my $trim3  = 47;
my $btopt = "--best -v 3 -k 1";
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
Output: Mapped read BAM files, unmapped reads FASTQ files, mapping stats

Usage: $0 --junc=FILE [options] FILE
Options:
         --junc      Bowtie index (base name) for junction library
        [--outDir]   Output directory [default: ./]
        [--bin]      Directory containing scripts [default: $path]
        [--cores]    Computing cores [default: $cores]
        [--trim5]    Nt to trim off the 5' end of the read [default: $trim5]
        [--trim3]    Nt to trim off the 6' end of the read [default: $trim3]
        [--btopt]    Bowtie mapping options (use double quotes) [default: $btopt]

";
}

# Complete input
$bin = $path unless (defined $bin);

# Create output dir if missing
$outDir =~ s/\/$//;
mkdir $outDir        unless (-e $outDir);
mkdir $outDir."/map" unless (-e $outDir."/map");

# Run bowtie on each file
foreach my $input (@ARGV) {
    system "$bin/bin/run_bowtie.sh $input $junc ".$outDir."/map $cores $trim5 $trim3 \"$btopt\" $bin/bin" and
	die "[2] Error during bowtie mapping of $input\n";

}



