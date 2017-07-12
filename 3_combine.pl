#!/usr/bin/env perl                                                                                                                                     
### U. Braunschweig, 04/2017
###
### Changes: New option --pattern

use strict;
use Getopt::Long;
use Cwd qw(abs_path);

my $helpFlag = 0;
my $juncBase;
my $pattern = "";
my $outDir;
my $cores  = 2;
my $bin;

# Initialize
my $path = abs_path($0);
$0 =~ s/^.*\///;
$path =~ s/\/$0$//;


GetOptions("help"        => \$helpFlag,
           "juncBase=s"  => \$juncBase,
	   "pattern:s"   => \$pattern,
           "outDir:s"    => \$outDir,
           "bin:s"       => \$bin,
           "cores:i"     => \$cores,
    );

### Input checking
# Usage message
if ($helpFlag | !defined($juncBase) | !scalar(@ARGV) | defined($ARGV[1])) {
  die "
*** Extract raw PSI, RPM and pseudocounts from aligned SPAR-seq data ***
Input:  Output directory containing folder /map with BAM files
        (*sorted.bam) produced by 2_align_reads.pl
Output: Tables with read counts per junction, one per sample

Usage: $0 --juncBase JUNCBASE DIR
Options:
         --juncBase  Base name of junction BED files
        [--pattern]  Pattern to look for in BAM files. Use this if there are e.g. multiple lanes
                     and run each lane separately.
        [--bin]      Directory containing scripts [default: $path]
        [--cores]    Computing cores [default: $cores]

";
}

# Check input
$bin = $path unless (defined $bin);

my $outDir = $ARGV[0];
die "Folder /map not found in $outDir\n" unless (-e $outDir."/map");
die "Fwd junction BED not found at ".$juncBase."_fwd.bed\n" unless (-e $juncBase."_fwd.bed");
die "Rev junction BED not found at ".$juncBase."_rev.bed\n" unless (-e $juncBase."_rev.bed");

# Create output dir if missing
$outDir =~ s/\/$//;
mkdir $outDir."/counts" unless (-e $outDir."/counts");
mkdir $outDir."/welldata" unless (-e $outDir."/welldata");
mkdir $outDir."/batchdata" unless (-e $outDir."/batchdata");

# Find BAM files
opendir(DIR, $outDir."/map") or die;
#my @files = grep(/\.sorted.bam$/, readdir(DIR));
my $re = qr/$pattern/;
my @files = grep(/${re}.*sorted\.bam/, readdir(DIR));
closedir(DIR);
    
my %FWF;
my %RVF;
my %UNQ;

foreach my $file (@files) {
    if ($file =~ /.+Fwd.sorted.bam/) {
        $file =~ m/.+(W[0-9]+)_FwBC.+/;
        $FWF{$1} = $file;
	$UNQ{$1} = 1
    } elsif ($file =~ /.+Rev.sorted.bam/) {
        $file =~ m/.+(W[0-9]+)_FwBC.+/;
        $RVF{$1} = $file;
	$UNQ{$1} = 1
    }
}

# Check if both fwd and rev are there and get counts
foreach my $well (sort keys %UNQ) {
    die "Fwd BAM file not found for $well\n" unless defined ($FWF{$well});
    die "Rev BAM file not found for $well\n" unless defined ($RVF{$well});
    if ($well eq "") {
	print "Skipping $FWF{$well} and $RVF{$well}\n";
	next;
    }
    
    system "$bin/bin/junction_counts.sh $outDir/map/$FWF{$well} $juncBase" and
	die "[3] Error when getting counts for $well\n";
}
print "[3] Done getting junction counts\n";

# Get PSI, RPM, pseudo counts etc. per well
system "$bin/R/compute_psi_rpm.R -c $cores $outDir" and
    die "[3] Error when calculating PSI etc.\n";
print "[3] Done computing PSI\n";

# Merge tables
system "$bin/R/merge_tables.R -c $cores $outDir" and
    die "[3] Error when merging tables\n";
print "[3] Done merging tables\n";
