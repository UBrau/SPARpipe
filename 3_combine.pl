#!/usr/bin/env perl                                                                                                                                     
### U. Braunschweig 2017-2020
###
### Changes: - Removed option --bin

use strict;
use Getopt::Long;
use Cwd qw(abs_path);

my $helpFlag = 0;
my $juncBase;
my $evTab;
my $treatTab;
my $outDir;
my $cores  = 1;

# Initialize
my $path = abs_path($0);
$0 =~ s/^.*\///;
$path =~ s/\/$0$//;


GetOptions("help"        => \$helpFlag,
           "juncBase=s"  => \$juncBase,
           "eventTab=s"  => \$evTab,
           "treatTab=s"  => \$treatTab,
           "outDir:s"    => \$outDir,
           "cores:i"     => \$cores,
    );

### Input checking
# Usage message
if ($helpFlag | !defined($evTab)| !defined($treatTab)| !defined($juncBase) | !scalar(@ARGV) | defined($ARGV[1])) {
  die "
*** Extract raw PSI, RPM and pseudocounts from aligned SPAR-seq data ***
Input:  Output directory containing folder /map with BAM files produced by 2_align_reads.pl
Output: Tables with read counts per junction, one per sample

Usage: $0 --juncBase JUNCBASE --eventTab FILE --treatTab FILE DIR
Options:
         --juncBase  Base name of junction BED files. This is a pair (fwd/rev) of simple BED files
                     derived from the Junctions library, where every junction has a line with
                     junction name as chromosome, 0 as start, junction length as end.
         --eventTab  Event table specifying which junctions to use etc. Must contain columns
                     Gene, Event, Label, Relative, JunctionsFw, and JunctionsRv.
                     This was generated using MatchPrimers.R.
         --treatTab  Treatment table specifying all samples in the project, including replicates. 
                     Must contain columns ID, Replicate, Batch, Barcode.
        [--cores]    Computing cores [default: $cores]

";
}

# Check input
my $outDir = $ARGV[0];
die "Folder /map not found in $outDir\n" unless (-e $outDir."/map");
die "Event table not found\n" unless (-e $evTab);
die "Treatment table not found\n" unless (-e $treatTab);
die "Fwd junction BED not found at ".$juncBase."_fwd.bed\n" unless (-e $juncBase."_fwd.bed");
die "Rev junction BED not found at ".$juncBase."_rev.bed\n" unless (-e $juncBase."_rev.bed");

$outDir =~ s/\/$//;


# Find BAM files
opendir(DIR, $outDir."/map") or die;
my @files = grep(/.*\.bam$/, readdir(DIR));
closedir(DIR);

my %BATCHES;
foreach my $file (@files) {
    $file =~ m/([^_]*)_?W.+/;
    $BATCHES{$1} = 1;
}

my %FWF;
my %RVF;
my %UNQ;

foreach my $file (@files) {
    if ($file =~ /.+Fwd.bam/) {
        $file =~ m/(.+W[0-9]+)_.+/;
        $FWF{$1} = $file;
	$UNQ{$1} = 1
    } elsif ($file =~ /.+Rev.bam/) {
        $file =~ m/(.+W[0-9]+)_.+/;
        $RVF{$1} = $file;
	$UNQ{$1} = 1
    }
}


# Check if both fwd and rev are there and get counts
mkdir $outDir."/counts" unless (-e $outDir."/counts");
foreach my $sample (sort keys %UNQ) {
    die "Fwd BAM file not found for $sample\n" unless defined ($FWF{$sample});
    die "Rev BAM file not found for $sample\n" unless defined ($RVF{$sample});
    if ($sample eq "") {
	print "Skipping $FWF{$sample} and $RVF{$sample}\n";
	next;
    }
    
    system "$path/bin/junction_counts.sh $outDir/map/$FWF{$sample} $juncBase" and
	die "[3] Error when getting counts for $sample\n";
}
print "[3] Done getting junction counts\n";


# Generate files for expression analysis
mkdir $outDir."/expression" unless (-e $outDir."/expression");
system "$path/R/combine_reads_expression.R -e $evTab -t $treatTab -d fwd -c $cores $outDir" and
    die "[3] Error when combining reads for expression\n";
print "[3] Done generating tables for expression\n";


# Get PSI, RPM, pseudo counts etc. per well
mkdir $outDir."/welldata" unless (-e $outDir."/welldata");
system "$path/R/compute_psi_rpm.R -c $cores $outDir -e $evTab" and
    die "[3] Error when calculating PSI etc.\n";
print "[3] Done computing PSI\n";


# Merge tables
mkdir $outDir."/batchdata" unless (-e $outDir."/batchdata");
foreach my $batch (sort keys %BATCHES) {
    system "$path/R/merge_tables.R -b $batch -c $cores $outDir" and
	die "[3] Error when merging tables\n";
    print "[3] Done merging tables for batch $batch\n";
}


# Merge batches and reconcile with expected samples and events
mkdir $outDir."/raw" unless (-e $outDir."/raw");
system "$path/R/combine_batches.R -e $evTab -t $treatTab $outDir" and
    die "[3] Error when creating output tables\n";
print "[3] Done creating output tables\n";
