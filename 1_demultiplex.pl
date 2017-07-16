#!/usr/bin/perl
### U. Braunschweig 02/2017 using some code from Nuno
### Changes: - This version accepts gzipped files for R1 and R4
###          - Added input checking and more convenient interface

use strict;
use warnings;
use Getopt::Long;


### Initialize
my $ulimit = `ulimit -n | head`;  # increase ## of open files; '| head' because 'ulimit' is
$ulimit > 2000 or die "Insufficient ulimit! Type 'ulimit -n 2048'\n";  

my $helpFlag = 0;
my $barcodefile;
my $prefix = "L1";
my $outDir = ".";

GetOptions("help"        => \$helpFlag,
           "barcodes=s"  => \$barcodefile,
	   "prefix:s"    => \$prefix,
           "outDir:s"    => \$outDir
    );
 
### Input checking
# Usage message
if ($helpFlag | !defined($barcodefile) | !defined($ARGV[3])) {
  die "
*** Demultiplexing of double-barcoded SPAR-seq screen ***
using bowtied mappings to barcode libraries (i.e., from separate
barcode reads.
Output: One FASTQ file for every permissible barcode combination (= well)
        and one (per direction) for unassigned reads.
        Also output a table with read counts per file.

Usage: $0 [options] --barcodes <allowed barcode table> <R1.fq[.gz]> <FwBC.sam> <RevBC.sam> <R4.fq[.gz]>
Options:
         --barcodes  Table with columns well, i5 BC, i5 BC name, i7 BC, i7 BC name
        [--prefix]   Batch prefix, e.g. 'Lane1'. Must not contain , or _ [default: $prefix]
        [--outDir]   Output directory [default: ./]

";
}

# Capture input
my $FQfileR1  = $ARGV[0];
my $SAMfileR2 = $ARGV[1];
my $SAMfileR3 = $ARGV[2];
my $FQfileR4  = $ARGV[3];

if ($prefix ne "") {
    $prefix = $prefix."_"; 
}

my %fqfilename0; # Hash of file name bases
my %FH1; # Hash of filehandles for forward sequences
my %FH2; # Hash of filehandles for reverse sequences
my %tally; # Read counts per file

# Create output dir if missing
$outDir =~ s/\/$//;
mkdir $outDir unless (-e $outDir);

# Check input files
die "File not found: $barcodefile\n" unless -e $barcodefile;
die "File not found: $FQfileR1\n"    unless -e $FQfileR1;
die "File not found: $SAMfileR2\n"   unless -e $SAMfileR2;
die "File not found: $SAMfileR3\n"   unless -e $SAMfileR3;
die "File not found: $FQfileR4\n"    unless -e $FQfileR4;

# Read allowable barcodes
open(BARCODES,$barcodefile);
<BARCODES>; #header
while(<BARCODES>){
    no strict 'refs';  # to enable symbolic refs to the file name hashes
    chomp($_);
    my($well,$i5bcseq,$i5bcname,$i7bcseq,$i7bcname) = split(/\t/,$_);
    $well = sprintf("%3d", $well);
    $well =~ tr/ /0/;
    my $filename = "W".$well."_".$i5bcname."_".$i7bcname; # File name base for each type
    my $combbarcode = $i5bcseq."_".$i7bcseq; # Combined barcode as unique identifier
    $fqfilename0{$combbarcode} = $filename; 
    my $filename1 = $filename."1";
    my $filename2 = $filename."2";
    $FH1{$combbarcode} = *$filename1;
    $FH2{$combbarcode} = *$filename2;
    open($FH1{$combbarcode}, ">", $outDir."/".$prefix.$filename."_Fwd.fq"); # Fastq files for forward sequences
    open($FH2{$combbarcode}, ">", $outDir."/".$prefix.$filename."_Rev.fq"); # Fastq files for reverse sequences
}
do { # Open files for reads with undefined barcode combinations
    no strict 'refs';  # to enable symbolic refs to the file name hashes
    my $filename = "Unassigned";
    $fqfilename0{"NA"} = $filename;
    my $filename1 = $filename."1";
    my $filename2 = $filename."2";
    $FH1{"NA"} = *$filename1;
    $FH2{"NA"} = *$filename2;
    open($FH1{"NA"},">", $outDir."/".$prefix.$filename."_Fwd.fq"); # Fastq files for forward sequences
    open($FH2{"NA"},">", $outDir."/".$prefix.$filename."_Rev.fq"); # Fastq files for reverse sequences
};
close BARCODES;

# Read data and distribute to barcoded files
open (R1, "gzip -dc $FQfileR1 |") or die "cannot open $FQfileR1\n"; 
open (R2, "<", $SAMfileR2) or die "cannot open $SAMfileR2\n"; 
open (R3, "<", $SAMfileR3) or die "cannot open $SAMfileR3\n"; 
open (R4, "gzip -dc $FQfileR4 |") or die "cannot open $FQfileR4\n"; 
my ($lineR2, $lineR3);
my ($nameR1, $nameR4);
my ($seqR1, $seqR4);
my @readR2;
my @readR3;
my $combbarcode;
my $FH;
my ($R2bc, $R3bc);
my $linecount = 1;  # counter for lines in FASTQ

while ($lineR2 = <R2>) {
    no strict 'refs';  # to enable symbolic refs to the file name hashes
    $lineR3 = <R3>;
    chomp $lineR2;
    chomp $lineR3;
    @readR2 = split(/\t/, $lineR2);
    @readR3 = split(/\t/, $lineR3);

    # read lines from barcode files and combine
    if ($readR2[2] =~ /[^_]+_([ACGT]+)/) {
	$R2bc = $1;
    }
    else {
	$R2bc = "NA";
    }

    if ($readR3[2] =~ /[^_]+_([ACGT]+)/) {
	$R3bc = $1;
    }
    else {
	$R3bc = "NA";
    }
    
    $combbarcode = $R2bc."_".$R3bc;
        if (!$fqfilename0{$combbarcode}){
	$combbarcode = "NA";
    }

    # Write 4 consecutive lines of FASTQ file R1 into appropriate file
    $tally{$combbarcode}++;
    $FH = $FH1{$combbarcode};
    while ($linecount <= 4) {
	$_ = <R1>;
	print $FH $_;
	$linecount++;
    }
    $linecount = 1;
	
    # Same for FASTQ file R4 
    $FH = $FH2{$combbarcode};
    while ($linecount <= 4) {
	$_ = <R4>;
	print $FH $_;
	$linecount++;
    }
    $linecount = 1
	
}  # end loop over barcode file lines

foreach my $combbc (sort keys %fqfilename0){
    my $FH = $FH1{$combbc};
    close $FH;
    $FH = $FH2{$combbc};
    close $FH;
}
close R1;
close R2;
close R3;
close R4;

## Save read count table
open TABLE, ">", $outDir."/"."$prefix"."ReadsPerWell.tab";
print TABLE "well\treads\n";
foreach my $ind (sort { $tally{$b} <=> $tally{$a} } keys %tally) {  # sort by # of reads
    print TABLE "$prefix$fqfilename0{$ind}\t$tally{$ind}\n";
}

close TABLE;
