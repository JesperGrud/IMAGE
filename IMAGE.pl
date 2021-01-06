#!/usr/bin/perl -w

# Import necessary packages
use POSIX;
use Scalar::Util qw(looks_like_number);
use Cwd 'abs_path';

# Define installation directory
my $INSTALLDIR = abs_path($0);
$INSTALLDIR = substr $INSTALLDIR, 0, -9;

# Initialize variables
my $start_run = time();
my $resizesize=200;
my $CPUCount = 1;
my @regionfields=();
my $regionlinecounter=0;
my $regionarraysize=0;
my $regionarraysizetotal=0;
my $inputfile="";
my $processors=2;
my @RNADesign=();
my @EnhancerDesign=();
my @exonfields=();
my $exonlinecounter=0;
my $exonarraysize=0;
my $exonfile="";
my $string;
my $pairing=0;
my $fasta="";
my $threshold=0;
my $RESULTNAME="Result";
my $targets="no";
my $TMP = "/tmp";
my $UTIL = "/utils";
my $TMPDIR = $INSTALLDIR . $TMP;
my $UTILDIR = $INSTALLDIR . $UTIL;

my @chars = ("A".."Z", "a".."z");

# Help functions
sub printHelp {
	print STDERR "\n\tUsage: IMAGE.pl -region [BED file] -expression [ GENE file ] -fasta [ FASTA file ] -RNADesign [ Grouping levels ] -EnhancerDesign [ Grouping levels ] -n [ Prefix of output files ] -p [ Number of processors ] (-paired) (-targets) \n";
	print STDERR "\n\t\tObligatory parameters\n";
	print STDERR "\t\t\t-region\t\t:\t Path to BED file with normalized enhancers counts. \n";
	print STDERR "\t\t\t-expression\t:\t Path to GENE file with raw gene expression: At least 2 samples per condition \n";
	print STDERR "\t\t\t-fasta\t\t:\t Path to FASTA file of the genome. \n";
	print STDERR "\t\t\t-RNADesign\t:\t The grouping factor provides information about the number of replicates.\n";
 	print STDERR "\t\t\t\t\t\t\t Example: 4 conditions, 3 sets of samples. Run: -RNADesign 1 1 1 2 2 2 3 3 3 4 4 4\n";
	print STDERR "\t\t\t-EnhancerDesign\t:\t The grouping factor provides information about the number of replicates.\n";
	print STDERR "\t\t\t\t\t\t\t Example: 4 conditions, 2 sets of samples. Run: -EnhancerDesign 1 1 2 2 3 3 4 4";
	print STDERR "\n";
	print STDERR "\n\t\tOptional parameters\n";
	print STDERR "\t\t\t-p\t\t:\t The number of processors to use. Needs to be larger than 1. Default = 2\n";
	print STDERR "\t\t\t-threshold\t:\t Threshold for removing too lowly expressed transcription factors.  Default = 0\n";
	print STDERR "\t\t\t-n\t\t:\t The prefix of result files. \n";
	print STDERR "\t\t\t-paired\t\t:\t Indicate if the gene expression samples are paired \n";
	print STDERR "\t\t\t-targets\t:\t Indicate if the predicted target genes of ALL motifs with TFs expressed above threshold should be determined.";
	print STDERR "\n\t\t\t\t\t\t Otherwise only target genes of putative regulators is determined. (Can take a long time and use a lot of memory!) \n";
	print STDERR "\n\t For examples and explanation of BED file and GENE file. See the examples folder.";
	print STDERR "\n\n";
	exit;
}

print STDOUT "\n#### Welcome to IMAGE ####";

# Check that the necessary programs are installed
my $HOMER_path = `which annotatePeaks.pl 2>/dev/null`;
print STDERR "\n\nHOMER is not executable, please see: homer.salk.edu/homer/ \n\n" unless ( $HOMER_path );
exit unless ( $HOMER_path );

my $R_path = `which R 2>/dev/null`;
print STDERR "\n\nR is not executable, please see: www.r-project.org/ \n\n" unless ( $R_path );
exit unless ( $R_path );

my $Package = `Rscript $UTILDIR/Check_Packages.R`;
if ($Package eq 'FALSE') {
print STDERR "\n\nThe required packages are not installed under R, please see: http://http://bioconductor.org/\n";
print STDERR "The required packages are:\n\tglmnet\n\tdata.table\n\t\n\tcluster\n\tdoParallel\n\tforeach\n\tMatrix\n\tGenomicRanges\n\tmethods\n\n";
exit;
}

## Check that users has supplied parameters otherwise exit
if (@ARGV == 0 && -t STDIN && -t STDERR) { 
	print STDOUT "\n\nPlease supply input parameters. Exitting\n\n";
    printHelp()
}

## Read user inputs
for (my $i=0;$i<@ARGV;$i++) {
		if ($ARGV[$i] eq '-region') {
			$inputfile = $ARGV[($i + 1)];
        }
		if ($ARGV[$i] eq '-expression') {
			$exonfile = $ARGV[($i + 1)];
        }
		if ($ARGV[$i] eq '-RNADesign') {
			$end = @ARGV;
			$m = $i + 1;
			while ($m < @ARGV) {
				if ($ARGV[$m] =~ /^\-/) { $end = $m; last; }
				$m = $m + 1;
			}
			$start = $i+1;
			for (my $m=$start;$m<$end;$m++) {
				push(@RNADesign, $ARGV[$m]);
			}
		}	
		if ($ARGV[$i] eq '-paired') {
			$pairing = 1;
        }

        if ($ARGV[$i] eq '-EnhancerDesign') {
			$end = @ARGV;
			$m = $i + 1;
            while ($m < @ARGV) {
                if ($ARGV[$m] =~ /^\-/) { $end = $m; last; }
                $m = $m + 1;
            }
            $start = $i+1;
            for (my $m=$start;$m<$end;$m++) {
                push(@EnhancerDesign, $ARGV[$m]);
            }
        }
		if ($ARGV[$i] eq '-threshold') {
			$threshold = $ARGV[($i + 1)];
        }
		if ($ARGV[$i] eq '-n') {
			$RESULTNAME = $ARGV[($i + 1)];
        }
		if ($ARGV[$i] eq '-p') {
			$processors = $ARGV[($i + 1)];
        }
		if ($ARGV[$i] eq '-fasta') {
			$fasta = $ARGV[($i + 1)];
        }
		if ($ARGV[$i] eq '-targets') {
			$targets = "yes";
        }
}

## Check the user supplied parameters for sanity
# Check that the selected number of processors is not larger than the number of the computer
if(open CPU, "/proc/cpuinfo" ) {
	$CPUCount = scalar (map /^processor/, <CPU>);  
	close CPU;
} else {
	print STDOUT "\n\nCould not open /proc/cpuinfo to get number of processors. Is this really linux? Exitting\n\n";
	exit;
}
if ($CPUCount < $processors) { print STDOUT "\n\nYou have selected more processors than your computer have. Exitting\n\n"; exit; }

# Check that genome is selected and can be opened
if ($fasta eq "" || ! -e $fasta) { print STDOUT "\n\nPlease select a valid fasta file. Exitting\n\n"; printHelp() }

# Check if region file exists and read into array. Use array to check the format and get number of conditions½
if(open FILE, $inputfile) {
	while(<FILE>) { 
		chomp; 
		@regionfields = split('\t', $_); 
		$regionlinecounter++;
	}
	close FILE;
} else {
	print STDOUT "\n\nRegion file could not be opened. Exitting\n\n";
	printHelp()
}
$regionarraysize = scalar @regionfields;
$regionarraysize = $regionarraysize - 4;
if (looks_like_number($regionfields[1]) == 0 || looks_like_number($regionfields[2]) == 0) { 
	print STDOUT "\n\nInput file does not look formatted correctly. Exitting\n\n"; 
	printHelp() 
}
for (my $i=4;$i<@regionfields;$i++) { 
        if (looks_like_number($regionfields[$i]) == 0) { 
			print STDOUT "\n\nInput file does not look formatted correctly. Exitting\n\n"; 
			printHelp()
		}
}

# Check if exon file exists and read into array. Use array to check the format and get number of conditions
if(open FILE, $exonfile) {
	while(<FILE>) { 
		chomp; 
		@exonfields = split('\t', $_); 
		$exonlinecounter++;
	}
	close FILE;
} else {
	print STDOUT "\n\nExon file could not be opened. Exitting\n\n";
	printHelp()
}
$exonarraysize = scalar @exonfields;
$exonarraysize = $exonarraysize - 8;
for (my $i=8;$i<@exonfields;$i++) { 
        if (looks_like_number($exonfields[$i]) == 0) { 
			print STDOUT "\n\nExon file does not look formatted correctly. Exitting\n\n"; 
			printHelp()
		}
}

# Check that RNA design array is defined. If it is defined get number of conditions and number of files
if (@RNADesign) { 
	my @unique_design_RNA = do { my %seen; grep { !$seen{$_}++ } @RNADesign };
	$RNADesignarrayunique = scalar @unique_design_RNA;
	$RNADesignarray = scalar @RNADesign;
} else {
	print STDOUT "\n\nRNA design not given. Please see help. Exitting\n\n"; 
	printHelp()
}

# Check that enhancer design array is defined. If it is defined get number of conditions and number of files
if (@EnhancerDesign) {
	my @unique_design_Enhancers = do { my %seen; grep { !$seen{$_}++ } @EnhancerDesign };
	$EnhancerDesignarrayunique = scalar @unique_design_Enhancers;
	$EnhancerDesignarray = scalar @EnhancerDesign;
} else {
	print STDOUT "\n\nEnhancer design not given. Please see help. Exitting\n\n";
printHelp()
}

#print STDOUT $RNADesignarray;
#print STDOUT $exonarraysize;

# Match to number of conditions in design to genes and regions
if ($RNADesignarray != $exonarraysize || $RNADesignarrayunique != $EnhancerDesignarrayunique) {
print STDOUT "\n\nThere is an unequal number of conditions in the RNA design array compared to expression or between RNA and enhancers. Exitting\n\n";
printHelp() 
}

## If all sanity checks are passed. Output summary of input data.
$exonlinecounter = $exonlinecounter - 1;
print STDOUT "\n\nYour region file contains $regionlinecounter regions and $EnhancerDesignarrayunique conditions across $EnhancerDesignarray files";
print STDOUT "\nYour exon file contains $exonlinecounter genes and $RNADesignarrayunique conditions across $exonarraysize files";
print STDOUT "\nYour fasta file is $fasta";
print STDOUT "\nYou are using $processors processors.";
if ($threshold != 0) { print STDOUT "\nYou are using a non-standard threshold for filtering transcription factors. Your threshold is $threshold";}
if ($pairing == 1) { print STDOUT "\nYour samples will be analyzed using paired statistics"; }

## Generate random strings
$string .= $chars[rand @chars] for 1..8;

## Setup the region file (using UNIX) and scan for motifs (using HOMER)
print STDOUT "\n\nStarting the analysis";
print STDOUT "\n\tSetting up for parallizing motif search";
my $pid = fork;
if (not $pid) {
	`sh $UTILDIR/Parallel.sh $processors $INSTALLDIR $string`;
	exit;
}

wait();

print STDOUT "\n\tPreparing input file for motif searching";
system("cut -f 1,2,3,4 $inputfile > $TMPDIR/$string.bed");
system("bed2pos.pl $TMPDIR/$string.bed > $TMPDIR/$string.pos");
system("resizePosFile.pl $TMPDIR/$string.pos $resizesize > $TMPDIR/$string.tmp");
system("mv $TMPDIR/$string.tmp $TMPDIR/$string.pos");
system("$UTILDIR/extractSequence $TMPDIR/$string.pos $fasta > $TMPDIR/$string.tab 2> /dev/null");

my $MOTIF = "_motifs";
my $MOTIFSTRING = $string . $MOTIF;
my $MOTIFDIR = $TMPDIR . "/" . $MOTIFSTRING;

opendir DIR, $MOTIFDIR;
my @files = grep { $_ ne '.' && $_ ne '..' } readdir DIR;
closedir DIR;

my $processorCount=$processors-1;

print STDOUT "\n\tScanning for motifs";
for (my $i=0; $i <= $processorCount; $i++) { 
	my $pid = fork;
	if (not $pid) {
		`homer2 find -s $TMPDIR/$string.tab -m $MOTIFDIR/$files[$i]  2> /dev/null > $TMPDIR/$string.sub.$i.count`;
		exit;
	}
}

for (my $i=0; $i <= $processorCount; $i++) {
	wait();
} 

`cat $TMPDIR/$string.sub.* > $TMPDIR/$string.count`;
`rm $TMPDIR/$string.sub.*`;
`rm -rf $MOTIFDIR`;

# Call Rscript to do all the stuffs in 3 stages
print STDOUT "\n\tRunning analysis in R";
my $stage=1;
open (R, "Rscript $UTILDIR/Regression.R $INSTALLDIR $string $inputfile $exonfile $pairing $RESULTNAME $processors $threshold @EnhancerDesign $targets $stage @RNADesign|");
while ( <R> ) {
        print STDOUT;
}

$stage=2;
open (R, "Rscript $UTILDIR/Regression.R $INSTALLDIR $string $inputfile $exonfile $pairing $RESULTNAME $processors $threshold @EnhancerDesign $targets $stage @RNADesign|");
while ( <R> ) {
        print STDOUT;
}

$stage=3;
open (R, "Rscript $UTILDIR/Regression.R $INSTALLDIR $string $inputfile $exonfile $pairing $RESULTNAME $processors $threshold @EnhancerDesign $targets $stage @RNADesign|");
while ( <R> ) {
        print STDOUT;
}

# Remove the tmp files
`rm $TMPDIR/$string*`;

# Print the total runtime
my $end_run = time();
my $run_time = $end_run - $start_run;
print STDOUT "IMAGE completed in $run_time seconds\n";
