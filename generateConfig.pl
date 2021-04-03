#!/usr/bin/env perl

# GOAL: Generate config file for panelCNVsim

use strict;
use warnings;
use Sort::Key::Natural qw(natsort);
use File::Basename;
use Getopt::Long;

my $bamDir;
my $ROI_file;
my $mode;
my $cnv_type;
my $autosomes = 0;

# Min max sizes for partial cnvs
my $minSize = "";
my $maxSize = "";
my $minExons = "";
my $maxExons = "";

Help () if (@ARGV < 4 or !GetOptions(
	'i|indir=s' =>\$bamDir,
	'b|bed=s' =>\$ROI_file,
	'm|mode=s' =>\$mode,
	'c|class=s' =>\$cnv_type,
	'minsize=i' =>\$minSize,
	'maxsize=i' =>\$maxSize,
	'minexons=i' =>\$minExons,
	'maxexons=i' =>\$maxExons,
	'autosomes'=>\$autosomes
	)
);

if ($mode ne 'partial' && $mode ne 'single' && $mode ne 'multiple' && $mode ne 'whole') {
	print " ERROR: $mode mode is not valid. Please, select one of the following modes: partial, single, multiple or whole\n";
	exit;
}

if ($cnv_type ne 'del' && $cnv_type ne 'dup'){
	print " ERROR: $cnv_type cnv type is not valid. Please, choose between: del or dup\n";
	exit;
} 

if (!$ROI_file || !-e $ROI_file) {
	print " ERROR: No input ROI bed file\n";
	exit;
}

# Default sizes for partial CNVs
if (!$minSize){
	if ($mode eq 'partial'){
		$minSize = 20;
	} 
}
if (!$maxSize){
	if ($mode eq 'partial'){
		$maxSize = 200;
	} 
}

# Default sizes for multiple CNVs
if (!$minExons){
	if ($mode eq 'multiple'){
		$minExons = 2;
	} 
}
if (!$maxExons){
	if ($mode eq 'multiple'){
		$maxExons = 8;
	} 
}

if ($minExons && $maxExons){
	if ($minExons > $maxExons){
		print " ERROR: --maxexons must be greater than --minexons\n";   
		exit;
	} 	
} 

# Checking dependencies
my $samtools = `which samtools`;
chomp $samtools;

if (!$samtools) {
	print " ERROR: samtools was not found on path\n"; 
	exit;
}

our @bams = glob ("$bamDir/*.bam");
if (!@bams) {
	print " ERROR: no BAM files were found on $bamDir directory!\n"; 
	exit;
}

my @ARRAY = ();
my %Exon  = ();
our %ExonCount = ();
our %GeneHash  = ();
our %seen_gene = ();
our %seen_whole= ();

open (IN, "<", $ROI_file) || die " ERROR: Unable to open $ROI_file\n";
while (my $line=<IN>) {
	chomp $line;

	my @tmp = split (/\t/, $line);
	my $chr = $tmp[0]; 
	if ($autosomes){
		next if $chr =~/X/;
		next if $chr =~/Y/;
	} 

	$Exon{"$tmp[0]\t$tmp[1]\t$tmp[2]"} = $tmp[3];
	my $gene = $tmp[3]; 

	if ($gene =~/.{1,}_.{1,}_.{1,}_.{1,};.{1,}$/) {
		my @tmpGene = split(/;/, $gene);
		$gene = $tmpGene[1]; 
	} 

	$ExonCount{$gene}++;
	$GeneHash{$gene}{$ExonCount{$gene}} = "$tmp[0]\t$tmp[1]\t$tmp[2]\n";
	push @ARRAY, $line;
}
close IN;


my %seen = ();
my $count = 1;

foreach my $bam (natsort@bams) {

	my ($chr, $start, $end, $gene);
	if ($mode eq 'partial') {
		($chr, $start, $end, $gene) = getPartialExon(@ARRAY);
	}
	if ($mode eq 'single') {
		($chr, $start, $end, $gene) = getSingleExon(@ARRAY);
	}
	if ($mode eq 'multiple') {
		($chr, $start, $end, $gene) = getMultipleExon(@ARRAY);
	}

	my $ploidy = $cnv_type eq 'DEL' ? '1' : '3';
	print "$bam\t$chr\t$start\t$end\t$gene\t$cnv_type\t$ploidy\t$mode\n";
	$count++;
}
#####################
sub getPartialExon {
	my @ARRAY = @_;

	# Selecting a random exon from the BED list
	my $index = int (rand(scalar @ARRAY));

	my ($chr, $start, $end, $info) = split (/\t/, $ARRAY[$index]);
	my $random = 1;
	my $new_start;
	my $new_end;
	my $x = $minSize + int(rand($maxSize - $minSize));
	  
	if( $random == 1) {
		$new_start = $start+10;
		$new_end   =$new_start+$x;  
	}
	return ($chr, $new_start, $new_end, $info);
}

#####################
sub getSingleExon {
	my @ARRAY = @_;

	my($chr, $start, $end, $info);
	while (1) {

		my $index = int (rand(scalar @ARRAY));
		($chr, $start, $end, $info) = split (/\t/, $ARRAY[$index]);

		my ($chrPrev, $startPrev, $endPrev, $infoPrev) = split (/\t/, $ARRAY[$index-1]);
		my ($chrNext, $startNext, $endNext, $infoNext) = split (/\t/, $ARRAY[$index+1]);

		if ($chr eq $chrPrev && $chr eq $chrNext){
			if ($start > $endPrev+500 && $startNext > $end+500){
				$seen{$ARRAY[$index]}++;
				last;
			}  
		}
		if ($chr eq $chrPrev && $chr ne $chrNext){
			if ($start > $endPrev+500){
				$seen{$ARRAY[$index]}++;
				last;				
			} 
		} 
		if ($chr eq $chrNext && $chr ne $chrPrev){
			if ($startNext > $end+500){
				$seen{$ARRAY[$index]}++;
				last;				
			} 
		} 
	}	

	$start = $start-500;
	$end += 500;
	return ($chr, $start, $end, $info);
}

#####################
sub getMultipleExon {

	my @ARRAY = @_;

	# Conditions: gene must have >= chosen_N
	my $chosen_N = int ($minExons+rand($maxExons-$minExons));

	my $selectedGene;
	my $counter = 0;

	while (1) {
		$counter++;
		if ($counter == scalar @ARRAY) {
			$chosen_N-- if $chosen_N >= 2;
			$counter = 0;
		}

		my $index = int (rand (scalar @ARRAY));
		my ($chr, $start, $end, $info) = split (/\t/, $ARRAY[$index]);
		my $gene = $info;
		if ($gene =~/.{1,}_.{1,}_.{1,}_.{1,};.{1,}$/) {
			my @tmpGene = split(/;/, $gene);
			$gene = $tmpGene[1]; 
		} 

		if ($ExonCount{$gene} >= $chosen_N) {
			$selectedGene = $gene;

			$seen_gene{$gene}++;
			my ($chrPrev, $startPrev, $endPrev, $infoPrev) = split (/\t/, $ARRAY[$index-1]);
			my ($chrNext, $startNext, $endNext, $infoNext) = split (/\t/, $ARRAY[$index+1]);

			if ($seen_gene{$gene} > 1) {
				next;
			}
			else {
				last;
			}

			if ($chr eq $chrPrev && $chr eq $chrNext){
				if ($start > $endPrev+500 && $startNext > $end+500){
					$seen{$ARRAY[$index]}++;
					last;
				}  
			}
		}	
	}

	my @tmp = (); 
	foreach my $exon (natsort keys %{ $GeneHash{$selectedGene} } ){
		push @tmp, $GeneHash{$selectedGene}{$exon}; 
	}

	my @sortedTmp = natsort @tmp[0..$chosen_N-1];

	my $chr = ( split /\t/, $sortedTmp[0] )[0];
	my $initial = ( split /\t/, $sortedTmp[0] )[1];
	my $final   = ( split /\t/, $sortedTmp[-1] )[2];   

	$initial = $initial-500;
	$final   = $final+500;

	return ($chr, $initial, $final, "$selectedGene\_1_to_$chosen_N");
}

##############
sub Help {
print "\n Usage: perl $0 <options>

 Options:
 -i,--indir  STRING   Input BAM directory
 -b,--bed    STRING   ROI bed
 -m,--mode   STRING   Mode (partial, single, multiple)
 -c,--class  STRING   CNV class (del, dup)
 --minsize   INTEGER  Minimum size for partial CNV (default=20)
 --maxsize   INTEGER  Maximum size for partial CNV (default=200)
 --minexons  INTEGER  Minimum number of exons for multiple CNVs (default=2)
 --maxexons  INTEGER  Maximum number of exons for multiple CNVs (default=2)
 --autosomes          Simulate only over autosomes\n";exit;
}