#!/usr/bin/env perl

use strict;
use warnings;
use Cwd;
use File::Basename;
use Sort::Key::Natural qw(natsort);
use Getopt::Long;
use Config;

my $version = "1.0";
my $cnv_list;
my $ROI_file;
my $genome;
my $outDir;
my $keepTmp;
my $breakpoint;

Help () if (@ARGV < 4 or !GetOptions(
	'b|bed=s' =>\$ROI_file,
	'g|genome=s' =>\$genome,
	'c|config=s' =>\$cnv_list,
	'o|outdir=s' =>\$outDir,
	'keep_tmp' =>\$keepTmp,
	'breakpoint' =>\$breakpoint, 
	)
);	
my $dirname = dirname(__FILE__);

# Checking dependencies
my $bwa;
my $samtools;
my $devNull = ">/dev/null 2>&1";

# Setting binaries for either mac (darwin) or Linux systems
if ($Config{osname} =~/darwin/) {
	$bwa = ( -x "$dirname/bin/darwin/bwa" )
	? "$dirname/bin/darwin/bwa"
	:  die " ERROR: Unable to execute bwa";

	$samtools = ( -x "$dirname/bin/darwin/samtools" )
	? "$dirname/bin/darwin/samtools"
	:  die " ERROR: Unable to execute samtools";
}
else {
	$bwa = ( -x "$dirname/bin/linux/bwa" )
	? "$dirname/bin/linux/bwa"
	:  die " ERROR: Unable to execute bwa";

	$samtools = ( -x "$dirname/bin/linux/samtools" )
	? "$dirname/bin/linux/samtools"
	:  die " ERROR: Unable to execute samtools";
}

my  $awk  = 
  ( -x "/usr/bin/awk" )
  ? "/usr/bin/awk"
  : die "ERROR: awk not found\n";
my  $split  = 
  ( -x "/usr/bin/split" )
  ? "/usr/bin/split"
  : die "ERROR: split not found\n";
my  $head  = 
  ( -x "/usr/bin/head" )
  ? "/usr/bin/head"
  : die "ERROR: head not found\n";
my  $grep  = 
  ( -x "/bin/grep" )
  ? "/bin/grep"
  : die "ERROR: grep not found\n";
my  $uniq  = 
  ( -x "/usr/bin/uniq" )
  ? "/usr/bin/uniq"
  : die "ERROR: uniq not found\n";
my  $paste  = 
  ( -x "/usr/bin/paste" )
  ? "/usr/bin/paste"
  : die "ERROR: paste not found\n";
my  $sort  = 
  ( -x "/usr/bin/sort" )
  ? "/usr/bin/sort"
  : die "ERROR: sort not found\n";
my $cat  = 
  ( -x "/bin/cat" )
  ? "/bin/cat"
  : die "ERROR: cat not found\n";
my $cut  = 
  ( -x "/usr/bin/cut" )
  ? "/usr/bin/cut"
  : die "ERROR: cut not found\n";
my $tr  = 
  ( -x "/usr/bin/tr" )
  ? "/usr/bin/tr"
  : die "ERROR: tr not found\n";
my $tail  = 
  ( -x "/usr/bin/tail" )
  ? "/usr/bin/tail"
  : die "ERROR: tail not found\n";

if (!-e $outDir){
	mkdir $outDir;
} 

if (!$ROI_file) {
	print " ERROR: No input ROI file was entered\n\n"; 
	Help();
}
if (!$genome) {
	print " ERROR: No input genome file was entered\n\n"; 
	Help(); 
}
if (!$cnv_list) {
	print " ERROR: No input cnv list was entered\n\n"; 
	Help(); 
}
my $dir = cwd(); # Current working dir

# Unmappable sequence to the human genome 
my $unmap = "tgccaacgatgctgttatcgctaattcagttgctcaggcacgcttttcaggcttgttgattgccaacgatgctgttatcgctaattcagttgctcaggcacgcttttcaggcttgttgat";

my %HoS = getSampleInfo($cnv_list, $outDir);
my $N = 0;
foreach my $bam ( natsort keys %HoS ) {
	$N++;
	next if $HoS{$bam}{CNV} eq "none";
	my @cnvs = @{$HoS{$bam}{CNV}};

	my @arrOfHashes = (); 
	my %ReadName = ();

	foreach my $cnv (@cnvs) {

		print " INFO:$bam Processing $cnv\n";
		my @tmp = split (/\t/, $cnv);
		my $chrCnv  = $tmp[0];
		my $startCnv= $tmp[1];
		my $endCnv  = $tmp[2];
		my $cnvType = $tmp[4];
		my $ploidy  = $tmp[5];
		my $mode    = $tmp[6];

		my $startOffset = $startCnv;
		my $endOffset   = $endCnv;

		# Extract reads overlapping the CNV coordinates	
		my $roiReads = `$samtools view $bam $chrCnv:$startOffset-$endOffset`;
		chomp $roiReads;

		# Classify reads depending on its position relative to breakpoints	
		my @tmpReads = split (/\n/, $roiReads);
		foreach my $read (@tmpReads) { 
			my @sam = split (/\t/, $read);  
			my $readname = $sam[0] ;
			my $flag= $sam[1]; 
			my $chr = $sam[2]; 
			my $pos = $sam[3]; 
			my $end = $pos+length($sam[9]);
			my $seq = $sam[9];
			my $qual= $sam[10];  

		    # First read 	
			if ($flag & 0x40 && $flag & 0x2){

				$ReadName{$readname}{LEFT}{CHR} = $chr;
				$ReadName{$readname}{LEFT}{POS} = $pos;
				$ReadName{$readname}{LEFT}{END} = $end;
				$ReadName{$readname}{LEFT}{SEQ} = $seq;
				$ReadName{$readname}{LEFT}{QUAL}= $qual;
				$ReadName{$readname}{LEFT}{FLAG}= $flag; 
				$ReadName{$readname}{LEFT}{FIRST} = 1;
				$ReadName{$readname}{CHANGED} = 0; 

				if (exists $ReadName{$readname}{RIGHT}){

					if ($ReadName{$readname}{RIGHT}{POS} < 
						$ReadName{$readname}{LEFT}{POS}){
							$ReadName{$readname}{CHANGED} = 1; 

							$ReadName{$readname}{LEFT}{CHR} = $ReadName{$readname}{RIGHT}{CHR};
							$ReadName{$readname}{LEFT}{POS} = $ReadName{$readname}{RIGHT}{POS};
							$ReadName{$readname}{LEFT}{END} = $ReadName{$readname}{RIGHT}{END};
							$ReadName{$readname}{LEFT}{SEQ} = $ReadName{$readname}{RIGHT}{SEQ};
							$ReadName{$readname}{LEFT}{QUAL}= $ReadName{$readname}{RIGHT}{QUAL};
							$ReadName{$readname}{LEFT}{FLAG}= $ReadName{$readname}{RIGHT}{FLAG};
							$ReadName{$readname}{LEFT}{FIRST} = 0;

							$ReadName{$readname}{RIGHT}{CHR} = $chr;
							$ReadName{$readname}{RIGHT}{POS} = $pos;
							$ReadName{$readname}{RIGHT}{END} = $end;
							$ReadName{$readname}{RIGHT}{SEQ} = $seq;
							$ReadName{$readname}{RIGHT}{QUAL}= $qual;
							$ReadName{$readname}{RIGHT}{FLAG}= $flag;							
							$ReadName{$readname}{RIGHT}{FIRST} = 1;
					}	
				}
				$ReadName{$readname}{COUNT}++;
			} 
		    # Second read 	
			if ($flag & 0x80 && $flag & 0x2){

				$ReadName{$readname}{RIGHT}{CHR} = $chr;
				$ReadName{$readname}{RIGHT}{POS} = $pos;
				$ReadName{$readname}{RIGHT}{END} = $end;
				$ReadName{$readname}{RIGHT}{SEQ} = $seq;
				$ReadName{$readname}{RIGHT}{QUAL}= $qual;
				$ReadName{$readname}{RIGHT}{FLAG}= $flag;
				$ReadName{$readname}{CHANGED} = 0; 

				$ReadName{$readname}{RIGHT}{FIRST} = 0;
				if (exists $ReadName{$readname}{LEFT}){

					if ($ReadName{$readname}{RIGHT}{POS} < 
						$ReadName{$readname}{LEFT}{POS}){
							$ReadName{$readname}{CHANGED} = 1; 

							$ReadName{$readname}{RIGHT}{CHR} = $ReadName{$readname}{LEFT}{CHR};
							$ReadName{$readname}{RIGHT}{POS} = $ReadName{$readname}{LEFT}{POS};
							$ReadName{$readname}{RIGHT}{END} = $ReadName{$readname}{LEFT}{END};
							$ReadName{$readname}{RIGHT}{SEQ} = $ReadName{$readname}{LEFT}{SEQ};
							$ReadName{$readname}{RIGHT}{QUAL}= $ReadName{$readname}{LEFT}{QUAL};
							$ReadName{$readname}{RIGHT}{FLAG}= $ReadName{$readname}{LEFT}{FLAG};

							$ReadName{$readname}{RIGHT}{FIRST} = 1;

							$ReadName{$readname}{LEFT}{CHR} = $chr;
							$ReadName{$readname}{LEFT}{POS} = $pos;
							$ReadName{$readname}{LEFT}{END} = $end;
							$ReadName{$readname}{LEFT}{SEQ} = $seq;
							$ReadName{$readname}{LEFT}{QUAL}= $qual;
							$ReadName{$readname}{LEFT}{FLAG}= $flag;
							$ReadName{$readname}{LEFT}{FIRST} = 0;
					} 	 
				}
				$ReadName{$readname}{COUNT}++;
			}	
		}

		foreach my $read (sort keys %ReadName){

			$ReadName{$read}{CNVTYPE} = $cnvType; 
			# Both reads outisde the variant 
			if ($ReadName{$read}{LEFT}{END} <= $startCnv && $ReadName{$read}{RIGHT}{POS} >= $endCnv){ 
				$ReadName{$read}{CLASS} = "OUTSIDE";
				$ReadName{$read}{OLAP} = 0;
			} 

			if ($ReadName{$read}{LEFT}{END} <= $startCnv && $ReadName{$read}{RIGHT}{END} <= $startCnv){ 
				$ReadName{$read}{CLASS} = "OUTSIDE";
				$ReadName{$read}{OLAP} = 0;
			} 

			if ($ReadName{$read}{LEFT}{POS} >= $endCnv && $ReadName{$read}{RIGHT}{POS} >= $endCnv){ 
				$ReadName{$read}{CLASS} = "OUTSIDE";
				$ReadName{$read}{OLAP} = 0;
			}

			# Both reads within the variant 
			if ($ReadName{$read}{LEFT}{POS} >= $startCnv && $ReadName{$read}{LEFT}{END} <= $endCnv
				&& $ReadName{$read}{RIGHT}{POS} >= $startCnv && $ReadName{$read}{RIGHT}{END} <= $endCnv){
				$ReadName{$read}{CLASS} = "INSIDE";
				$ReadName{$read}{OLAP} = 0;
			} 
			# Left read outside, right read inside
			if ($ReadName{$read}{LEFT}{END} <= $startCnv && $ReadName{$read}{RIGHT}{POS} >= $startCnv
				&& $ReadName{$read}{RIGHT}{END} <= $endCnv ){

				$ReadName{$read}{CLASS} = "LEFT_OUTSIDE_RIGHT_INSIDE";
				$ReadName{$read}{OLAP} = 0;
			} 
			# Right read outside, left read inside
			if ($ReadName{$read}{RIGHT}{POS} >= $endCnv && $ReadName{$read}{LEFT}{POS} >= $startCnv
				&& $ReadName{$read}{LEFT}{END} <= $endCnv ){

				$ReadName{$read}{CLASS} = "LEFT_INSIDE_RIGHT_OUTSIDE";
				$ReadName{$read}{OLAP} = 0;
			} 

			# Left read overlapping 5', right read inside
			if ($ReadName{$read}{LEFT}{POS} < $startCnv && $ReadName{$read}{LEFT}{END} > $startCnv
				&& $ReadName{$read}{RIGHT}{POS} >= $startCnv && $ReadName{$read}{RIGHT}{END} <= $endCnv ) {
				$ReadName{$read}{CLASS} = "LEFT_OVERLAP_RIGHT_INSIDE";
				$ReadName{$read}{OLAP}  = 1;	
				$ReadName{$read}{"5PRIME"} = 1;
				$ReadName{$read}{"3PRIME"} = 0;
			} 

			# Left read inside, right read overlapping 3'
			if ($ReadName{$read}{LEFT}{POS} >= $startCnv && $ReadName{$read}{LEFT}{END} <= $endCnv
				&& $ReadName{$read}{RIGHT}{POS} < $endCnv && $ReadName{$read}{RIGHT}{END} > $endCnv ) {
				$ReadName{$read}{CLASS} = "LEFT_INSIDE_RIGHT_OVERLAP";
				$ReadName{$read}{OLAP}  = 1;	
				$ReadName{$read}{"5PRIME"} = 0;
				$ReadName{$read}{"3PRIME"} = 1;
			} 

			# Left overlapping 3', right inside 
			if ($ReadName{$read}{LEFT}{POS} < $endCnv && $ReadName{$read}{LEFT}{END} > $endCnv 
			   && $ReadName{$read}{RIGHT}{POS} >= $endCnv ){

				$ReadName{$read}{CLASS} = "LEFT_OVERLAP_RIGHT_OUTSIDE";
				$ReadName{$read}{OLAP} = 1;	
				$ReadName{$read}{"5PRIME"} = 0;
				$ReadName{$read}{"3PRIME"} = 1;				
			} 
			# Left outside, right olap 5' 
			if ($ReadName{$read}{LEFT}{END} <= $startCnv 
				&& $ReadName{$read}{RIGHT}{POS} <= $startCnv && $ReadName{$read}{RIGHT}{END} > $startCnv ){
				$ReadName{$read}{CLASS} = "LEFT_OUTSIDE_RIGHT_OVERLAP";
				$ReadName{$read}{OLAP} = 1;	
				$ReadName{$read}{"5PRIME"} = 1;
				$ReadName{$read}{"3PRIME"} = 0;
			} 

			# Left outside, right olap 3' 
			if ($ReadName{$read}{LEFT}{END} <= $startCnv 
				&& $ReadName{$read}{RIGHT}{POS} < $endCnv && $ReadName{$read}{RIGHT}{END} > $endCnv ){
				$ReadName{$read}{CLASS} = "LEFT_OUTSIDE_RIGHT_OVERLAP";
				$ReadName{$read}{OLAP} = 1;	
				$ReadName{$read}{"5PRIME"} = 0;
				$ReadName{$read}{"3PRIME"} = 1;
			} 

			# Right outside, left olap 3' 
			if ($ReadName{$read}{LEFT}{POS} < $startCnv && $ReadName{$read}{LEFT}{END} > $startCnv 
			   && $ReadName{$read}{RIGHT}{POS} >= $endCnv ){
				$ReadName{$read}{CLASS} = "LEFT_OVERLAP_RIGHT_OUTSIDE";
				$ReadName{$read}{OLAP} = 1;	
				$ReadName{$read}{"5PRIME"} = 1;
				$ReadName{$read}{"3PRIME"} = 0;
			} 

			# Both reads overlapping 3' 
			if ($ReadName{$read}{LEFT}{POS} < $startCnv && $ReadName{$read}{LEFT}{END} > $startCnv 
			   && $ReadName{$read}{RIGHT}{POS} < $startCnv && $ReadName{$read}{RIGHT}{END} > $startCnv){
				$ReadName{$read}{CLASS} = "LEFT_OVERLAP_RIGHT_OVERLAP";
				$ReadName{$read}{OLAP} = 1;	
				$ReadName{$read}{"5PRIME"} = 1;
				$ReadName{$read}{"3PRIME"} = 0;
				next;
			} 

			# Both reads overlapping 5' 
			if ($ReadName{$read}{LEFT}{POS} < $endCnv && $ReadName{$read}{LEFT}{END} > $endCnv 
			   && $ReadName{$read}{RIGHT}{POS} < $endCnv && $ReadName{$read}{RIGHT}{END} > $endCnv){

				$ReadName{$read}{CLASS} = "LEFT_OVERLAP_RIGHT_OVERLAP";
				$ReadName{$read}{OLAP} = 1;	
				$ReadName{$read}{"5PRIME"} = 0;
				$ReadName{$read}{"3PRIME"} = 1;
				next;
			} 
			# Both reads overlapping 5' and 3' breakpoints
			if ($ReadName{$read}{LEFT}{POS} < $startCnv && $ReadName{$read}{LEFT}{END} > $startCnv
				&& $ReadName{$read}{RIGHT}{POS} < $endCnv && $ReadName{$read}{RIGHT}{END} > $endCnv){
				$ReadName{$read}{CLASS}    = "LEFT_OVERLAP_RIGHT_OVERLAP";
				$ReadName{$read}{OLAP}     = 1;
				$ReadName{$read}{"5PRIME"} = 1;
				$ReadName{$read}{"3PRIME"} = 1;
				next;
			} 

		} 
		my $href = simulateCnv(\%ReadName, $chrCnv, $startCnv, $endCnv, $cnvType, $mode, $HoS{$bam}{READ_LENGTH}, $genome );
		%ReadName = %$href;
	}

	if (!-e $HoS{$bam}{FQ1} && !-e $HoS{$bam}{FQ1}){
		open (FQ1, ">", $HoS{$bam}{FQ1}) || die " ERROR: Unable to open $HoS{$bam}{FQ1}\n"; 
		open (FQ2, ">", $HoS{$bam}{FQ2}) || die " ERROR: Unable to open $HoS{$bam}{FQ2}\n"; 

		# Rewriting new FASTQ files
		my %seen = ();
		my $duplicatedCount = 0;
		open BAM,"$samtools view -h $bam |";
		while (my $line=<BAM>) {
			next if $line=~/^\@/;
			chomp $line;
			my @sam = split (/\t/, $line);
			my $rname = $sam[0] ; 
			my $flag = $sam[1];
			my $seq  = $flag & 0x10 ? revcomp($sam[9]) : $sam[9];
			my $qual = $sam[10];

			next if $sam[6] ne "=";
			if (!checkFlag($flag)){
				next;
			}
			if ($ReadName{$rname}{CLASS} ){
				my $cnvType = $ReadName{$rname}{CNVTYPE};

				# First read 
				if ($flag & 0x40 && $flag & 0x2){
					$ReadName{$rname}{LEFT}{SEQ} = $ReadName{$rname}{LEFT}{FLAG} & 0x10 
						? revcomp($ReadName{$rname}{LEFT}{SEQ}) : $ReadName{$rname}{LEFT}{SEQ};
					$qual = "E" x length($ReadName{$rname}{LEFT}{SEQ});
					print FQ1 "\@$rname\n$ReadName{$rname}{LEFT}{SEQ}\n+\n$qual\n";
					if ($cnvType =~/dup/i) {
						if ($ReadName{$rname}{DUPLICATED} ) {
							print FQ1 "\@dup$rname\n$ReadName{$rname}{LEFT}{SEQ}\n+\n$qual\n";
						}
					}   
				}
				# Second read 
				if ($flag & 0x80 && $flag & 0x2){
					$ReadName{$rname}{RIGHT}{SEQ} = $ReadName{$rname}{RIGHT}{FLAG} & 0x10 
						? revcomp($ReadName{$rname}{RIGHT}{SEQ}) : $ReadName{$rname}{RIGHT}{SEQ};
					$qual = "E" x length($ReadName{$rname}{RIGHT}{SEQ});
					print FQ2 "\@$rname\n$ReadName{$rname}{RIGHT}{SEQ}\n+\n$qual\n";
					if ($cnvType =~/dup/i) {
						if ($ReadName{$rname}{DUPLICATED}) {
							print FQ2 "\@dup$rname\n$ReadName{$rname}{RIGHT}{SEQ}\n+\n$qual\n";						
						} 
					}				
				}
			} 
			else {
				if ($flag & 0x40 && $flag & 0x2){
					print FQ1 "\@$rname\n$seq\n+\n$qual\n";
				}
				if ($flag & 0x80 && $flag & 0x2){
					print FQ2 "\@$rname\n$seq\n+\n$qual\n";
				}			
			} 
		}
		close BAM;
	} 

	my $sortedFq1 = $HoS{$bam}{FQ1};
	$sortedFq1 =~s/.fastq/.sorted.fastq/; 

	my $sortedFq2 = $HoS{$bam}{FQ2};
	$sortedFq2 =~s/.fastq/.sorted.fastq/; 
	my $cmd;
	if (!-e $sortedFq1 && $sortedFq2 ){
		$cmd = "$cat $HoS{$bam}{FQ1} | $paste - - - - | $sort -V -t \" \" | $tr \"\t\" \"\n\" > $sortedFq1";
		system $cmd;

		$cmd = "$cat $HoS{$bam}{FQ2} | $paste - - - - | $sort -V -t \" \" | $tr \"\t\" \"\n\" > $sortedFq2";
		system $cmd;
	} 

	$cmd  = "$bwa mem -t 7 -R \'\@RG\\tID:$HoS{$bam}{NAME}\\tSM:$HoS{$bam}{NAME}\' $genome $sortedFq1 $sortedFq2";
	$cmd .= " | $samtools sort -O BAM -o $outDir/$HoS{$bam}{NAME}.simulated.bam - $devNull";
	system $cmd;

	$cmd = " $samtools index $outDir/$HoS{$bam}{NAME}.simulated.bam";
	system $cmd;

	if (!$keepTmp){
		unlink($HoS{$bam}{FQ1}, $HoS{$bam}{FQ2}, $sortedFq1, $sortedFq2 );
	} 

	#exit;
# my $cmd = "$cat $HoS{$bam}{NAME}_1.fq | $paste - - - - | $sort -k1,1 -t \" \" | $tr \"\t\" \"\n\" > $HoS{$bam}{NAME}_1.sorted.fq";
# system $cmd;

# $cmd = "$cat $HoS{$bam}{NAME}_2.fq | $paste - - - - | $sort -k1,1 -t \" \" | $tr \"\t\" \"\n\" > $HoS{$bam}{NAME}_2.sorted.fq";
# system $cmd;

	#`bwa mem -t 4 ~/Programes/reference/bwa/ucsc.hg19.fasta $HoS{$bam}{NAME}_1.sorted.fq $HoS{$bam}{NAME}_2.sorted.fq | samtools view -Shu - | samtools sort  - second_test`;
	#`samtools index second_test.bam`;
}

##########################
sub generateCnvSequence {
	my $chrCnv  = shift;
	my $startCnv= shift;
	my $endCnv  = shift;
	my $cnvType = shift;
	my $readLength = shift;
	my $genomeFasta= shift;
	my $seq;

	if ($cnvType =~/del/i ){

		# For segment A 	
		my $startA = $startCnv-$readLength-1;
		my $endA   = $startCnv-1;
		my $seqA   = `$samtools faidx $genomeFasta $chrCnv:$startA-$endA | $grep -v '>' `;
		$seqA =~s/\n//g; 

		# For segment B
		my $startB = $endCnv;
		my $endB   = $endCnv+$readLength; 
		my $seqB   = `$samtools faidx $genomeFasta $chrCnv:$startB-$endB  | $grep -v '>'`;
		$seqB =~s/\n//g; 
		$seq = $seqA.$seqB;
	} 
	if ($cnvType =~/dup/i ){
		my $startA = $startCnv-$readLength-1;
		my $endA   = $startCnv-1;
		my $seqA   = `$samtools faidx $genomeFasta $chrCnv:$startA-$endA | $grep -v '>' `;
		$seqA =~s/\n//g; 

		# For segment B
		my $startB = $startCnv;
		my $endB   = $endCnv; 
		my $seqB   = `$samtools faidx $genomeFasta $chrCnv:$startB-$endB  | $grep -v '>'`;
		$seqB =~s/\n//g; 

		# For segment C
		my $startC = $endCnv;
		my $endC   = $endCnv+$readLength; 
		my $seqC   = `$samtools faidx $genomeFasta $chrCnv:$startC-$endC  | $grep -v '>'`;
		$seqC =~s/\n//g; 
		$seq = $seqA.$seqB.$seqB.$seqC;	
	} 
	return $seq;
} 

##########################
sub simulateCnv {
	my $href    = shift;
	my $chrCnv  = shift;
	my $startCnv= shift;
	my $endCnv  = shift;
	my $cnvType = shift;
	my $mode    = shift;
	my $readLength = shift;
	my $genomeFasta= shift;

	my $seq =  generateCnvSequence($chrCnv, $startCnv, $endCnv, $cnvType, $readLength, $genomeFasta);
	my $cnvLength = $endCnv-$startCnv;
	my %readHash = %$href;
	my $totalReads = 0;
	my $totalEdited = 0;
	
	foreach my $read (sort keys %readHash) {

		next if !$readHash{$read}{CLASS};  

		# Skipping reads that fall outise of the variant	
		next if $readHash{$read}{CLASS} eq 'OUTSIDE';  
		$totalReads++;

		# edit reads at a 50% chance
		my $edit = int(rand(1)+0.5);
		$totalEdited++ if $edit;
		my $duplicated = int(rand(1)+0.5);
		if ($edit){
			$readHash{$read}{EDITED} = 1; 
		}
		else{
			$readHash{$read}{EDITED} = 0; 
		} 
		if ($duplicated) {
			$readHash{$read}{DUPLICATED} = 1; 
		}
		else{
			$readHash{$read}{DUPLICATED} = 0; 
		} 
		if ($edit) {
			# For breakpoint non-overlapping paired reads
			if ($readHash{$read}{OLAP} == 0) {
				if ($cnvType =~/del/i) { 
					if ($readHash{$read}{CLASS} eq 'INSIDE'){
						$readHash{$read}{LEFT}{SEQ}  
							= substr($unmap, 0, length($readHash{$read}{LEFT}{SEQ}) ); 
						$readHash{$read}{RIGHT}{SEQ} 
							= substr($unmap, 0, length($readHash{$read}{RIGHT}{SEQ}) ); 
					} 
					elsif ($readHash{$read}{CLASS} eq 'LEFT_OUTSIDE_RIGHT_INSIDE'){
						$readHash{$read}{RIGHT}{SEQ} 
							= substr($unmap, 0, length($readHash{$read}{RIGHT}{SEQ}) ); 
					} 
					elsif ($readHash{$read}{CLASS} eq 'LEFT_INSIDE_RIGHT_OUTSIDE'){
						$readHash{$read}{LEFT}{SEQ}  
							= substr($unmap, 0, length($readHash{$read}{LEFT}{SEQ}) ); 
					}
				} 
				#print "INSIDE\n";
			} 
			# For breakpoint-overlapping paired reads
			if ($readHash{$read}{OLAP} == 1) {
				if ($readHash{$read}{CLASS} eq "LEFT_OVERLAP_RIGHT_INSIDE"){
					# Modify left read
					if ($cnvType =~/del/i) { 
						my $segaLength = $startCnv-$readHash{$read}{LEFT}{POS};
						my $segmentA   = substr($readHash{$read}{LEFT}{SEQ}, 0, $segaLength);

						my $segbLength = $readHash{$read}{LEFT}{END}-$startCnv; 
						my $segmentB   = substr($seq, $readLength, $segbLength);

						my $newSeq = $segmentA . $segmentB;
						$readHash{$read}{LEFT}{SEQ} = $newSeq;
						$readHash{$read}{RIGHT}{SEQ}
							= substr($unmap, 0, length($readHash{$read}{RIGHT}{SEQ}));
					}
					if ($cnvType =~/dup/i){
						my $segaLength = $startCnv-$readHash{$read}{LEFT}{POS};
						my $segmentA = substr($seq, $readLength+$cnvLength-$segaLength, $segaLength);

						my $segbLength = $readHash{$read}{LEFT}{END}-$startCnv; 
						my $segmentB = substr($readHash{$read}{LEFT}{SEQ}, $segaLength, $segbLength);

						my $newSeq = $segmentA . $segmentB;
						$readHash{$read}{LEFT}{SEQ} = $newSeq;
						print "1\t$newSeq\n";			
					} 
				} 
				if ($readHash{$read}{CLASS} eq "LEFT_INSIDE_RIGHT_OVERLAP"){
					print "LEFT_INSIDE_RIGHT_OVERLAP\n";
					print "$cnvType\n";
					if ($cnvType =~/del/i) { 
						# Modify right read
						my $segaLength = $endCnv-$readHash{$read}{RIGHT}{POS};
						my $segmentA = reverse(substr(reverse($seq), $readLength, $segaLength));

						my $segbLength = $readHash{$read}{RIGHT}{END}-$endCnv; 
						my $segmentB = substr($readHash{$read}{RIGHT}{SEQ}, $segaLength+1, $segbLength);
						my $newSeq = $segmentA . $segmentB;
							print "$newSeq\n";
						$readHash{$read}{RIGHT}{SEQ} = $newSeq;
						$readHash{$read}{LEFT}{SEQ}
							= substr($unmap, 0, length($readHash{$read}{LEFT}{SEQ}));
					}
					if ($cnvType =~/dup/i) { 
						# Modify right read
						my $segaLength = $endCnv-$readHash{$read}{RIGHT}{POS};
						my $segmentA = substr($readHash{$read}{RIGHT}{SEQ}, 0, $segaLength);
						#my $segmentA = reverse(substr(reverse($seq), $readLength, $segaLength));

						my $segbLength = $readHash{$read}{RIGHT}{END}-$endCnv; 
						my $segmentB = substr($seq, $readLength, $segbLength);
						my $newSeq = $segmentA . $segmentB;
						$readHash{$read}{RIGHT}{SEQ} = $newSeq;
						print "2\t$newSeq\n";			

					} 
				}
				if ($readHash{$read}{CLASS} eq "LEFT_OUTSIDE_RIGHT_OVERLAP"){

					if ($readHash{$read}{"5PRIME"}){ 

						if ($cnvType =~/del/i) {
							my $segaLength = $startCnv-$readHash{$read}{RIGHT}{POS};
							my $segmentA = substr($readHash{$read}{RIGHT}{SEQ}, 0, $segaLength);
							
							my $segbLength = $readHash{$read}{RIGHT}{END}-$startCnv; 
							my $segmentB = substr($seq, $readLength+1, $segbLength);
							my $newSeq = $segmentA . $segmentB;
							$readHash{$read}{RIGHT}{SEQ} = $newSeq;
						} 
						if ($cnvType =~/dup/i) {
							my $segaLength = $startCnv-$readHash{$read}{RIGHT}{POS};
							my $segmentA = substr($seq, $readLength+$cnvLength-$segaLength, $segaLength);

							my $segbLength = $readHash{$read}{RIGHT}{END}-$startCnv; 
							my $segmentB = substr($readHash{$read}{RIGHT}{SEQ}, $segaLength, $segbLength);

							my $newSeq = $segmentA . $segmentB;
							$readHash{$read}{RIGHT}{SEQ} = $newSeq;
							print "3\t$newSeq\n";		
						} 
					}
					if ($readHash{$read}{"3PRIME"}){
						if ($cnvType =~/del/i) {
							my $segaLength = $endCnv-$readHash{$read}{RIGHT}{POS};
							my $segmentA = reverse(substr(reverse($seq), $readLength, $segaLength));
							#print length($readHash{$read}{RIGHT}{SEQ}) . "\n";
							
							my $segbLength = $readHash{$read}{RIGHT}{END}-$endCnv; 
							my $segmentB = substr($readHash{$read}{RIGHT}{SEQ}, $segaLength+1, $segbLength);
							my $newSeq = $segmentA . $segmentB;
							$readHash{$read}{RIGHT}{SEQ} = $newSeq;
						}
						if ($cnvType =~/dup/i) {
							# Modify right read
							my $segaLength = $endCnv-$readHash{$read}{RIGHT}{POS};
							my $segmentA = substr($readHash{$read}{RIGHT}{SEQ}, 0, $segaLength);

							my $segbLength = $readHash{$read}{RIGHT}{END}-$endCnv; 
							my $segmentB = substr($seq, $readLength, $segbLength);
							my $newSeq = $segmentA . $segmentB;
							$readHash{$read}{RIGHT}{SEQ} = $newSeq;
							print "4\t$newSeq\n";			
						}	
					 	#print "$newSeq\n";
					} 
				# 	# print length($newSeq). "\n\n";
				}	
				if ($readHash{$read}{CLASS} eq "LEFT_OVERLAP_RIGHT_OUTSIDE"){

					if ($readHash{$read}{"5PRIME"}){
						if ($cnvType =~/del/i) { 
							my $segaLength = $startCnv-$readHash{$read}{LEFT}{POS};
							my $segmentA = substr($readHash{$read}{LEFT}{SEQ}, 0, $segaLength);
							#print length($readHash{$read}{RIGHT}{SEQ}) . "\n";
							
							my $segbLength = $readHash{$read}{LEFT}{END}-$startCnv; 
							my $segmentB = substr($seq, $readLength+1, $segbLength);
							my $newSeq = $segmentA . $segmentB;
							$readHash{$read}{LEFT}{SEQ} = $newSeq;
							#print "$newSeq\n";
						}
						if ($cnvType =~/dup/i) { 
							my $segaLength = $startCnv-$readHash{$read}{LEFT}{POS};
							my $segmentA = substr($seq, $readLength+$cnvLength-$segaLength, $segaLength);

							my $segbLength = $readHash{$read}{LEFT}{END}-$startCnv; 
							my $segmentB = substr($readHash{$read}{LEFT}{SEQ}, $segaLength, $segbLength);

							my $newSeq = $segmentA . $segmentB;
							$readHash{$read}{LEFT}{SEQ} = $newSeq;
							print "5\t$newSeq\n";	
						}	
					}	
					if ($readHash{$read}{"3PRIME"}){
						if ($cnvType =~/del/i) { 
							# Modify right read
							my $segaLength = $endCnv-$readHash{$read}{LEFT}{POS};
							my $segmentA = reverse(substr(reverse($seq), $readLength, $segaLength));

							my $segbLength = $readHash{$read}{LEFT}{END}-$endCnv; 
							my $segmentB = substr($readHash{$read}{LEFT}{SEQ}, $segaLength, $segbLength);
							my $newSeq = $segmentA . $segmentB;
							$readHash{$read}{LEFT}{SEQ} = $newSeq;
							#print "$newSeq\n";
						}
						if ($cnvType =~/dup/i) {
							my $segaLength = $endCnv-$readHash{$read}{LEFT}{POS};
							my $segmentA = substr($readHash{$read}{LEFT}{SEQ}, 0, $segaLength);

							my $segbLength = $readHash{$read}{LEFT}{END}-$endCnv; 
							my $segmentB = substr($seq, $readLength, $segbLength);
							my $newSeq = $segmentA . $segmentB;
							$readHash{$read}{LEFT}{SEQ} = $newSeq;
							print "6\t$newSeq\n";	
						} 
					}	
				}
				if ($readHash{$read}{CLASS} eq "LEFT_OVERLAP_RIGHT_OVERLAP"){
					if ($readHash{$read}{"5PRIME"} && !$readHash{$read}{"3PRIME"} ){
						if ($cnvType =~/del/i) {
							my $segaLength = $startCnv-$readHash{$read}{LEFT}{POS};
							my $segmentA = substr($readHash{$read}{LEFT}{SEQ}, 0, $segaLength);
						
							my $segbLength = $readHash{$read}{LEFT}{END}-$startCnv; 
							my $segmentB = substr($seq, $readLength+1, $segbLength);
							my $newSeq = $segmentA .  $segmentB;
							$readHash{$read}{LEFT}{SEQ} = $newSeq;

							$segaLength = $startCnv-$readHash{$read}{RIGHT}{POS};
							$segmentA = substr($readHash{$read}{RIGHT}{SEQ}, 0, $segaLength);

							$segbLength = $readHash{$read}{RIGHT}{END} - $startCnv; 
							$segmentB = substr($seq, $readLength+1, $segbLength);			
							$newSeq = $segmentA . $segmentB;
							$readHash{$read}{RIGHT}{SEQ} = $newSeq;
						}
						if ($cnvType =~/dup/i) {
							# my $segaLength = $startCnv-$readHash{$read}{LEFT}{POS};
							# my $segmentA = substr($seq, $readLength+$cnvLength-$segaLength, $segaLength);

							# my $segbLength = $readHash{$read}{LEFT}{END}-$startCnv; 
							# my $segmentB = substr($readHash{$read}{LEFT}{SEQ}, $segaLength, $segbLength);

							# my $newSeq = $segmentA . $segmentB;
							# $readHash{$read}{LEFT}{SEQ} = $newSeq;
							# print "7 $newSeq\n";

							# $segaLength = $endCnv-$readHash{$read}{RIGHT}{POS};
							# $segmentA = substr($readHash{$read}{RIGHT}{SEQ}, 0, $segaLength);

							# $segbLength = $readHash{$read}{RIGHT}{END}-$endCnv; 
							# $segmentB = substr($seq, $readLength, $segbLength);

							# $newSeq = $segmentA . $segmentB;
							# $readHash{$read}{RIGHT}{SEQ} = $newSeq;
							# 						print "8  $readHash{$read}{RIGHT}{END} $endCnv $newSeq\n";

						} 
					}
					if ($readHash{$read}{"3PRIME"} && !$readHash{$read}{"5PRIME"} ){
						if ($cnvType =~/del/i) {
							my $segaLength = $endCnv-$readHash{$read}{LEFT}{POS};
							my $segmentA = reverse(substr(reverse($seq), $readLength, $segaLength));
							#print length($readHash{$read}{RIGHT}{SEQ}) . "\n";
							
							my $segbLength = $readHash{$read}{LEFT}{END}-$endCnv; 
							my $segmentB = substr($seq, $readLength, $segbLength);
							my $newSeq = $segmentA . $segmentB;
							$readHash{$read}{LEFT}{SEQ} = $newSeq;
							#print "$newSeq\n";

							$segaLength = $endCnv-$readHash{$read}{RIGHT}{POS};
							$segmentA = reverse(substr(reverse($seq), $readLength, $segaLength));

							$segbLength = $readHash{$read}{RIGHT}{END}-$endCnv; 
							$segmentB = substr($seq, $readLength, $segbLength);
							$newSeq = $segmentA . $segmentB;
							$readHash{$read}{RIGHT}{SEQ} = $newSeq;
						} 
					}
					if ($readHash{$read}{"5PRIME"} && $readHash{$read}{"3PRIME"} ){
						if ($cnvType =~/del/i) {  
							# my $segaLength;
							# my $segmentA;
							# my $segbLength;
							# my $segmentB;
							# my $newSeq;
							# my $seq =  generateCnvSequence($chrCnv, $startCnv, $endCnv, $cnvType, $readLength, $genomeFasta);
							my $lenA = length($readHash{$read}{LEFT}{SEQ});
							my $lenB = length($readHash{$read}{RIGHT}{SEQ});

							my $segaLength = $startCnv-$readHash{$read}{LEFT}{POS};
							my $segmentA = substr($readHash{$read}{LEFT}{SEQ}, 0, $segaLength);
							my $segbLength = $readHash{$read}{LEFT}{END}-$startCnv; 
							my $segmentB = substr($seq, $readLength+1, $segbLength);
							my $newSeq = $segmentA . $segmentB;
							$readHash{$read}{LEFT}{SEQ} = $newSeq;
							#print "left $readHash{$read}{LEFT}{SEQ}\t$newSeq\n\n"; 

							my $sizeSegA1 = $segaLength;
							my $sizeSegB1 = $segbLength;

							$segaLength = $endCnv-$readHash{$read}{RIGHT}{POS};
							$segmentA = reverse(substr( reverse($seq), $readLength, $segaLength));
							$segbLength = $readHash{$read}{RIGHT}{END}-$endCnv; 
							$segmentB = substr($readHash{$read}{RIGHT}{SEQ}, $segaLength+1, $segbLength);
							$newSeq = $segmentA . $segmentB;
							$readHash{$read}{RIGHT}{SEQ} = $newSeq;
							#print "right $readHash{$read}{RIGHT}{SEQ}\t$newSeq\n\n";

							my $sizeSegA2 = $segaLength;
							my $sizeSegB2 = $segbLength;

							if ($sizeSegA2 >= $sizeSegA1 or $sizeSegB1 >= $sizeSegB2){

								my $segLeft = $readHash{$read}{LEFT}{SEQ};
								my $segRight= $readHash{$read}{RIGHT}{SEQ};
								$readHash{$read}{LEFT}{SEQ} = $segRight;
								$readHash{$read}{RIGHT}{SEQ} = $segLeft;
							}
						}  
					}	
				} 
			}
		}	
	} 
	print "$totalReads\t$totalEdited\n";

	return\%readHash; 
} 


##############
sub revcomp {
my $seq = shift;
my $revcomp = reverse $seq;
$revcomp =~ tr/ATGCatgc/TACGtacg/;
return $revcomp;
}


##############
sub checkFlag {
	my $flag = shift;
	my $value = 1;
	if ($flag & 0x4){
		$value = 0;
	}
	if ($flag & 0x8){
		$value = 0;
	}
	if ($flag & 0x100){
		$value = 0;
	}
	if ($flag & 0x200){
		$value = 0;
	}
	if ($flag & 0x400){
		$value = 0;
	}
	if ($flag & 0x800){
		$value = 0;
	}
	return $value;
}


##############
sub getSampleInfo {

	my $config = shift;
	my $outDir = shift;

	my %HoS = ();

	my $str = `$cat $config`; chomp $str;
	my @tmpStr = split (/\n/, $str);

	foreach my $line (@tmpStr) {
		next if $line =~/^#/;   
		my @tmp = split (/\t/, $line);
		my $bam = $tmp[0];
		if (!-e $bam) {
			print " ERROR: $bam does not exist\n";exit;
		}
		if (!-e "$bam.bai") {
			print " INFO: Indexing $bam\n";
			my $cmd = "$samtools index $bam";
			system $cmd;
		}

		#my $metrics=`$samtools view -q 20 -f2 $bam | $cut -f9 | $head -1000| $awk '{if (\$1<0){\$1=-\$1}else{\$1=\$1} sum+=\$1; sumsq+=\$1*\$1} END {print sum/NR, sqrt(sumsq/NR - (sum/NR)^2)}'`;
		#my ($mean, $stdev)= split (/\s/,$metrics);
		my $header = `$samtools view -H $bam`; chomp $header;

		$HoS{$bam}{NAME}       = basename($bam);
		$HoS{$bam}{NAME}       =~s/.bam//;   
		$HoS{$bam}{HEADER}     = $header;
		#$HoS{$bam}{MEAN_ISIZE} = $mean;
		#$HoS{$bam}{STD_ISIZE}  = $stdev;
		$HoS{$bam}{READ_LENGTH}= `$samtools view $bam | $awk '{print length(\$10)}' | $head -1000 | $sort -u | $tail -1`; 
		chomp $HoS{$bam}{READ_LENGTH};
		$HoS{$bam}{FQ1}        = "$outDir/$HoS{$bam}{NAME}" . "_R1_simulated.fastq";
		$HoS{$bam}{FQ2}        = "$outDir/$HoS{$bam}{NAME}" . "_R2_simulated.fastq";
		my $str = `$grep "$HoS{$bam}{NAME}.bam" $cnv_list`;
		if (!$str) {
			print "WARNING: bam $HoS{$bam}{NAME} has no CNV associated!\n";
			$HoS{$bam}{CNV} = "none";			
			next;
		}
		my @arr = split (/\n/, $str);
		my @arrCNV = ();
		foreach my $item (@arr) {
			my @tmp = split (/\t/, $item);
			my $cnv = "$tmp[1]\t$tmp[2]\t$tmp[3]\t$tmp[4]\t$tmp[5]\t$tmp[6]\t$tmp[7]";
			$HoS{$bam}{"$tmp[1]\t$tmp[2]\t$tmp[3]"}{PLOIDY} = $tmp[6];
			$HoS{$bam}{"$tmp[1]\t$tmp[2]\t$tmp[3]"}{MODE} = $tmp[7];

			push @arrCNV, $cnv;
		}
		@arrCNV = natsort @arrCNV;
		$HoS{$bam}{CNV}   =   \@arrCNV;
	}
	return %HoS;
}

##############
sub Help {
print "\n Usage: perl $0 <options>

 Options:
 -b,--bed      STRING   ROI bed
 -g,--genome   STRING   Genome file in FASTA format
 -c,--config   STRING   Config file with a CNV list
 -o,--outdir   STRING   Output directory\n\n";
 exit;
}

