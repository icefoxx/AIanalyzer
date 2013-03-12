#!/usr/bin/perl --

use Getopt::Long qw(:config bundling auto_help);
use Pod::Usage;

=head1 NAME

samPEFilter.pl

=head1 AUTHOR

Cheng Jia
cheng.jia@outlook.com

=head1 SYNOPSIS

samFilter.pl 

This script checks the validity of the input .sam file under the assumption that the fragments are paired-end sequencing results mapped to the chromosome.

It reads blocks of fragments with the same RNAME, and perform pairwise comparison between these fragments to determine the correctly mapped pairs.

If you want to use it on BOWTIE mapped .sam file, sort the .sam file first to allow reading of grouped reads with the same RNAME.

The script employs the following criteria:

1. paired reads should be mapped onto the same chr;

2. paired reads should be mapped onto different strands;

3. paired reads should be within n bp of each other;

4. at least 0/1/2 read(s) of the pair should pass QC; (this info is extracted from the bitwise FLAG of .sam file)

5. MAPQ score should be higher than m;

6. At least 0/1/2 reads(s) of the pair should be a unique hit in the genome; (this info is extracted from the optional field of .sam file)

In the case that more than 2 possible pairs of reads left after filtering, they will be deemed non-uniquely mapped reads and discarded. 

Parameters:

=over 4

-V, --verbose	Output extra warnings and messages. (Default:no verbose);

-N, --distance	Largest distance separating a pair of reads; (This parameter is inclusive, which means that it contains the reads' lengths themselves. Default=200000)

-Q, --passqc	this parameter can take 3 values (Default=1)
				0-no restrictions on QC status of the PE reads;
				1-at least 1 read from the pair should pass QC;
				2-both reads from the pair should pass QC;

-M,	--mapqscore	MAPQ score for both reads should be equal to or higher than m; (Default=40)

-U, --uniquehit	this parameter can take 3 values (Default=1)
				0-do not consider the uniqueness of mapping;
				1-at least 1 read from the pair should be uniquely mapped to the chr;
				2-both reads from the pair should be uniquely mapped to the chr;

-A, --accepted	filename for output of the accepted reads; (Default=samfile_accepted.sam)

-D,	--discarded	filename for output of the discarded reads;(Default:do not output discarded reads unless this name is indicated.)
				
-H,	--help		Display help message.

Warning: This program is free software. It comes without warranty, to the extent permitted by applicable law. 
Use it at your own risk. The authors and contributors are NOT responsible for any negative situations,
including but not limited to physical, psychological, pecuniary, social, emotional, and romantic losses, incurred
by usage of this program. When in doubt, refer to http://sam.zoy.org/wtfpl/COPYING for licensing information.

=back

=cut


#TODO: Add a flag --onlybestpair to output only one pair per RNAME group 
#right now, since it can output multiple pairs of results, so if you run this on a filtered .sam file.
#you get an exponetially larger file. I need to figure out what's happening there. Could be permanent loop?
#this phenomenon disappears once I add command to break out the loop immediately after discovering the FIRST CORRECT pair.
#TODO: Add a flag --sortreads to sort the reads before piping them into this program.

#Process command line input;
$distance=200000;
$passqc=1;
$mapqscore=40;
$uniquehit=1;
$accepted="samfile_accepted.sam";
$discarded="";
$help=0;
$verbose=0;

GetOptions(
	'distance|N=i'	=>	\$distance,		
	'passqc|Q=i'			=>	\$passqc,		
	'mapqscore|M=i'			=>	\$mapqscore,	
	'uniquehit|U=i'		=>	\$uniquehit,	
	'accepted|A=s'		=>	\$accepted,
	'discarded|D=s'		=>	\$discarded,
	'verbose|v!'		=>	\$verbose,
	'help|H!'				=>	\$help

);

if($help){
my $podHandle=\*STDERR;
pod2usage(-output=>$podHandle);
return;
}

$totalLine=0;
$numOutput=0;


sub outputLines{
	
	my @lineOut=@_;
	
	if (scalar(@lineOut)==2){
		print ACPTOUT @lineOut[0];
		print ACPTOUT @lineOut[1];

		$numOutput+=2;
	}else
	{
		
	unless($discarded eq ""){
		for(my $linenum=0;$linenum<scalar(@lineOut);$linenum++){
				print DISDOUT @lineOut[$linenum];
				
				}}}}


sub pairwiseComp{
	
	
	
	my @block=@_;
	
	my @outBlock;
	if(scalar(@block) <2){
		#DISCARD IF A BLOCK HAS ONLY ONE ELEMENT;
		unless($discarded eq ""){
				print DISDOUT $block[0]."\n";}
	
	}else{
		#OUTERLOOP:#
		 for(my $i=0;$i<scalar(@block);$i++){
			for(my $j=$i+1;$j<scalar(@block);$j++){
				
				my $readOne=$block[$i];
				my $readTwo=$block[$j];
				my @lineOne=split("\t",$readOne);
				my @lineTwo=split("\t",$readTwo);
				
				if($lineOne[1] & 0x10){$readOneStrand="-";}
					else{$readOneStrand="+";}
				if($lineOne[1] & 0x20){$readOneNextStrand="-";}
					else{$readOneNextStrand="+";}
				if($lineTwo[1] & 0x10){$readTwoStrand="-";}
					else{$readTwoStrand="+";}
				if($lineTwo[1] & 0x20){$readTwoNextStrand="-";}
					else{$readTwoNextStrand="+";}
					
					if($verbose){
				if ($readOneNextStrand ne $readTwoStrand){print STDERR "mate negative strand flag does not match read negative strand flag of mate.\nUse own strand info from flag.\n"}
				if (($lineOne[3] != $lineTwo[7]) or ($lineOne[7] != $lineTwo[3])){
					print STDERR "mate position does not match position info of the mate.\nUse own position info.\n";
					}}
					
				my $readOnePosi=$lineOne[3];
				my $readTwoPosi=$lineTwo[3];
				my $passuniqhits=0;
				my $passmapqtest=0;
				my $passqctest=0;
				my $samechr=0;
				my $diffstrand=0;
				my $closeenough=0;
				
				
				my $correct=0;
######################START PAIRWISE COMPARISON
				if($lineOne[2] eq $lineTwo[2]){		#whether two PE reads are on the same Chr
					$samechr=1;
					if($readOneStrand ne $readTwoStrand){	#whether two PE reads are on the different strand
						$diffstrand=1;
						$readdist=abs($readOnePosi-$readTwoPosi);
						if($readdist<=$distance){		#whether two PE reads are closer than $distance bp.
							$closeenough=1;
							$numberpassqc=(1-($lineOne[1] & 0x200))+(1-($lineTwo[1] & 0x200));
							if($numberpassqc>=$passqc){	#whether sufficient number of reads in a pair pass QC;
								$passqctest=1;
								if(($lineOne[4]>=$mapqscore)&&($lineTwo[4]>=$mapqscore)){ #whether both reads have MAPQ score larger than $mapqscore
									$passmapqtest=1;
									if($lineOne[11] =~ /.*NH\:i:(\d+).*/){
										if($1 eq "1"){
											$uniqueOne=1;
											}else{
											$uniqueOne=0;
										}
									}
									if($lineTwo[11] =~ /.*NH\:i:(\d+).*/){
										if($1 eq "1"){
										$uniqueTwo=1;
									}else{
										$uniqueTwo=0;
									}
									}
										my $uniquetotal=$uniqueOne+$uniqueTwo;
										if($uniquehit>=$uniquetotal){	#whether number of unique hits of the pair are larger than or equal to $uniquehit
										$passuniqhits=1;
										
										push(@outBlock, $readOne."\n");
										push(@outBlock, $readTwo."\n");
										
										
								
							}}}}}
				
				
	
				}}
		 
		 }}
	
unless(scalar(@outBlock)==0){
			outputLines(@outBlock);}
	
}


if((scalar(@ARGV) == 0) ||
        ((scalar(@ARGV) == 1) && (@ARGV[0] eq "-"))) {
        $samInCMD = "<-";
} else {
	
		$listARGS = join (' ', @ARGV);
		$samInCMD = "cat $listARGS | ";
}


open(SAMIN,$samInCMD)
or die("Cannot open file\(s\) $samInCMD.\n");

open ACPTOUT, ">", $accepted or die $!;
unless($discarded eq ""){
open DISDOUT, ">", $discarded or die $!;}
while(<SAMIN>){
unless ($_ =~ m/^\@SQ/){
	$previousRead=$_;
	last;
}}
close(SAMIN);

open(SAMIN,$samInCMD)
or die("Cannot open file\(s\) $samInCMD.\n");

while(<SAMIN>){
	chomp;
	unless ($_ =~ m/^\@SQ/){
	$currentRead=$_;}else{next;}
	@currentLine=split("\t",$currentRead);
	@previousLine=split("\t",$previousRead);
	
	unless($currentLine[0] eq $previousLine[0]){
	pairwiseComp(@currentBlock);
	undef @currentBlock;	
	}
	
	
	$totalLine++;
	push(@currentBlock, $currentRead);
	$previousRead=$currentRead;
	}
	
	pairwiseComp(@currentBlock);

close (SAMIN);
close (ACPTOUT);
unless($discarded eq ""){
close (DISDOUT);}

print STDERR "Processed $totalLine lines\n";
print STDERR "Output $numOutput lines\n";
