#!/usr/bin/perl --

use Getopt::Long qw(:config bundling auto_help);
use Pod::Usage;

=head1 NAME

samSEFilter.pl

=head1 AUTHOR

Cheng Jia
cheng.jia@outlook.com

=head1 SYNOPSIS

samSEFilter.pl 

This script checks the validity of the input .sam file under the assumption that the fragments are single-end sequencing results mapped to the chromosome.

It reads in the .sam file line by line, and outputs the reads that pass the standard specified by the user.

The script employs the following criteria:

1. Reads should pass QC; (this info is extracted from the bitwise FLAG of .sam file)

2. MAPQ score should be higher than m;

3. Reads should be uniquely mapped to the genome.

Parameters:

=over 4

-V, --verbose	Output extra warnings and messages. (Default:no verbose);

-M,	--mapqscore	MAPQ score for both reads should be equal to or higher than m; (Default=40)

-A, --accepted	filename for output of the accepted reads; (Default=samfile_accepted.sam)

-D,	--discarded	filename for output of the discarded reads;(Default:do not output discarded reads unless this name is indicated.)
				
-H,	--help		Display help message.

Warning: This program is free software. It comes without warranty, to the extent permitted by applicable law. 
Use it at your own risk. The authors and contributors are NOT responsible for any negative situations,
including but not limited to physical, psychological, pecuniary, social, emotional, and romantic losses, incurred
by usage of this program. When in doubt, refer to http://sam.zoy.org/wtfpl/COPYING for licensing information.

=back

=cut

#TODO: Add a flag --sortreads to sort the reads before piping them into this program.

#Process command line input;
$mapqscore=30;
$accepted="samfile_accepted.sam";
$discarded="";
$help=0;
$verbose=0;

GetOptions(
	'mapqscore|M=i'			=>	\$mapqscore,	
	'accepted|A=s'		=>	\$accepted,
	'discarded|D=s'		=>	\$discarded,
	'verbose|V!'		=>	\$verbose,
	'help|H!'				=>	\$help

);

if($help){
my $podHandle=\*STDERR;
pod2usage(-output=>$podHandle);
return;
}

$totalLine=0;
$numOutput=0;

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

$totalLine=0;
$numOutput=0;
while(<SAMIN>){
	$totalLine++;
	chomp;
	my $cur_read=$_;
	my @cur_line=split("\t",$cur_read);
	
	my $cur_qcflag=1-($cur_line[1] & 0x200);

	if($cur_line[14] =~ /.*NH\:i:(\d+).*/){
										if($1 eq "1"){
										$cur_uniqflag=1;
									}else{
										$cur_uniqflag=0;
									}
	my $cur_mapq=$cur_line[4];

	
	if(($cur_qcflag == 1)&&($cur_uniqflag == 1)&&($cur_mapq>=$mapqscore)){
		$numOutput++;
		print ACPTOUT $cur_read."\n";		
	}	else
	{
		unless($discarded eq "")
		{
		print DISDOUT $cur_read."\n";}
	}

}}

close (SAMIN);
close (ACPTOUT);
unless($discarded eq ""){
close (DISDOUT);}

print STDERR "Processed $totalLine lines\n";
print STDERR "Output $numOutput lines\n";
