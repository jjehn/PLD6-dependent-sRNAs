# intention of this script: remove sequences bearing TGTGATGCATGTA (part of PLD6-siRNA-sequence #3) from collapsed and dusted FASTQ file

######################################################################################################################################
# scipt was adapted from TBr2_duster.pl (NGStoolbox, small RNA group Uni Mainz) Version: 2.1 LAST MODIFIED: 13. March 2015
######################################################################################################################################

#!/usr/bin/perl
use Getopt::Long;
$|=1;

GetOptions
	(
	"i=s"=>\$input_file,
	);

# check input format (FASTA/FASTQ)
print"\nChecking input file format...";
open(IN,$input_file)||die print"\nCould not open input file $input_file.\n$!\n\n";
$i=0;
$fasta_heads=0;
$fastq_heads=0;
while(<IN>)
	{
	$i++;
	last if($i==1000);
	if($_=~/^>/)
		{
		$fasta_heads++;
		}
	elsif($_=~/^@/||$_=~/^\+/)
		{
		$fastq_heads++;
		}
	}
close IN;
if($fasta_heads>$fastq_heads)
	{
	$format='FASTA';
	}
elsif($fastq_heads>$fasta_heads)
	{
	$format='FASTQ';
	}
else
	{
	die print"\nUnknown file format of input file $input_file.\n\n";
	}
print" done. -> $format";



# process input file
print"\nProcess sequences in $input_file...";
open(IN,$input_file)||die print"\nCould not open input file $input_file.\n$!\n\n";
open(OK,">$input_file.non_sirna")||die print"\nCould not create output file $input_file.non_sirna.\n$!\n\n"; # OK means non-sirna
open(SIRNA,">$input_file.sirna")||die print"\nCould not create output file $input_file.sirna.\n$!\n\n"; # DUST is SIRNA

$i=0;
$sirna=0; # simple is sirna
$non_sirna=0; # complex is non_sirna
$reads_sirna=0;
$reads_non_sirna=0;
@seq_data=();
while(<IN>)
	{
	$_=~s/\s*$//;
	push(@seq_data,$_);
	if($format eq'FASTQ')
		{
		$i+=0.25;
		}
	elsif($format eq'FASTA')
		{
		$i+=0.5;
		}
	if($i==1)
		{
		$i=0;
		$ok=1;
		$fragment="TGTGATGCATGTA"; # exchange sequence here, if other siRNA was transfected
		$check=$seq_data[1];
		
		$readcount = $seq_data[0]; # get rid of @ or > before readcount
		@readcount = split('',$readcount);
		shift @readcount;
		$reads = join('',@readcount);
		
		
		if($check=~s/$fragment//g)
			{
			$ok=0; # goes into siRNA
			}
			
		if($format eq'FASTQ')
			{
			$output="$seq_data[0]\n$seq_data[1]\n$seq_data[2]\n$seq_data[3]\n";
			}
		else
			{
			$output="$seq_data[0]\n$seq_data[1]\n";
			}
		if($ok==1)
			{
			$non_sirna++;
			print OK$output;
			$reads_non_sirna=($reads_non_sirna+$reads);
			}
		else
			{
			$sirna++;
			print SIRNA$output;
			$reads_sirna=($reads_sirna+$reads);
			}
		@seq_data=();
		$readcount=0;
		$reads=0;
		@readcount=();
		}
	}

$total_reads=$reads_non_sirna+$reads_sirna;
$perc_reads_non_sirna=int(100*$reads_non_sirna/$total_reads);
$perc_reads_sirna=int(100*$reads_sirna/$total_reads);
print"\nFound $sirna siRNA-derived (anti-PLD6 #3) non-identical sequences (Their reads account for $perc_reads_sirna % of the input reads).\n$non_sirna non-identical sequences are not related to the transfected siRNA ($perc_reads_non_sirna % of total reads).\n\n";
exit;