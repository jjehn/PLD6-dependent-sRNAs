use Getopt::Long;
use Iterator::Array::Jagged;
$|=1;
$date="2018-05-29";
$version="0.0.0";

### USAGE / INFORMATION ###
$logo="


                          voXZMOBBBBGZBBBkvriiiiii iuM                         
                       kBBBBNILiiiirYLvISBMOOOBBBBBBBU                         
                    iBBMFvi       iriYoLvrUPSPXGBBBB                           
                   MBBvi       i iiiiivLYvLYJuJjVMBj                           
                 IBBkrLriii  iiiiii iiivLviiYvJuuLXv                           
                BBBuvSFkYiivvvLivYrirrii iiiiirriivFUv                         
              MBBMoriii iuXGMSuuiLYLLIUvvijri iiiirruXOEki            vOBB     
           LMBBBViii    iuFXUoFOvXYNFLriiuiuiiiiiiLivvjooqENNBF     qBBBBBv    
        jMBMBNSVviiivirXMNFMNLYForUMkLvrrriviiiiiii irirLvuPMBE  iGMBMBBBBi    
       BBBBNBBVNMNoVXoXBBViMBBXoLrLujuYvrLrrjvivriiiiiiiiiirvVuGBBBBBOMBBu     
       iBBBVBMBqXYMMZkoNBVYMBNNNULivUvLviiiirLiiiiririiLrivvrrvFPVSOBBBBB      
     vvi NBqIBMkirIENMGGPqNMMNVPVLirvuvVUi vJuvvYYLrvvriLLJIrYoIUuVZBBMBMB     
    BBBBNvLIGOEFqNNIjqMNZMBBMGPPIoVuPSGjirrNuuUoNILvrvPuFSNXPOYiFOMOBBBBBBB    
     oMBBBoYVBBBMGkXGBENOBBBOVFOLkkPUSVPPFkGGGGPkkYSrjZGZPkSMBB   uBBBGBMBB    
        vGBBMZMMMZOMONNBBMGPoPEjojJVkvrUEPuYZVuUUSuEMEGNNoMBBBBi    NBBMBBB    
           rqBMBBBNPXZXkqNVYIOXPZPVVUrYoYurIuJvJPBEEqkoPOBBBBI        uBBBM    
              iYGBBBBMGkqkNGGPGMBBBBBBBMIuuZXvuXqNqSkNBBBBBNi           ii     
                  vBBBBBBOVZVuuVVNqNNNNGPkJFSSVZNquIEBBBBG                     
                     vSOBBBBMqXVkouLJjISqNNZMBBEqujBBBBMi                      
                         ivoSZOMMBMMNEZMBBMBBNPBBBBBMBB                        
                                   ii             rrri                         
                                                                                 
      iirrii     ir    iirrii        iii      ii    ir    ii    iri      ii    
     GBBBBBMB   BBBi  MBBBBBBMB    iBMBM     BBBB  qBBJ  BBBi  EBBq    kBBBM   
     BBO  BBB  YBBE   BBM  NBBY   uBBFBB    vBBBMv BMB  iBBBrYLBBB    OBEOBB   
    BMBBBBBG   BBB   BBBBBBBo    qBB iBB    BBOOBBEBB   BBBBBBBBB    MBN qBB   
   vMBN i     PBMu  vBBk BBB    BBBBMBBBi  PBB  BBBBS  jBBV  iBBG  iBBMBBBMB   
   BBB       iBBB   BBB  OBMB iBBBrrvBBBv iBBB  BBBB   MBB   BBB  PBBNrvLBBB   
   uri       ivv    Yr    rrr iUr     rvi  Yr    ri    Yr    jri  uLi    irv   
                                                                                                  
  
";
$usage="
piRanha: Transcriptome-wide prediction of piRNA target sites.
Version:  $version
Date:     $date

USAGE:
perl piRanha.pl -i piRNAs.fasta -r transcripts.fasta [-option (value)]


OPTIONS (default values):
-i OR -input           Name of the file that contains piRNA sequences in
                       FASTA format.
-r OR -reference       Name of the file that contains transcripts in FASTA
                       format.
-o OR -output          Name of the output folder (./piRanha_results)
-h OR -help            Prints this information and quits.
-keep_long_headers     Do not trim FASTA headers in transcript file after
                       first white space (off).
-skip_N                Skip sequences in piRNA input file that contain 'N'
                       bases (off).
-min_length            Minimum length [nt] of a piRNA to search for target
                       sites (24). Useful to skip e.g. miRNAs etc.
-max_length            Maximum length [nt] of a piRNA to search for target
                       sites (32). Useful to skip e.g. long tRFs etc.
-merge_results         Merge alignments for all transcripts into one text file
                       (off).
-SAM                   Generate alignment output in SAM format (off).
-use_scores            FASTA headers in piRNA input file refer to score
                       value, e.g. rpm count (off).
-max_score_decimals    Maximum number of decimals for a score in piRNA input
                       FASTA header (4). Longer numbers will be rounded.
-score_table           Create a file that contains score values for each
                       transcript (off).
-seed_size             Seed size (7). Seed starts with position 2.
-max_total_mm          Maximum number of total mismatch allowed for the entire
                       alignment (6). Note option -ignore_mm_beyond_pos.
-max_GU_mm_seed        Maximum number of GU mismatch allowed in seed (2).
-max_nonGU_mm_seed     Maximum number of non-GU mismatch allowed in seed (0).
-max_GU_mm_tail        Maximum number of GU mismatch allowed beyond seed (3).
-max_nonGU_mm_tail     Maximum number of non-GU mismatch allowed beyond seed
                       (2).
-ignore_mm_beyond_pos  Ignore any mismatch beyond this position (21). This is
                       to avoid bias towards alignments of short piRNAs.
-ignore_free_energy    Do not calculate free energy values (off).
-max_free_energy       Maximum free energy for an alignment (-25) [kJ/mol].
-silent                Do not show progress (off).


CONTACT:
David Rosenkranz
Institute of Organismic and Molecular Evolution
Johannes Gutenberg University Mainz, Germany
Email: rosenkranz\@uni-mainz.de
Web: www.smallRNAgroup.uni-mainz.de

";

### SET DEFAULTS ###

# input/output settings
$input_piR="piRNAs.fasta";
$input_transcript="transcripts.fasta";
$min_length=24;
$max_length=32;
$max_score_decimals=4;

# settings for piRNA targeting according to Wu et al. 2018 (option for nonGU mismatch in seed is not available yet)
$seed_size=7;
$max_overall_mm=6;
$max_GU_mm_seed=2;
$max_nonGU_mm_seed=0;
$max_GU_mm_3p=3;
$max_nonGU_mm_3p=2;
$ignore_mm_beyond_pos=21;
$max_free_energy=-25;


# collect command line arguments
GetOptions
	(
	# mandatory
	"i|input=s"=>\$input_piR,
	"r|reference=s"=>\$input_transcript,
	
	# change input/output settings
	"o|output=s"=>\$output_folder, # define an output folder
	"keep_long_headers"=>\$keep_long_headers, # do not trim FASTA headers in transcript file after first white space
	"skip_N"=>\$skip_N, # skip sequences in piRNA input file that contain N bases
	"min_length=i"=>\$min_length, # minimum length [nt] of a piRNA to search for target sites (skip e.g. miRNAs)
	"max_length=i"=>\$max_length, # maximum length [nt] of a piRNA to search for target sites (skip e.g. long tRFs)
	"merge_results"=>\$merge_results, # merge results for all transcripts into one text file
	"SAM"=>\$make_SAM, #Generate alignment output in SAM format
	"use_scores"=>\$use_scores, # FASTA headers in piRNA input FASTA file refer to score  (e.g. rpm count)
	"max_score_decimals=i"=>\$max_score_decimals, # maximum number of decimals for a score in piRNA input FASTA header (e.g. rpm count)
	"score_table"=>\$tr_score_table, # create a table file that contains score values for each transcript
	
	# change  piRNA targeting parameters
	"seed_size=i"=>\$seed_size, # change seed size of the targeting piRNA
	"max_total_mm=i"=>\$max_overall_mm, # maximum number of mismatches
	"max_GU_mm_seed=i"=>\$max_GU_mm_seed, # maximum number of GU pairs in seed
	"max_nonGU_mm_seed=i"=>\$max_nonGU_mm_seed, # maximum number of non-GU mm pairs in seed
	"max_GU_mm_tail=i"=>\$max_nonGU_mm_3p, # maximum number of GU pairs in tail
	"max_nonGU_mm_tail=i"=>\$max_nonGU_mm_3p, # maximum number of non-GU mm pairs in tail
	"ignore_mm_beyond_pos=i"=>\$ignore_mm_beyond_pos, # ignore mismatches beyond this position (be fair with long piRNAs)
	"ignore_free_energy"=>\$ignore_free_energy, # ignore free energy of piRNA-target alignments (marginal speed improvement)
	"max_free_energy=f"=>\$max_free_energy, # maximum allowed free energy for an alignment [kJ/mol]

	# other settings
	"t|threads=i"=>\$threads, # number of parallel threads
	"silent"=>\$silent, # do not print progress information and time estimation
	"h|help"=>\$print_help, # print help desk and quit
	)||die print"\nError in command line arguments.\n\n$usage\n\n";

if($print_help)
	{
	print$usage;
	exit;
	}

# define energy rules (Turner, 37\B0C). Last pair = 1st and 2nd base (target-piRNA), current pair = 3rd and 4th base (target-piRNA). Note: 1st and 3rd base are complementary
%TERrules=
	(
	"TTTT"=>-0.9,"TTGG"=>-1.8,"TTCC"=>-2.3,"TTAA"=>-1.1,"TTCT"=>-0.5,"TTAG"=>-0.7,
	"GGTT"=>-2.1,"GGGG"=>-2.9,"GGCC"=>-3.4,"GGAA"=>-2.3,"GGCT"=>-1.5,"GGAG"=>-1.5,
	"CCTT"=>-1.7,"CCGG"=>-2.0,"CCCC"=>-2.9,"CCAA"=>-1.8,"CCCT"=>-1.3,"CCAG"=>-1.9,
	"AATT"=>-0.9,"AAGG"=>-1.7,"AACC"=>-2.1,"AAAA"=>-0.9,"AACT"=>-0.7,"AAAG"=>-0.5,
	"CTTT"=>-0.9,"CTGG"=>-1.7,"CTCC"=>-2.1,"CTAA"=>-0.9,"CTCT"=>-0.5,"CTAG"=>-0.5,
	"AGTT"=>-0.9,"AGGG"=>-1.7,"AGCC"=>-2.1,"AGAA"=>-0.9,"AGCT"=>-0.6,"AGAG"=>-0.5,
	"0"=>0.0,"1"=>0.0,"2"=>0.8,"3"=>1.3,"4"=>1.7,"5"=>2.1,"6"=>2.5,"7"=>2.6,"8"=>2.8,
	);

# define match/mismatch characters for textual alignments
%alignment_char=
	(
	'AA'=>'|','TT'=>'|','GG'=>'|','CC'=>'|','GA'=>':','TC'=>':',
	'AT'=>' ','AC'=>' ','AG'=>' ','TA'=>' ','TG'=>' ','GT'=>' ','GC'=>' ','CA'=>' ','CT'=>' ','CG'=>' ',
	);

# check if input/output files can be opened for reading/writing
open(PIR,$input_piR)||die print"\nERROR! Unable to open input file $input_piR.\n$!\n\n";
open(TRANSCRIPT,$input_transcript)||die print"\nERROR! Unable to open input file $input_transcript.\n$!\n\n";
if(!$output_folder)
	{
	($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
	$year+=1900;
	$mon++;
	$mday=substr("0".$mday,-2);
	$mon=substr("0".$mon,-2);
	$fID=0;
	while(1)
		{
		$fID++;
		$output_folder="piRanha.$year-$mon-$mday.#$fID";
		if(!-d$output_folder)
			{
			mkdir($output_folder)||die print"\nERROR! Unable to create output folder $output_folder.\n$!\n\n";
			last;
			}
		}
	}
if($tr_score_table)
	{
	open(TRSCORE,">$output_folder/transcript.scores")||die print"\nERROR! Unable to create output file $output_folder/transcript.scores.\n$!\n\n";
	}
print$logo;
print"
piRNA sequence file: ... $input_piR
reference file: ........ $input_transcript

min. input sequences length: ................ $min_length
max. input sequences length: ................ $max_length
max. decimals for scores: ................... $max_score_decimals
seed size: .................................. $seed_size
max. overall mismatch: ...................... $max_overall_mm
max. GU-whobbles in seed: ................... $max_GU_mm_seed
max. non-GU mismatch in seed:................ $max_nonGU_mm_seed
max. GU-whobbles beyond seed: ............... $max_GU_mm_3p
max. non-GU mismatch beyond seed: ........... $max_nonGU_mm_3p
ignore any mismatch beyond position: ........ $ignore_mm_beyond_pos
max. free energy for piR-target alignment: .. $max_free_energy
";

# load and index piRNA input data 
if($max_nonGU_mm_seed==0) # only GU mismatches = fast
	{
	print"\nLoad and index sequence data from $input_piR (.=10k sequences) ";
	$i=0;
	%piR=();
	%score_per_seq=();
	while($input=<PIR>)
		{
		$input=~s/\s*$//;
		if($input=~s/^>//)
			{
			$i++;
			if($input=~/^\s*\d+\.*\d*\s*$/) # reduce number of decimals
				{
				$input=(int(($input*(10**$max_score_decimals))+0.5))/(10**$max_score_decimals);
				}
			$score=$input;
			if($i=~/0000$/)
				{
				print".";
				}
			}
		else
			{
			next if(length$input<$seed_size);
			next if($skip_N&&$input=~/N/);
			next if(length$input<$min_length||length$input>$max_length);
			$input=~s/U/T/g;
			$seed=substr($input,1,$seed_size);
			alt_seeds();
			foreach$var(@vars)
				{
				$GU=0;
				foreach$p(0..$seed_size-1)
					{
					if(substr($var,$p,1) ne substr($seed,$p,1))
						{
						$GU++;
						last if($GU>$max_GU_mm_seed);
						}
					}
				if($GU<=$max_GU_mm_seed)
					{
					# structure of the array: 3' tail sequence, original seed, 5' base, GU mismatch in seed, non-GU mismatch in seed
					push(@{$piR{$var}},substr($input,$seed_size+1)."#".$seed."#".substr($input,0,1)."#".$GU."#"."0");
					}
				}
			$score_per_seq{$input}=$score;
			}
		}
	close PIR;
	$n=keys%piR;
	print" done.\nCreated $n seed-indexes for $i sequences.";
	}
else # allow non-GU mismatches = slow
	{
	%seeds2target=();
	@vars=();
	@data=();
	foreach(1..$seed_size)
		{
		push(@data,[qw / A T G C /]);
		}
	$iterator=Iterator::Array::Jagged->new(data=>\@data);
	while(@set=$iterator->next)
		{
		$perm=join('',@set);
		push(@vars,$perm);
		}
	
	print"
Generate seed-to-target matrix for seeds containing non-GU mismatches.

0%                                   50%                                   100%
|-------------------------------------|-------------------------------------|
";
	$i=0;
	$dots=0;
	foreach$seedX(@vars)
		{
		$i++;
		if($i/@vars>=1/77)
			{
			$i=0;
			$dots++;
			print".";
			}
		foreach$targetX(@vars)
			{
			$ok=1;
			$GU=0;
			$nonGU=0;
			foreach$p(0..$seed_size-1)
				{
				$pair=substr($seedX,$p,1).substr($targetX,$p,1);
				if(substr($seedX,$p,1) ne substr($targetX,$p,1))
					{
					if($pair eq'GA'||$pair eq'TC')
						{
						$GU++;
						}
					else
						{
						$nonGU++;
						}
					}
				if($GU>$max_GU_mm_seed||$nonGU>$max_nonGU_mm_seed||$GU+$nonGU>$max_overall_mm)
					{
					$ok=0;
					last;
					}
				}
			if($ok==1)
				{
				$seeds2target{$seedX}{$targetX}=[$GU,$nonGU];
				}
			else
				{
				$bad++;
				}
			}
		}
	foreach($dots..76)
		{
		print".";
		}
	
	print"\n\nLoad and index sequence data from $input_piR (.=10k sequences) ";
	$i=0;
	%piR=();
	%score_per_seq=();
	while($input=<PIR>)
		{
		$input=~s/\s*$//;
		if($input=~s/^>//)
			{
			$i++;
			if($input=~/^\s*\d+\.*\d*\s*$/) # reduce number of decimals
				{
				$input=(int(($input*(10**$max_score_decimals))+0.5))/(10**$max_score_decimals);
				}
			$score=$input;
			if($i=~/0000$/)
				{
				print".";
				}
			}
		else
			{
			next if(length$input<$seed_size);
			next if($skip_N&&$input=~/N/);
			next if(length$input<$min_length||length$input>$max_length);
			$input=~s/U/T/g;
			$seed=substr($input,1,$seed_size);
			foreach$var(keys%{$seeds2target{$seed}})
				{
				# structure of the array: 3' tail sequence, original seed, 5' base, GU mismatch in seed, non-GU mismatch in seed
				push(@{$piR{$var}},substr($input,$seed_size+1)."#".$seed."#".substr($input,0,1)."#".$seeds2target{$seed}{$var}[0]."#".$seeds2target{$seed}{$var}[1]);
				}
			$score_per_seq{$input}=$score;
			}
		}
	undef%seeds2target;
	close PIR;
	$n=keys%piR;
	print" done.\nCreated $n seed-indexes for $i sequences.";
	}

# read transcript input data
print"\n\nStart searching target sites in $input_transcript.\nProcessed bp: ";
$num_id=0;
$transcript="";
$processed_bp=0;
while(<TRANSCRIPT>)
	{
	$_=~s/\s*$//;
	if($_=~s/^>//)
		{
		if($transcript ne "")
			{
			$num_id++;
			screen();
			}
		
		# show progress
		if(!$silent)
			{
			print"\b" x $backspace;
			$procesed_bp+=length$transcript;
			print$procesed_bp;
			$backspace=length$procesed_bp;
			}
		
		$head=$_;
		$transcript="";
		if(!$keep_long_headers)
			{
			$head=~s/ .+$//
			}
		}
	else
		{
		$_=~s/U/T/g;
		$transcript.=uc$_;
		$total_ref_length+=length$_;
		}
	}
screen();
close TRSCORE;
close TRANSCRIPT;

sub screen
	{
	open(OUTPUT,">$output_folder/$num_id.txt")||print"\nERROR! Unable to create output file $output_folder/$num_id.txt.\n$!\nContinue nevertheless...\n\n";
	$valid_alignments_found=0;
	$transcript=reverse$transcript;
	$transcript=~tr/ATGC/TACG/;
	$p=-1;
	$score=0;
	$alignments_in_this_transcript=0;
	while($p<length$transcript)
		{
		$p++;
		if($piR{substr($transcript,$p,$seed_size)})
			{
			foreach$data(@{$piR{substr($transcript,$p,$seed_size)}})
				{
				@data=split('#',$data);
				$overall_mm=$data[3]+$data[4];
				$GU_mm_seed=$data[3];
				$nonGU_mm_seed=$data[4];
				$GU_mm_3p=0;
				$nonGU_mm_3p=0;
				$valid_alignment=1;
				foreach$p3(0..(length$data[0])-1)
					{
					if(substr($data[0],$p3,1) ne substr($transcript,$p+$seed_size+$p3,1))
						{
						$overall_mm++;
						if(substr($data[0],$p3,1)eq'G'&&substr($transcript,$p+$seed_size+$p3,1)eq'A')
							{
							$GU_mm_3p++;
							}
						elsif(substr($data[0],$p3,1)eq'A'&&substr($transcript,$p+$seed_size+$p3,1)eq'G')
							{
							$GU_mm_3p++;
							}
						else
							{
							$nonGU_mm_3p++;
							}
						}
					if($overall_mm>$max_overall_mm||$GU_mm_3p>$max_GU_mm_3p||$nonGU_mm_3p>$max_nonGU_mm_3p)
						{
						$valid_alignment=0;
						last;
						}
					$export_p3=$p3;
					last if($seed_size+$p3+2==$ignore_mm_beyond_pos);
					}
				if($valid_alignment==1)
					{
					$output="";
					$GU_X=0;
					$nonGU_X=0;
					$overall_mm_plusX=$overall_mm;
					$alignments_in_this_transcript++;
					$piR_seq=$data[2].$data[1].$data[0];
					$target_seq=substr($transcript,$p-1,length$piR_seq);
					$alignment="";
					$loop_size=0;
					$free_energy=0;
					$pairs_in_structure=0;
					foreach$p_int(0..(length$piR_seq)-1)
						{
						$alignment.=$alignment_char{substr($piR_seq,$p_int,1).substr($target_seq,$p_int,1)};
						
						# mismatch counts of mismatches occurring in neglected 3' part (default: beyond position 24)
						if($p_int>=$ignore_mm_beyond_pos) 
							{
							if($alignment=~/:$/)
								{
								$GU_X++;
								$overall_mm_plusX++;
								}
							elsif($alignment=~/ $/)
								{
								$nonGU_X++;
								$overall_mm_plusX++;
								}
							}
						
						if(!$ignore_free_energy)
							{
							if($alignment_char{substr($piR_seq,$p_int,1).substr($target_seq,$p_int,1)}ne' ')
								{
								$pairs_in_structure++;
								if($pairs_in_structure>1)
									{
									$free_energy+=$TERrules{$prev_pair.substr($target_seq,$p_int,1).substr($piR_seq,$p_int,1)};
									$free_energy+=$TERrules{$loop_size};
									}
								$loop_size=0;
								$prev_pair=substr($target_seq,$p_int,1).substr($piR_seq,$p_int,1);
								}
							else
								{
								$loop_size++;
								}
							}
						}
					$target_seq=~tr/ATGC/TACG/;
					$start_in_target=(length$transcript)-$p+2-length$piR_seq;
					$end_in_target=($start_in_target+length$piR_seq)-1;
					
					$output.="#$alignments_in_this_transcript $head vs. $score_per_seq{$piR_seq} ($start_in_target-$end_in_target). GU-S:$GU_mm_seed nonGU-S:$nonGU_mm_seed GU-T:$GU_mm_3p nonGU-T:$nonGU_mm_3p GU-X:$GU_X nonGU-X:$nonGU_X TOTAL:$overall_mm TOTAL+X:$overall_mm_plusX.";
					if(!$ignore_free_energy)
						{
						$output.=" Free energy:$free_energy kJ/mol.";
						}
					$output.="\nQ: 5'-$piR_seq-3'\n      $alignment\nR: 3'-$target_seq-5'\n\n";
					if($ignore_free_energy)
						{
						$valid_alignments_found++;
						print OUTPUT$output;
						}
					elsif(!$ignore_free_energy&&$free_energy<=$max_free_energy)
						{
						$valid_alignments_found++;
						print OUTPUT$output;
						}
					
					# save read counts per transcript (works only if FSATA header of sRNA data refers to read counts /diff-expr. score)
					if($use_scores)
						{
						$score+=$score_per_seq{$piR_seq};
						}
					}
				}
			}
		}
	if($use_scores)
		{
		print OUTPUT"\nTotal score for $head: $score";
		if($tr_score_table){print TRSCORE"$head\t$score\n";}
		}
	close OUTPUT;
	if($valid_alignments_found==0) # remove empty output files
		{
		unlink"$output_folder/$num_id.txt";
		}
	}

if($make_SAM)
	{
	print"\nCreate SAM output file...";
	open(SAM,">$output_folder/piRanha.sam")||die print"\nERROR! Unable to create SAM output file $output_folder/piRanha.sam.\n$!\n\n";
	print SAM"\@HD\tVN:1.5\n";
	foreach$id(sort{$a cmp $b}keys%transcripts)
		{
		$LN=length$transcripts{$id};
		print SAM"\@SQ\tSN:$id\tLN:$LN\n";
		}
	$num_id=0;
	while(1)
		{
		$num_id++;
		open(SPLIT,"$output_folder/$num_id.txt")||last;
		while(<SPLIT>)
			{
			chomp$_;
			if($_=~s/^#\d+ //)
				{
				$_=~s/^[^ ]+//;
				$ref=$&;
				$_=~s/ \(\d+-\d+\)//;
				$coordinates=$&;
				$coordinates=~s/[ \(\)]//g;
				@coordinates=split('-',$coordinates);
				$_=~s/^.* vs. //;
				$_=~s/.*\. GU-S://;
				$query=$&;
				$query=~s/\. GU-S://;
				$_=~s/^.*Free energy://;
				$_=~s/ kJ\/mol.$//;
				$FE=$_;
				}
			elsif($_=~s/^Q: 5'-//)
				{
				$_=~s/-3'$//;
				$seq=reverse$_;
				$seq=~tr/ATGC/TACG/;
				}
			elsif($_=~s/^      //)
				{
				$CIGAR="";
				while(1)
					{
					if($_=~s/[ :]+$//)
						{
						$i=length$&;
						$CIGAR.=$i."=";
						}
					elsif($_=~s/\|+$//)
						{
						$i=length$&;
						$CIGAR.=$i."X";
						}
					last if(length$_==0);
					}
				}
			elsif($_=~s/^R: 3'-//)
				{
				$_=~s/-5'$//;
				$Tlen=length$seq;
				print SAM"$query\t16\t$ref\t$coordinates[0]\t255\t$CIGAR\t*\t0\t$Tlen\t$seq\t*\tFE:f:$FE\n";
				}
			}
		close SPLIT;
		}
	close SAM;
	print" done.";
	}

if($merge_results)
	{
	print"\nMerge output files...";
	open(MERGE,">$output_folder/piRanha.results")||die print"\nERROR! Unable to create merge output file $output_folder/piRanha.results.\nNote: Target prediction has finished without errors but results were not merged into one single file.\n\n";
	opendir(DIR,$output_folder)||die print"\nERROR! Unable to read output folder $output_folder.\nNote: Target prediction has finished without errors but results were not merged into one single file.\n\n";
	while($file=readdir DIR)
		{
		if($file=~/\.txt$/)
			{
			open(SPLIT,"$output_folder/$file")||print"\nERROR! Unable to read output file $output_folder/$_.\n$!\nContinue nevertheless.\n\n";
			while(<SPLIT>)
				{
				print MERGE$_;
				}
			close SPLIT;
			unlink"$output_folder/$file" or warn "\n$file: $!";
			}
		}
	closedir DIR;
	close MERGE;
	print" done.";
	}

sub alt_seeds
	{
	@vars=();
	@data=();
	@seed=split('',$seed);
	foreach$nt(@seed)
		{
		if($nt eq "A"){push(@data,[qw / A /]);}
		elsif($nt eq "T"){push(@data,[qw / T C /]);} # allow UG whoble
		elsif($nt eq "G"){push(@data,[qw / G A /]);} # allow GU whoble
		elsif($nt eq "C"){push(@data,[qw / C /]);}
		elsif($nt eq "N"){push(@data,[qw / A T G C /]);}
		}
	$iterator=Iterator::Array::Jagged->new(data=>\@data);
	while(@set=$iterator->next)
		{
		$perm=join('',@set);
		push(@vars,$perm);
		}
	}
