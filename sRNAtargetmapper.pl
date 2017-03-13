#!/usr/bin/perl
use Getopt::Long qw(:config no_ignore_case);
use Data::Dumper;
$|=1;
$script=$0;
$script=~s/^.+[\\\/]//;

$information="
 ============================= sRNAtargetmapper ===============================
 VERSION 1.0                                        LAST MODIFIED: 21. Dec 2016

 Please cite:
 Roovers EF, Rosenkranz D, Mahdipour M, Chung-Ting H, He N, Chuva de Sousa
 Lopes SM, van der Westerlaken LAJ, Zischler H, Butter F, Roelen BAJ, Ketting
 RF. 2015. Piwi proteins and piRNAs in mammalian oocytes and early embryos.
 Cell Rep in press.

 SCOPE:
 sRNAtargetmapper is specifically designed to detect potential  small RNA targets 
 in genomes.  To this end it uses a specialized mapping algorithm that requires a 
 perfect 5' seed match (default: 10nt) and optionally allows non-template 3' ends
 as well as internal mismatches in the part of the sequence that follows
 the seed region. Allowing non-template 3' ends will ensure the mapping of
 3' modified (adenylated/uridylated) small RNAs while allowing internal mis-
 matches can enhance sensitivity considering degressive read quality towards
 3' ends.
 
 ------------------------------------------------------------------------------
 
 EXEMPLARY ALIGNMENT OF VALID HIT USING DEFAULT PARAMETERS:
 
 genome  5'-TGCGTACTAGCTCGCATGACTCGCATAGCTACGTGGTAGC-3'
                  x|||||||||||||||||x||||||xx
 sRNA          5'-TTAGCTCGCATGACTCGCTTAGCTAAA-3'
                  C                        NN 
                  |----S----||-----T------|
 
 C - clipped nucleotides from the seed. 
 S - seed region
 T - template region
 N - non-template nucleotides
 ------------------------------------------------------------------------------
 
 USAGE:
 perl $script -i input.fasta -g genome.fasta (-option [value])
 
 *You may also use FASTAQ-formatted files.

 ------------------------------------------------------------------------------

 OPTIONS:
 -f [value]
 -format [value]
             Specifies the output format. Valid values are <sam>, <eland> and
             <compact> <eeland>. The default value is <eeland>. Extended eland
             displays mismatch values for clipped mismatches, seed mismatches,
             template mismatches and non-template mismatches.

 -g [value]
 -genome [value]
             Name of the reference/genome file (FASTA). *REQUIRED
 -i [value]
 -input [value]
             Name of the input file (FASTA). Multiple input files are allowed.
             Input files will be processed in parallel. *REQUIRED
 -o [value]
 -output [value]
             Name of the output file. By default the output file will be named
             input.fasta.map (input file: input.fasta).
 -m [value]
 -mismatch [value]
             Defines the maximum internal mismatch in the part of the sequence
             that follows the seed region. The default value is <1>.
 -n [value]
 -nontemplates [value] 
             Defines the number of allowed non-template 3' nucleotides. The
             default value is <2> which means that the last 2 nucleotides do
             not need to match the reference to produce a valid hit.
 -a [value]
 -alignments [value]
             Strictness of the output. Valid values are <all> and <best>. <all>
             will output all hits while <best> will output only the best hits.
             Can be still more than one hit per sequence if the hits have equal
             quality. The default value is <all>.
 -p
 -printalignment
             Print the alignment (eeland format only)
 -s [value]
 -seedmatch [value]
             Defines the number of 5' nucleotides that must produce a perfect
             hit (seed match). The default value is <10>.
 
 -S [value]
 -seedmismatch [value]
             Defines the number of nucleotides within the seedmatch which may
             differ from the genomic sequence. Note increasing this value
             has an exponential run time increase.
 -c 
 -clipfirst
             Defines the 5' clip on the seed match. This is intended to allow a
             1-base 5' mismatch as is commonly seen in target complementarity.
 -r [value]
 -replaceheads [value]
             Titles of the reference sequence (FASTA headers) can optionally be
             replaced by a consecutively numbering to reduce output size and
             save diskspace. Valid values are <0> (keep original names) and <1>
             (replace headers by numbers). The default value is <0>.
 -d [value]
 -direction [value]
             Defines if mapping is performed on plus strand only or on both
             strands. Valid values are <1> for plus strand only, <-1> for
             antisense strand (best for target searching of a transcriptome)
             or <0> to search both strands. The default value is <-1>.
 -G [value]
 -GUpenalty
             Penalty for GU pairs in the template region. Default is 0.5.
 -v [value]
 -verbosity [value]
             Defines verbosity of the program during the mapping process. Valid
             values are <0> (silent), <1> (concise) and <2> (verbose). The
             default value is <2>.
 -h
 -help
             Will print this information.

 ------------------------------------------------------------------------------

 (c) David Rosenkranz
 Institute of Anthropology, small RNA group
 Johannes Gutenberg University Mainz, Germany
 Contact: rosenkranz\@uni-mainz.de
 Visit our homepage at:
 www.smallRNAgroup-mainz.de
 or
 www.smallRNAgroup.uni-mainz.de
 
 ##############################################################################


";

$outfile="";
GetOptions
	(
	"h"=>\$print_info,
	"help"=>\$print_info,
	"i=s"=>\@input,
	"input=s"=>\@input,
	"o=s"=>\$outfile,
	"output=s"=>\$outfile,
	"g=s"=>\$genome_file,
	"genome=s"=>\$genome_file,
	"s=s"=>\$perfect_5p,
	"seedmatch=s"=>\$perfect_5p,
	"S=s"=>\$seed_mm,
	"seedmismatch=s"=>\$seed_mm,
	"c"=>\$clip,
	"clipfirst"=>\$clip,
	"G=f"=>\$GUpenalty,
	"GUpenalty=f"=>\$GUpenalty,
	"p"=>\$printalignment,
	"printalignment"=>\$printalignment,
	"n=s"=>\$untemplate_3p,
	"nontemplates=s"=>\$untemplate_3p,
	"m=s"=>\$max_internal_mm,
	"mismatch=s"=>\$max_internal_mm,
	"d=s"=>\$direction,
	"direction=s"=>\$direction,
	"r=s"=>\$replace_titles_by_id,
	"replaceheads=s"=>\$replace_titles_by_id,
	"f=s"=>\$format,
	"format=s"=>\$format,
	"a=s"=>\$output_strictness,
	"alignments=s"=>\$output_strictness,
	"v=s"=>\$verbosity,
	"verbosity=s"=>\$verbosity,
	);

$format=lc$format;
if($print_info){print$information;exit;}
if($clip){$clip=1;}else{$clip=0;}
if($printalignment){$printalignment=1;}else{$printalignment=0;}
if($seed_mm!~/^\d+$/){$seed_mm=0;}
if($perfect_5p!~/^\d+$/){$perfect_5p=10;}
if($untemplate_3p!~/^\d+$/){$untemplate_3p=2;}
if($max_internal_mm!~/^\d+$/){$max_internal_mm=1;}
if($direction!~/^-?[01]$/){$direction=-1;}
if($replace_titles_by_id!~/^[01]$/){$replace_titles_by_id=0;}
if($format ne'eland'&&$format ne'compact'&&$format ne'sam'){$format="eeland";}
if($output_strictness ne 'all'&&$output_strictness ne 'best'){$output_strictness="all";}
if($GUpenalty!~/^\d+(\.\d+)?$/){$GUpenalty=0.5;}
if($verbosity!~/^[012]$/){$verbosity=2;}
if(@input==0){print"\n\n STOP: No probe/input file(s) specified.\n\n";exit;}
unless(-e$genome_file){print"\n\n STOP: No valid reference/genome file specified.\n\n";exit;}
if($direction == 1) {$strands = "+";} elsif($direction == -1) {$strands = "-";} else {$strands = "+/-";}

print"

             ------------------------------------------------
             sRNAtargetmapper Version 1.0        21. Dec 2016
             ------------------------------------------------
             Starting sRNAmapper with the following settings:


             Seed size: ................ $perfect_5p
             5' clip: .................. $clip
             Max. non-template 3' bases: $untemplate_3p
             Max. internal mismatch: ... $max_internal_mm
             Max. seed mismatch: ....... $seed_mm
             G:U penalty: .............. $GUpenalty
             Search strands: ........... $strands
             Output format: ............ $format 
             Alignments: ............... $output_strictness
     
             (c) David Rosenkranz
             Institute of Anthropology, small RNA group
             Johannes Gutenberg University Mainz, Germany
             Contact: rosenkranz\@uni-mainz.de
             ------------------------------------------------


";

sub index_combinations 
	{
	my ($n,$r) = @_;
	return if $r > $n;
	my @indices = (0..($r-1));
	my $ret = [];
	push(@$ret,[@indices]);
	while(1) 
		{
		my $index = -1;
		for my $i (reverse(0..($r-1))) 
			{
			if($indices[$i] != $i+$n-$r) 
				{
				$index = $i; last; 
				}
			}
		return @$ret if($index == -1);
		$indices[$index] ++;
		for my $i ($index+1..($r-1)) 
			{
			$indices[$i] = $indices[$i-1] + 1;
			}
		push(@$ret,[@indices]);
		}
	return @$ret;
	}

sub cartesian_product 
	{
	my $sets = shift @_;

	my $products = [[]];
	for my $set (reverse @$sets) 
		{
		my $partial_products = $products;
		$products = [];
		for my $item (@$set) 
			{
			for my $product (@$partial_products) 
				{
				push @$products, [ $item, @$product ];
				}
			}
		}

	return $products;
	}

# Generates all sequences of a at most the given edit distance from the 
# original sequence
#
# Returns hash containing all the sequences as keys mapping to the edit 
# distance of each sequence from the original sequence.

sub generate_edited_seqs 
	{
	my ($orig_seq,$n_edit_sites) = @_;
	my @site_list;
	for my $i (1..$n_edit_sites) 
		{
		push(@site_list,index_combinations(length($orig_seq),$i));
		}
	my  %seqs = ($orig_seq => 0);
	for my $edit_site_set (@site_list) 
		{
		my @repl_list;
		for my $p (@$edit_site_set) 
			{
			my @repl = grep { $_ ne substr($orig_seq,$p,1) } ("A","C","G","U");
			push(@repl_list,\@repl);
			}
		my $repl_prod = cartesian_product(\@repl_list);
		for my $i (0..$#$repl_prod) 
			{
			my $cur_seq = $orig_seq;
			for my $j (0..$#$edit_site_set) 
				{
				my ($c,$p) = ($repl_prod->[$i]->[$j],$edit_site_set->[$j]);
				$cur_seq = substr($cur_seq,0,$p).$c.substr($cur_seq,$p+1);
				}
			$seqs{$cur_seq} = scalar(@$edit_site_set);
			}
		}
	return \%seqs;
	}

@childs=();
foreach$probe(@input)
	{
	$pid=fork();
	if($pid) # parent
		{
		push(@childs,$pid);
		}
	elsif($pid==0) # child
		{
		# CHECK INPUT FILE FORMAT
		if($verbosity!=0)
			{
			print"\nChecking input file format...";
			}
		open(IN,$probe)||die print"\nCould not open input file $probe.\n$!\n\n";
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
			$iformat='FASTA';
			}
		elsif($fastq_heads>$fasta_heads)
			{
			$iformat='FASTQ';
			}
		else
			{
			die print"\nUnknown file format of input file $probe.\n\n";
			}
		if($verbosity!=0)
			{
			print" done. -> $iformat";
			}
		
		# INDEX PROBES (SENSE)
		$min_length=-1;
		$max_length=0;
		$n_probeseqs=0;
		%probes=();
		%probes_info=();
		%probes_seeds=();
		%probes_distance=();
		$t_start=time;
		if($direction != -1) 
		{
		open(PROBE,$probe)||die print"\nCould not open sequencefile $probe.\n$!\n\n";
		if($verbosity!=0)
			{
			print"\nIndexing sequences from $probe...";
			}
		$i=0;
		@seq_data=();
		while(<PROBE>)
			{
			$_=~s/\s*$//; # chomp depends on the operating system and can cause problems when processing files from different platforms (Unix vs. Windows)
			push(@seq_data,$_);
			if($iformat eq'FASTQ')
				{
				$i+=0.25;
				}
			elsif($iformat eq'FASTA')
				{
				$i+=0.5;
				}

			if($i==1) 
				{
				$seq_data[1]=~tr/atgcuT/AUGCUU/;
				my $primary_seed = substr($seq_data[1],$clip,$perfect_5p-$clip);
				if($seq_data[1]!~/^\s*$/&&$primary_seed!~/N/)
					{
					$n_probeseqs++;
					
					# hash of arrays: key = 5p seq / array = all 3p seqs that belong to 5p seq
					my $seeds = generate_edited_seqs($primary_seed,$seed_mm);
					for my $seed (keys %$seeds) 
						{
						push(@{$probes{$seed}},$seq_data[1]);	
						$seq_data[0]=~s/^[>\@]//; $seq_data[0]=~s/\s.*$//;
						push(@{$probes_info{$seed}},$seq_data[0]);
						push(@{$probes_distance{$seed}},$$seeds{$seed});
						}
					
					# check minimum and maximum sequence length
					if(length$seq_data[1]<$min_length||$min_length==-1)
						{
						$min_length=length$seq_data[1];
						}
					if(length$seq_data[1]>$max_length)
						{
						$max_length=length$seq_data[1];
						}
					}
				$i=0;
				@seq_data=();
				}
			}
		close PROBE;
		$n_indexes=keys%probes;
		print" done.\nCreated $n_indexes indexes for $n_probeseqs sequences ($min_length-$max_length nt)";
		}
		
		# INDEX PROBES (FOR REVERSE STRAND SEARCH)
		if($direction!=1)
			{
			%probes_rc=();
			%probes_rc_info=();
			%probes_rc_seeds=();
			%probes_rc_distance=();

			open(PROBE,$probe)||die print"\nCould not open sequencefile $probe.\n$!\n\n";
			if($verbosity!=0)
				{
				print"\nIndexing sequences for reverse strand search...";
				}
			$i=0;
			while(<PROBE>)
				{
				$_=~s/\s*$//;
				push(@seq_data,$_);
				if($iformat eq'FASTQ')
					{
					$i+=0.25;
					}
				elsif($iformat eq'FASTA')
					{
					$i+=0.5;
					}
				if($i==1)
					{
					$seq_data[1]=~tr/atgcuT/AUGCUU/;
					if($seq_data[1]!~/^\s*$/&&substr($seq_data[1],$clip,$perfect_5p-$clip)!~/N/)
						{
						$rc=$seq_data[1];
						$rc=~tr/AUGC/UACG/;
						$rc=reverse$rc;
						my $primary_seed = substr($rc,-($perfect_5p),$perfect_5p-$clip);
						my $seeds = generate_edited_seqs($primary_seed,$seed_mm);
						for my $seed (keys %$seeds) 
							{
							push(@{$probes_rc{$seed}},$rc);
							$seq_data[0]=~s/^[>\@]//; $seq_data[0]=~s/\s.*$//;
							push(@{$probes_rc_info{$seed}},$seq_data[0]);
							push(@{$probes_rc_distance{$seed}},$$seeds{$seed});
							}
						# check minimum and maximum sequence length
						if(length$seq_data[1]<$min_length||$min_length==-1)
							{
							$min_length=length$seq_data[1];
							}
						if(length$seq_data[1]>$max_length)
							{
							$max_length=length$seq_data[1];
							}
						}
					$i=0;
					@seq_data=();
					}
				}
			close PROBE;
			print" done.";
			}
		
		# READ GENOME, SAVE INFORMATION AND OPTIONALLY REPLACE FASTA TITLES BY IDs
		if($verbosity!=0)
			{
			print"\nRead genome file...";
			}
		open(GENOME,$genome_file)||die print"\nCould not open genome file $genome_file.\n$!\n\n";
		open(RAW,">$probe.raw");
		if($replace_titles_by_id==1)
			{
			print RAW"# original fasta header\t->\treplaced by ID\n";
			}
		
		$replace_id=0;
		$total_genome_size=0;
		@sam=();
		%sam=();
		while(<GENOME>)
			{
			$_=~s/\s*$//;
			if($_=~/>/)
				{
				$_=~s/^>//;
				$title=$_;
				push(@sam,$title);
				if($replace_titles_by_id==1)
					{
					$replace_id++;
					$replace_titles{$_}=$replace_id;
					print RAW"# $_\t->\t$replace_id\n";
					}
				}
			else
				{
				$sam{$title}+=length$_;
				$total_genome_size+=length$_;
				}
			}
		close GENOME;
		
		# START MAPPING
		if($verbosity!=0)
			{
			print" done. Total size: $total_genome_size"."bp"."\nStart mapping...";	
			}
		$chromosome="";
		$position=-1;
		$processed_bp=0;
		$raw_lines=0;
		$next_line="#";
		$stretch="##";
		$lagging_stretch="";
		%best_hit=();
		$t0=time;
		$t_start=time;
		open(GENOME,$genome_file)||die print"\nCould not open genome file $genome_file.\n$!\n\n";
		while(length$next_line>0)
			{
			# progress
			if($verbosity==2)
				{
				$t1=time;
				if($t1-$t0>=10)
					{
					$t0=time;
					$speed=int(($processed_bp/(time-$t_start))+0.5);
					$est_remaining="?";
					if($speed>0)
						{
						$est_remaining=int((($total_genome_size-$processed_bp)/$speed)+0.5);
						}
					print"\nPROBE:.. $probe\nREF:.... $chromosome\nPOS:.... $position\nScanned: $processed_bp bp\nSpeed:.. $speed bp/s\nEstimated remaining time: $est_remaining seconds\n";
					}	
				}
			
			if($concatenate_index==0)
				{
				sliding_window_increment();
				sub sliding_window_increment
					{
					# sliding window +1bp
					$stretch=~s/.//;
					$lagging_bp=$&;
					$position++;
					
					# save passed sequence for reverse strand search (extension of seed match)
					#if($direction!=1)
					#{
						$lagging_stretch.=$lagging_bp;
						if(length$lagging_stretch>$max_length)
							{
							$lagging_stretch=~s/.//;
							}
						#}
					}
				}
			
			if(length$stretch>=$perfect_5p-$clip)
				{
				$concatenate_index=0;
				search_seed_matches();
				sub search_seed_matches
					{
					# check for seed match in sense orientation
					my $seed = substr($stretch,$clip,$perfect_5p-$clip);
					if($probes{$seed})
						{
						foreach$i(0..scalar(@{$probes{$seed}})-1)
							{
							my $mapped_seq = $probes{$seed}->[$i];
							my $p3 = substr($mapped_seq,$perfect_5p); 
							$count_internal_mm=0;
							$count_untemplate_3p=0;
							$count_clip_mm = 0;
							foreach$pos(0..$clip-1) 
								{
									if(substr($stretch,$pos,1) ne substr($mapped_seq,$pos,1)) 
									{
										$count_clip_mm++;
									}
								}
							foreach$pos(0..((length$p3)-$untemplate_3p)-1)
								{
								$g=substr($stretch,$pos+$perfect_5p,1);
								$m=substr($p3,$pos,1);
								if($g ne $m)
									{
									if(($g eq "C" && $m eq "U") || ($g eq "A" && $m eq "G"))
										{
										$count_internal_mm += $GUpenalty;
										}
									else
										{
										$count_internal_mm++;
										}
									last if($count_internal_mm>$max_internal_mm);
									}
								}
							if($count_internal_mm<=$max_internal_mm)
								{
								$locus_seq=substr($stretch,0,$perfect_5p+length$p3);
								foreach$pos(0..$untemplate_3p-1)
									{
									if(substr($p3,-$untemplate_3p+$pos,1)=~/[ATGCN]/&&substr($p3,-$untemplate_3p+$pos,1) ne substr($locus_seq,-$untemplate_3p+$pos,1))
										{
										$count_untemplate_3p++;
										}
									}
								$hit_mm=$count_internal_mm+$count_untemplate_3p+$count_clip_mm+$probes_distance{$seed}->[$i];

								my $probe_name = $probes_info{substr($stretch,$clip,$perfect_5p-$clip)}->[$i];

								$output_line = "$replace_titles{$chromosome}\t$position\t$locus_seq\t$probe_name\t$mapped_seq\t$hit_mm\t+\n";
								if($format eq 'eeland') 
									{
									$output_line = substr($output_line,0,-1)."\t$count_clip_mm\t$probes_distance{$seed}->[$i]\t$count_internal_mm\t$count_untemplate_3p\n";
									$m = "";
									if($printalignment) 
										{
										foreach$pos(0..length($locus_seq)-1) 
											{
												$l = substr($mapped_seq,$pos,1); $p = substr($locus_seq,$pos,1);
												if($l eq $p) { $m .= "|" }
												elsif(($l eq "G" && $p eq "A") || ($l eq "U" && $p eq "C")) { $m .= ":" }
												else { $m .= " " }
											}
											$output_line .= "${locus_seq}\n${m}\n${mapped_seq}\n";
										}
									}

								
								# do not output hit if output_strictness is 'best' and it is worse than a hit before
								unless(exists($best_hit{$mapped_seq}))
									{
									$best_hit{$mapped_seq}=$hit_mm;
									print RAW $output_line;
									$raw_lines++;
									}
								else
									{
									if($hit_mm<=$best_hit{$mapped_seq})
										{
										$best_hit{$mapped_seq}=$hit_mm;
										print RAW $output_line;
										$raw_lines++;
										}
									elsif($output_strictness eq'all')
										{
										print RAW $output_line;
										$raw_lines++;
										}
									}
								}
							}
						}
					
					if($direction!=1)
						{
            #TODO: missing last alignment for minus strand
						# check for seed match in antisense orientation
						my $seed = substr($stretch,0,$perfect_5p-$clip);
						if($probes_rc{$seed}) {
							foreach$i(0..scalar(@{$probes_rc{$seed}})-1)
								{
								my $mapped_seq = $probes_rc{$seed}->[$i];
								my $p3 = substr($mapped_seq,0,-$perfect_5p);
								$count_internal_mm=0;
								$count_untemplate_3p=0;
								$count_clip_mm=0;
								foreach $pos($perfect_5p-$clip..$perfect_5p-1) 
									{
									if(substr($mapped_seq,length($p3)+$pos,1) ne substr($stretch,$pos,1)) 
										{
										$count_clip_mm++;
										}
									}
                $compares = "";
								foreach$pos(0..((length$p3)-$untemplate_3p)-1)
									{
									$g=substr($lagging_stretch,-1-$pos,1);
									$m=substr($p3,-1-$pos,1);
                  $compares .= " $g,$m";
									if($g ne $m) 
										{
										if(($g eq "G" && $m eq "A") || ($g eq "U" && $m eq "C")) 
											{
											$count_internal_mm += $GUpenalty;
                      $compares .= "GU";
											}
										else 
											{
											$count_internal_mm++;
                      $compares .= "m";
											}
										last if($count_internal_mm>$max_internal_mm);
										}
                    else { $compares .= "M" }
									}
								if($count_internal_mm<=$max_internal_mm)
									{
									$locus_seq=substr($stretch,0,$perfect_5p);
									if(length$p3>0)
										{
										$locus_seq=substr($lagging_stretch,-length$p3).$locus_seq;
										}
									#$locus_seq=~tr/AUGC/UACG/;
									#$locus_seq=reverse$locus_seq;
									
									$outside=0;
									foreach$pos(0..$untemplate_3p-1)
										{
										if(substr($lagging_stretch,-(length$p3)+$pos,1)eq'#') # check if 3' end of hit goes beyond start of the scaffold/chromosome sequence
											{
											$outside=1;
											last;
											}
										if(substr($p3,$pos,1)=~/[ATGCN]/&&substr($lagging_stretch,-(length$p3)+$pos,1) ne substr($p3,$pos,1))
											{
											$count_untemplate_3p++;
											}
										}
									next if($outside==1);
									$hit_mm=$count_internal_mm+$count_untemplate_3p+$count_clip_mm+$probes_rc_distance{$seed}->[$i];
									
                  $reverse_mapped_seq = reverse$mapped_seq;
                  $mapped_seed = reverse(substr($reverse_mapped_seq,0,length$seed));
                  $mapped_p3 = reverse(substr($reverse_mapped_seq,length$seed));

									$mapped_seq_rc=$mapped_seq;
									$mapped_seq=~tr/AUGC/UACG/;
									$mapped_seq_c=$mapped_seq;
									$mapped_seq=reverse$mapped_seq_c;
									$rc_position=($position-length$p3);

									my $probe_name = $probes_rc_info{$seed}->[$i];
									$output_line = "$replace_titles{$chromosome}\t$rc_position\t$locus_seq\t$probe_name\t$mapped_seq\t$hit_mm\t-\n";
									if($format eq 'eeland') 
										{
										$output_line = substr($output_line,0,-1)."\t$count_clip_mm\t$probes_rc_distance{$seed}->[$i]\t$count_internal_mm\t$count_untemplate_3p\n";
										if($printalignment)
											{
											$m = "";
											foreach$pos(0..length($locus_seq)-1) 
												{
													$l = substr($mapped_seq_rc,$pos,1); $p = substr($locus_seq,$pos,1);
													if($l eq $p) { $m .= "|" }
													elsif(($l eq "C" && $p eq "U") || ($l eq "A" && $p eq "G")) { $m .= ":" }
													else { $m .= " " }
												}
											$output_line .= "${locus_seq}\n${m}\n${mapped_seq_c}\n";
											}
										} 
									# do not output hit if output_strictness is 'best' and it is worse than a hit before
									unless(exists($best_hit{$mapped_seq}))
										{
										$best_hit{$mapped_seq}=$hit_mm;
										print RAW $output_line;
										$raw_lines++;
										}
									else
										{
										if($hit_mm<=$best_hit{$mapped_seq})
											{
											$best_hit{$mapped_seq}=$hit_mm;
											print RAW $output_line;
											$raw_lines++;
											}
										elsif($output_strictness eq'all')
											{
											print RAW $output_line;
											$raw_lines++;
											}
										}
									}
								}
							}
						}
					}
				}
			
			# elongate the sequence within the sliding windows
			else
				{
				$next_line=<GENOME>;
				$next_line=~s/\s*$//;
				if($next_line=~/^>/)
					{
					# finish processing previous reference sequence
					while(length$stretch>0)
						{
						search_seed_matches();
						sliding_window_increment();
						}
					
					$chromosome=$next_line;
					$chromosome=~s/>//;
					if($replace_titles_by_id==0)
						{
						$replace_titles{$chromosome}=$chromosome;
						}
					$position=-1;
					$stretch="##";
					$lagging_stretch="";
					}
				else
					{
					$processed_bp+=length$next_line;
					$stretch.=uc$next_line;
					$stretch=~tr/atgcuT/AUGCUU/;
					$concatenate_index=1;
					}
				}
			}
		$t_end=time;
		$time_used=$t_end-$t_start;
    if($verbosity==2)
			{
			print"\nPROBE:.. $probe\nREF:.... $chromosome\nPOS:.... $position\nScanned: $processed_bp bp\nSpeed:.. $speed bp/s\nEstimated remaining time: $est_remaining seconds.\n\nDone. Time used for mapping: $time_used seconds.\n\n"
			}
		elsif($verbosity==1)
			{
			print" done.\nTime used for mapping: $time_used seconds.\n\n"
			}

		if($format eq "eeland") 
			{
			if($outfile!~/^\s*$/) 
				{
				rename("$probe.raw","$outfile")||print"\nCould not rename RAW file '$probe.raw' -> '$outfile'.\n$!\n\n";
				}
			else
				{
				rename("$probe.raw","$probe.map")||print"\nCould not rename RAW file '$probe.raw' -> '$probe.map'.\n$!\n\n";
				}
				exit 0;
			}

		# FILTER AND/OR SORT
			if($verbosity == 2)
				{
				print "Filtering/Sorting raw file:\n\n0%                                50%                                100%\n|----------------------------------|----------------------------------|\n";
				}
			elsif($verbosity == 1)
				{
					print "Filtering/Sorting raw file:\n\n0%                                50%                                100%\n|----------------------------------|----------------------------------|\n";
				}


		$t_start=time;
		open(RAW,"$probe.raw")||die print"\nCould not open RAW file $probe.raw.\n$!\n\n";
		open(MAP,">$probe.map")||die print"\nCould not create MAP file $probe.map.\n$!\n\n";
		
		if($format eq 'sam')
			{
			print MAP"\@HD VN: 1.5 SO:coordinate\n";
			foreach$rname(@sam)
				{
				print MAP"\@SQ SN:$rname LN:$sam{$rname}\n";
				}
			}
		
		@unsorted=();
		$buffer_max=100000;
		$buffer_release=90000;
		$prev_chromosome="#";
		$printed_dots=0;
		$processed_raw_lines=0;
		while(<RAW>)
			{
			# show progress
			if($verbosity!=0)
				{
				$processed_raw_lines++;
				if($processed_raw_lines>=$raw_lines/71)
					{
					print".";
					$printed_dots++;
					$processed_raw_lines=0;
					}
				}
			
			if($_=~/^#/&&$format ne 'sam')
				{
				print MAP$_;
				}
			else
				{
				@d=split("\t",$_);
				if($prev_chromosome ne "#"&&$d[0] ne $prev_chromosome)
					{
					sort_map();
					sub sort_map
						{
						@sorted=sort@unsorted;
						foreach$element(@sorted)
							{
							$element=~s/^.{13}//; # delete the index tag from entry (index tag is used for sorting only)
							print MAP$element;
							}
						undef@sorted;
						@sorted=();
						undef@unsorted;
						@unsorted=();
						}
					}
				elsif(@unsorted==$buffer_max)
					{
					@sorted=sort@unsorted;
					foreach$i(0..$buffer_release-1)
						{
						$sorted[$i]=~s/^.{13}//; # delete the index tag from entry (index tag is used for sorting only)
						print MAP$sorted[$i];
						}
					undef@unsorted;
					@unsorted=();
					foreach$i($buffer_release..$buffer_max-1)
						{
						push(@unsorted,$sorted[$i]);
						}
					undef@sorted;
					@sorted=();
					}
				
				if($d[5]==$best_hit{$d[4]}||$output_strictness eq'all')
					{
					$index="000000000000".$d[1];
					$index=substr($index,-12);
					if($format eq'compact') # chromosome | coordinate | strand | length | divergence | title
						{
						$length_hit=length$d[4];
						if($d[5]>0)
							{
							@gen_seq=split('',$d[2]);
							@hit_seq=split('',$d[4]);
							$seq_div="";
							foreach$pos(1..$length_hit)
								{
								if($gen_seq[$pos-1]ne$hit_seq[$pos-1])
									{
									$seq_div.="$pos$hit_seq[$pos-1]";
									}
								}
							}
						else
							{
							$seq_div="-";
							}
						$d[6]=~s/\s*$//;
						push(@unsorted,"$index\t$d[0]\t$d[1]\t$d[6]\t$length_hit\t$seq_div\t$d[3]\n");
						}
					elsif($format eq'eeland')
						{
						push(@unsorted,"$index\t$d[0]\t$d[1]\t$d[2]\t$d[3]\t$d[4]\t$d[5]\t$d[6]\t$d[7]\t$d[8]\t$d[9]\t$d[10]");
						}
					elsif($format eq'eland')
						{
						push(@unsorted,"$index\t$d[0]\t$d[1]\t$d[2]\t$d[3]\t$d[4]\t$d[5]\t$d[6]");
						}
					elsif($format eq'sam')
						{
						# generate CIGAR string
						$cigar_pre="";
						$stretch=0;
						foreach$pos(0..(length$d[4])-1)
							{
							if(substr($d[2],$pos,1) eq substr($d[4],$pos,1))
								{
								$cigar_pre.="=";
								}
							else
								{
								$cigar_pre.="X";
								}
							}
						$cigar="";
						while(1)
							{
							if($cigar_pre=~s/^=+//)
								{
								$n=length$&;
								$cigar.="$n"."=";
								}
							elsif($cigar_pre=~s/^X+//)
								{
								$n=length$&;
								$cigar.="$n"."X";
								}
							else
								{
								last;
								}
							}

						if($d[6]=~/-/)
							{
							$d[4]=reverse$d[4];
							$d[4]=~tr/AUGC/UACG/;
							}
						push(@unsorted,"$index\t$d[3]\t0\t$d[0]\t$d[1]\t255\t$cigar\t$d[3]\t*\t0\t0\t$d[4]\t*\n");
						}
					if($d[1]>$max_coordinate)
						{
						$max_coordinate=$d[1];
						}
					}
				$prev_chromosome=$d[0];
				}
			}
		sort_map();
		close RAW;
		close MAP;
		unlink"$probe.raw"||print"\nCould not remove RAW file $probe.raw.\n$!\n\n";
		$t_end=time;
		$time_used=$t_end-$t_start;
		if($verbosity!=0)
			{
			foreach($printed_dots..70)
				{
				print".";
				}
			print"\nTime used for filtering/sorting: $time_used seconds.\n";
			}
		if($outfile!~/^\s*$/)
			{
			rename("$probe.map","$outfile")||print"\nCould not rename MAP file '$probe.map' -> '$outfile'.\n$!\n\n";
			}
		exit 0;
		}
	else
		{
		die print"\nCould not fork mapping process: $!\n";
		}
	}
foreach(@childs)
	{
	$tmp=waitpid($_,0);
	}
exit;
