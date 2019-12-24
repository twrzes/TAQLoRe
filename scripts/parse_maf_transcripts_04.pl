#!/usr/bin/perl
use strict;
use warnings;

my ($maf, $gtf, $min_prop, $min_prop_align, $min_insert, $min_exon_distance, $threads, $fasta_2D, $genomeDB, $outpath)=@ARGV[0 .. $#ARGV];

#######################
#
# $maf: last alignment file of the reads to the transcripts
# $gtf: annotation file in GTF format
# $min_prop: minimum proportion of the read aligned to trigger the potential exon search
# $min_prop_align: minimum proportion of the read aligned to report the alignment
# $min_insert: minimum insert size within the read
# $min_exon_distance: minimum distance to existing exons
# $threads: number of threads for the last alignments
# $fasta_2D: reads file in fasta format
# $genomeDB: genome database for the last alignment
# $outpath: path for output 
#
######################



my %models;
my %loci;
my %transcripts;
my %trans_models;
open(OUT4, ">$outpath/exons.bed");
#parse the annotation file
FLATTEN_GTF($gtf);
print "DONE GTF\n";
close OUT4;

my %align;
my $temp_pos;
my $genome_s;
my $genome_e;
my $test=0;


my %selected_reads;
open(OUT1, ">$outpath/potential_missing_exons.out");
open(OUT2, ">$outpath/potential_microexons.out");
open(OUT3, ">$outpath/read_transcripts.out");

#Parse the maf alignment
PARSE_MAF($maf);
close OUT1;
close OUT2;
close OUT3;


# prepare a fasta file with the reads of interest and align the reads to the genome

if(scalar keys %selected_reads > 0){
	open(OUT5, ">$outpath/selected_reads_potential_exons.fa");
	GET_SEQ($fasta_2D);
	close OUT5;
	system("lastal -P $threads $genomeDB $outpath/selected_reads_potential_exons.fa > $outpath/selected_reads_potential_exons.maf \n")==0||die"FAIL LAST $genomeDB \n";
}
########################################################
#parse the annotation file to generate a bed file of the exon positions and an hash associating positions within transcripts to genomic positions 
########################################################
sub FLATTEN_GTF{
	my ($file)=(@_);
	my %temp_trans;
	my %temp_chr;
	open(FILE, "$file")||die"FILE $file \n";
	while (<FILE>){
		if($_=~/\texon\t/){
			chomp;
			my $line=$_;
			$line=~/gene_id \"([0-9a-zA-Z_]+)\";/;
			my $geneid=$1;
			$line=~/transcript_id \"([0-9a-zA-Z_]+)\";/;
			my $transid=$1;
			$_=~s/\"//g;
			$_=~s/;//g;
			$_=~s/[\s]+/\t/g;
			my @split=split /\t/, $_;
			$transid=$split[11];
			$geneid=$split[9];
			$transcripts{$transid}=$geneid;
			push @{$temp_trans{$transid}}, [$split[0], $split[3], $split[4], $split[6], $transid];
			push @{$temp_chr{$split[0]}}, [$split[0], $split[3], $split[4], $split[6], $geneid, $transid];
		}
	}
	close FILE;

	foreach my $ele (keys %temp_chr){
		@{$temp_chr{$ele}}=sort{$a->[1]<=>$b->[1]||$a->[2]<=>$b->[2]}@{$temp_chr{$ele}};
		for my $i (0 .. $#{$temp_chr{$ele}}){
			print OUT4 join ("\t", @{$temp_chr{$ele}[$i]}),"\n";
		}
	}
	foreach my $cdna (keys %temp_trans){
		GENOME_POSITIONS(\@{$temp_trans{$cdna}});
	}
	print "DONE GENOME POSITIONS\n";
}
########################################################
# Associate positions within a transcript (1 .. n) to genomic positions
########################################################
sub GENOME_POSITIONS{
	my ($array)=(@_);
	@{$array}=sort{$a->[1]<=>$b->[1]||$a->[2]<=>$b->[2]}@{$array};
	if($array->[0][3]=~/-/){
		@{$array}=reverse @{$array};
	}
	my $s;
	my $e;
	for my $o (0 .. $#{$array}){
		if($o==0){
			$s=1;
			$e=$s+($array->[$o][2]-$array->[$o][1]);
			push @{$trans_models{$array->[0][4]}}, [$array->[$o][0], $s, $e, $array->[$o][1], $array->[$o][2], $array->[$o][3]];
		}
		else{
			$e=$s+($array->[$o][2]-$array->[$o][1]);
			push @{$trans_models{$array->[0][4]}}, [$array->[$o][0], $s, $e, $array->[$o][1], $array->[$o][2], $array->[$o][3]];
		}
		$s=$e+1;
	}

}
#########################################################
# Parse the maf alignment to identify insertions within the reads of size >= $min_insert that could indicate potential novel exons
#########################################################
sub PARSE_MAF{
	my ($file)=(@_);
	my $score;
	my $c=0;
	my $previous;
	my $first=0;
	my $test=0;
	my @inserts;
	my $s1;
	my $e1;
	my $size1;
	my $trans_name;
	my $name;
	my $s2;
	my $e2;
	my $size2;
	my $test2=0;
	my @in;
	my %micro;
	open(SEQ, "$file")||die"SEQ $file\n";
	while(<SEQ>){
		if($c==0){
			if($_=~/^a/){
				chomp;
				$_=~s/[\s]+/\t/g;
				my @split=split /\t/, $_;
				$score=$split[1];
				$score=~s/score=//;
				$c=1;
				$_="";
			}
		}
		if($c==1){
			$test=0;
			$test2=0;
			unless($_=~/^a/){
				if($_=~/^s/){
					push @in, $_;
				}
			}
			if($_=~/^a/ or eof(SEQ)){
				$in[0]=~s/[\s]+/\t/g;
				my @temp1=split /\t/, $in[0];
				$trans_name=$temp1[1];
				if($transcripts{$trans_name}){
					$s1=$temp1[2];
					$e1=($s1+$temp1[3])-1;
					$size1=$temp1[5];
					$test2=0;
					while($temp1[6]=~/-{$min_insert,}/g){#test for the possibility of microexons
						push @inserts, join("..", $-[0], $+[0]);
						$test2=1;
					}
					if(scalar @inserts == 0){
						push @inserts, "NA";
					}
					$in[1]=~s/[\s]+/\t/g;
					my @temp2=split /\t/, $in[1];
					$name=$temp2[1];
					$s2=$temp2[2];
					$e2=($s2+$temp2[3])-1;
					$size2=$temp2[5];
					if($test2==1){
						push @{$micro{$name}{$transcripts{$trans_name}}}, [$name, $transcripts{$trans_name}, $trans_name, join(";", @inserts), $temp2[6]];
					}
					push @{$align{$name}}, [$score, $transcripts{$trans_name}, $trans_name, $s1, $e1, $size1, $s2, $e2, $size2, $name];
				}
				@inserts=();
				@in=();
				$s1="";
				$e1="";
				$size1="";
				$trans_name="";
				$name="";
				$s2="";
				$e2="";
				$size2="";
				$test2=0;
				if($_=~/^a/){
					chomp;
					$_=~s/[\s]+/\t/g;
					my @split=split /\t/, $_;
					$score=$split[1];
					$score=~s/score=//;
					$c=1;
					$_="";
				}
				else{
					$c=0;
				}
			}
		}
	}
	close SEQ;
	foreach my $read (keys %align){
		@{$align{$read}}=sort{$b->[0]<=>$a->[0]}@{$align{$read}};
		if($micro{$read}{$align{$read}[0][1]}){
			print OUT2 join ("\t", @{$micro{$read}{$align{$read}[0][1]}[0]}),"\n";
		}
		PARSE_ALIGN(\@{$align{$read}});
	}
}
#########################################################
# Parse the quality of the alignment (more than $min_prop_align of the read needs to be aligned for it to be recoreded in the $outpath/read_transcripts.out output file and more than $min_prop of the read to be aligned for the search to continue)
#########################################################
sub PARSE_ALIGN{
	my ($array)=(@_);
	@{$array}=sort{$b->[0]<=>$a->[0]}@{$array};
	my $target=$array->[0][1];
	my $target_trans=$array->[0][2];
	my $start=$array->[0][6];
	my $end=$array->[0][7];
	if((($end-$start)+1)/$array->[0][8] >= $min_prop_align){
		print OUT3 join ("\t", $array->[0][9], $array->[0][2], $array->[0][3], (($end-$start)+1)/$array->[0][8], (($array->[0][4]-$array->[0][3])+1)/$array->[0][5]),"\n";
	}
	if((($end-$start)+1)/$array->[0][8] > $min_prop){

		my @positions;
		for my $j (0 .. $#{$array}){
			if($array->[$j][1]=~/$array->[0][1]/){
				push @positions, [@{$array->[$j]}];
			}
		}
		if(scalar @positions > 1){
			my %done;
			my @models;
			@positions=sort{$a->[6]<=>$b->[6]||$a->[7]<=>$b->[7]}@positions;
			for my $k (0 .. $#positions){
				if(!$done{$k}){
					my $s3=$positions[$k][6];
					my $e3=$positions[$k][7];
					my $target_s=$positions[$k][3];
					my $target_e=$positions[$k][5];
					my $name_s=$positions[$k][2];
					my $name_e=$positions[$k][2];
					for (my $l=$k+1; $l<=$#positions; $l++){
						my $s4=$positions[$l][6];
						my $e4=$positions[$l][7];
						if($s4>$e3){
							last;
						}
						elsif(($s3>=$s4 && $e3<=$e4) || ($s3<=$s4 && $e3>=$e4) || ($s3>=$s4 && $s3<$e4 && $e3>=$e4) || ($s3<=$s4 && $e3 > $s4 && $e3<=$e4)){
							if($s3>$s4){
								$s3=$s4;
								$target_s=$positions[$l][3];
								$name_s=$positions[$l][2];
							}
							if($e3<$e4){
								$e3=$e4;
								$target_e=$positions[$l][4];
								$name_e=$positions[$l][2];
							}
							$done{$l}=1;
						}
					}
					$done{$k}=1;
					push @models, [$s3, $e3, $target_s, $target_e, $name_s, $name_e, $positions[$k][1], $positions[$k][9]];
				}
			}
			if(scalar @models > 1){
				PARSE_MODELS(\@models);
			}
			@models=();
		}
	}
}
#########################################################
# test for contiguity
#########################################################
sub PARSE_MODELS{
	my ($reads)=(@_);
	for my $m (0 .. $#{$reads}){
		if($m>0){
			if($reads->[$m][0]-$reads->[$m-1][1] > $min_exon_distance){
				$test=1;
				$temp_pos="NA";
				GET_POS($reads->[$m-1][3], $reads->[$m-1][5]);
				$genome_s=$temp_pos;
				$temp_pos="NA";
				GET_POS($reads->[$m][2], $reads->[$m][4]);
				$genome_e=$temp_pos;
				$temp_pos="";
				print OUT1 join ("\t", $reads->[$m][7], $reads->[$m][6], $reads->[$m][4], $genome_s, $genome_e, $reads->[$m-1][1], $reads->[$m][0]),"\n";
				$selected_reads{$reads->[$m][7]}=1;
				$genome_s="";
				$genome_e="";
			}
		}
	}
}
#########################################################
#########################################################
sub GET_POS{
	my($p, $q)=(@_);
	for my $s (0 .. $#{$trans_models{$q}}){
		if($p >= $trans_models{$q}[$s][1] && $p <= $trans_models{$q}[$s][2]){
			$temp_pos=$trans_models{$q}[$s][3]+($p-$trans_models{$q}[$s][1]);
		}
	}
}
#########################################################
# Collect the nucleotide sequencine of the reads of interest
#########################################################
sub GET_SEQ{
	my ($gene)=(@_);
	my $seq="";
	my $head;
	my $c=0;
	open(SEQ, "$gene")||die"SEQ $gene\n";
	while(<SEQ>){
		if($c==0){
			if($_=~/>/){
				chomp;
				$_=~s/[\s]+/\t/g;
				$_=~s/>//;
				my @temp4=split /\t/, $_;
				if($selected_reads{$temp4[0]}){
					$head=$temp4[0];
					$c=1;
					$_="";
				}
				else{
					$c=0;
				}
			}
		}
		if($c==1){
			unless($_=~/>/){
				$_=~s/\s+//g;
				$seq.=$_;
			}
			if($_=~/>/ or eof(SEQ)){
				print OUT5 join ("\n", ">$head", $seq),"\n";
				$head="";
				$seq="";
				if($_=~/>/){
					chomp;
					$_=~s/[\s]+/\t/g;
					$_=~s/>//;
					my @temp4=split /\t/, $_;
					if($selected_reads{$temp4[0]}){
						$head=$temp4[0];
						$c=1;
						$_="";
					}
					else{
						$c=0;
					}
				}
				else{
					$c=0;
				}
			}
		}
	}
	close SEQ;
}
