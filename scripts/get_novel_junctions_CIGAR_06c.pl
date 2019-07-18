#!/usr/bin/perl
use strict;
use warnings;

my ($exon_positions, $align, $min_cov, $min_exon, $outpath)=@ARGV;

########################################################
#
# $exon_positions: exon coordinates in metagene
# $align: sam file
# $min_cov: miniumn number of nucleotide overlap
# $min_exon: minimum propotion of the exon covered 
# $outpath: path to the output directory
#
########################################################

my @meta_models;
my $temp_pattern;
my %read_patterns;
my %read_positions;
my %gene_positions;
my %full;
my $size_covered;
my %exon_coverage;
my %non_coding;

my %splicing_pattern;
my $file_name;
my @temp6=split /\//, $align;
$file_name=$temp6[$#temp6];
$file_name=~s/.gmap_metagene_all_exons.sam//;

# Collect all the positions both within the metagene and genomic positions

META_ANNOTATIONS($exon_positions);
print "DONE META ANNOTATIONS\n";
my $test_cov=0;

# retrieve the mapping pattern of the reads
PARSE_SAM($align);
print "DONE SAM\n";

open (OUT, ">$outpath/$file_name\_splicing_patterns_cds.tmp");
open(OUT2, ">$outpath/$file_name\_genomic_positions_reads_cds.bed");
my $count_reads=0;
my %done_read;
foreach my $k (keys %gene_positions){
	my @all_positions;
	if(scalar @{$gene_positions{$k}} > 1){ #merging of reads found in multiple entries of the sam file
		@{$full{$k}}=sort{$a->[0]<=>$b->[0]||$a->[1]<=>$b->[1]}@{$full{$k}};
		my %done;
		for my $i (0 .. $#{$full{$k}}){
			my $s1=$full{$k}[$i][0];
			my $e1=$full{$k}[$i][1];
			for (my $j=$i+1; $j<=$#{$full{$k}}; $j++){
				my $s2=$full{$k}[$j][0];
				my $e2=$full{$k}[$j][1];
				if(($s1>=$s2 && $e1<=$e2) || ($s1<=$s2 && $e1>=$e2) || ($s1>=$s2 && $s1<$e2 && $e1>=$e2) || ($s1<=$s2 && $e1 > $s2 && $e1<=$e2)){
					print join ("\t", $s1, $e1, $s2, $e2),"\n";
					$done{$i}=1;
					$done{$i+1}=1;
				}
				else{
					for my $j (0 .. $#{$gene_positions{$k}[$i]}){
						my @var=split /_/, $gene_positions{$k}[$i][$j];
						push @all_positions, [$var[0], $var[1], $var[2]];
					}
					$done{$i}=1;
				}
			}
			if($i==$#{$full{$k}}){
				if(!$done{$i}){
					for my $j (0 .. $#{$gene_positions{$k}[$i]}){
						my @var=split /_/, $gene_positions{$k}[$i][$j];
						push @all_positions, [$var[0], $var[1], $var[2]];
					}
				}
			}
		}
	}
	else{
		for my $i (0 .. $#{$gene_positions{$k}[0]}){
			my @var=split /_/, $gene_positions{$k}[0][$i];
			push @all_positions, [$var[0], $var[1], $var[2]];
		}
	}
	if(scalar @all_positions > 0){
		$count_reads++;
		GET_PATTERN(\@all_positions);
		$read_patterns{$k}=$temp_pattern;
		if(!$done_read{$k}){
			GET_GENOMIC_POSITIONS($k, \@all_positions, $meta_models[0][0]);
			$done_read{$k}=1;
		}
	}
}
close OUT;
close OUT2;

open(OUT4, ">$outpath/$file_name\_exon_counts_$min_exon.out");
for my $n (0 .. $#meta_models){
	if($exon_coverage{$n+1}){
		 print OUT4 join ("\t", $n+1, $exon_coverage{$n+1}),"\n";
	}
	else{
		print OUT4 join ("\t", $n+1, "0"),"\n";
	}
}
close OUT4;

print join ("\t", "TOTAL NUMBER OF READS:", $count_reads),"\n";
open(OUT3, ">$outpath/$file_name\_count_patterns.out");
foreach my $z (keys %splicing_pattern){
	print OUT3 join ("\t", $z, scalar @{$splicing_pattern{$z}}, join (",", @{$splicing_pattern{$z}})),"\n";
}
close OUT3;

########################################################
# Collect all the positions of the exons both within the metagene and genomic positions
########################################################
sub META_ANNOTATIONS{
	my ($file)=(@_);
	my $temp_count=0;
	open(FILE, "$file")||die"FILE $file \n";
	while (<FILE>){
		chomp;
		$temp_count++;
		my @split=split /\t/, $_;
		push @meta_models, [@split];
		if($_=~/UTR/){
			$non_coding{$temp_count}=1;
		}
	}
	close FILE;
}
########################################################
# Parse the alignment file to retrieve the mapping positions, collect the CIGAR string and parse it
########################################################
sub PARSE_SAM{
	my ($file)=(@_);
	open(SEQ, "$file")||die"SEQ $file\n";
	while(<SEQ>){
		if($_!~/^\@/){
			chomp;
			my @temp=split /\t/, $_;
			if($temp[5]!~/\*/){
				if($_=~/XS:A:-/){
					GET_POSITIONS($temp[0], $temp[3], $temp[5], $temp[2], length $temp[9], "-1");
				}
				else{
					GET_POSITIONS($temp[0], $temp[3], $temp[5], $temp[2], length $temp[9], "1");
				}
			}
		}
	}
	close SEQ;
}
########################################################
# Parse the CIGAR string to identify the mapping positions
########################################################
sub GET_POSITIONS{
	my ($head, $query, $cigar, $name, $read_size, $strand)=(@_);
	# adapted from https://davetang.org/muse/2011/01/28/perl-and-sam/
	my $position=$query;
	my $query_s=$query;
	my $query_e=0;
	my $read_s=0;
	my $read_e=0;
	my @query_positions;
	my $full_s=$query;
	my $full_e=0;
	my $temp_cigar=$cigar;
	my $counter1=0;
	while ($temp_cigar !~ /^$/){
		if ($temp_cigar =~ /^([0-9]+[MIDSNH])/){
			$counter1++;
			my $cigar_part=$1;
			$temp_cigar =~ s/$cigar_part//;
		}
	}
	my $counter2=0;
	while ($cigar !~ /^$/){
		if ($cigar =~ /^([0-9]+[MIDSNH])/){
			$counter2++;
			my $cigar_part = $1;
			if ($cigar_part =~ /(\d+)M/){
				$query_s=$query_e+1;
				$query_e=($query_s + $1)-1;
				$read_s=$read_e+1;
				$read_e=($read_s+$1)-1;
				push @query_positions, join ("_", $query_s, $query_e, $name, $read_s, $read_e, $1, $cigar_part);
				$full_e=$query_e;
				$size_covered+=$1;
			}
			elsif ($cigar_part =~ /(\d+)I/){
				$read_s=$read_e+1;
				$read_e=($read_s + $1) - 1;
				$query_s=$query_e;
				$query_e=$query_s;
			}
			elsif ($cigar_part =~ /(\d+)S/){
				$read_e=$read_s + $1;
					$query_e=$query_s-1;
			}
			elsif ($cigar_part =~ /(\d+)H/){
					$read_e=$read_s + $1;
					$query_e=$query_s-1;
			}
			elsif ($cigar_part =~ /(\d+)D/){
				$read_s=$read_e;
				$read_e=$read_s;
				$query_s=$query_e+1;
				$query_e=($query_s + $1)-1;
				$full_e=$query_e;
				$size_covered+=$1;
			}
			elsif ($cigar_part =~ /(\d+)N/){
				$read_s=$read_e;
				$read_e=$read_s;
				$query_s=$query_e+1;
				$query_e=($query_s + $1)-1;
				$full_e=$query_e;
				$size_covered+=$1;
			}
			$cigar =~ s/$cigar_part//;
		}
	}
	#positions of a read within a transcript use of a muti-dimentional array to take into account reads split in multiple entries of the sam file
	push @{$full{$head}}, [$full_s, $full_e];
	#mapped size and mapped positions within a read
	push @{$read_positions{$head}}, [($read_e-$read_s)+1, $read_s, $read_e];
	#projection of the reads to the transcript positions
	push @{$gene_positions{$head}}, [@query_positions];
	@query_positions=();
}
########################################################
# Report the mapping pattern 
########################################################
sub GET_PATTERN{
	my ($query)=(@_);
	my %done4;
	my @temp_junc;
	@{$query}=sort{$a->[0]<=>$b->[0]||$a->[1]<=>$b->[1]}@{$query};
	for my $k(0 .. $#{$query}){
		my $s3=$query->[$k][0];
		my $e3=$query->[$k][1];
		for my $l (0 .. $#meta_models){
			my $s4=$meta_models[$l][0];
			my $e4=$meta_models[$l][1];
			if(($s3>=$s4 && $e3<=$e4) || ($s3<=$s4 && $e3>=$e4) || ($s3<=$s4 && $e3>=$s4 && $e3<=$e4)|| ($s3>=$s4 && $s3<=$e4 && $e3>=$e4)){
				$test_cov=0;
				if($s3<=$s4 && $e3>=$e4){
					$test_cov=1;
				}
				else{
					TEST_COV($s3, $e3, $s4, $e4);
				}
				if($test_cov==1){
					if(!$done4{$l+1}){
						push @temp_junc, $l+1;
						$done4{$l+1}=1;
					}
				}
				$test_cov=0;
			}
			if($s4>$e3){
				last;
			}
		}
	}
	$temp_pattern=join ("_", @temp_junc);
}
########################################################
########################################################
sub TEST_COV{
	my ($s5, $e5, $s6, $e6)=(@_);
	if($s5==$s6 && $e5==$e6){
		$test_cov=1;
	}
	elsif($s5<=$s6 && $e5>=$s6 && $e5<=$e6){
		if(($e5-$s6) >= $min_cov ){
			$test_cov=1;
		}
	}
	elsif($s5>=$s6 && $s5<=$e6 && $e5>=$e6){
		if(($e6-$s5)  >= $min_cov ){
			$test_cov=1;
		}
	}
}
########################################################
# Retrieve the genomic positions of the different alignment blocks
########################################################
sub GET_GENOMIC_POSITIONS{
	my ($name2, $array_pos, $chr)=(@_);
	my %exons_pos;
	my %exons_blocks;
	my @g_pos;
	my @genome_pos;
	my @blocks;
	@{$array_pos}=sort{$a->[0]<=>$b->[0]||$a->[1]<=>$b->[1]}@{$array_pos};
	@meta_models=sort{$a->[1]<=>$b->[1]||$a->[2]<=>$b->[2]}@meta_models;
	my %exons;
	for my $i (0 .. $#{$array_pos}){
		my $s1=$array_pos->[$i][0];
		my $e1=$array_pos->[$i][1];
		my $g_start=$meta_models[0][2];
		for my $j (0 .. $#meta_models){
			my $s2=$meta_models[$j][0];
			my $e2=$meta_models[$j][1];
			if($s1 >= $s2 && $s1 <= $e2){
				if($e1 <= $e2){
					push @{$exons_blocks{$j+1}}, ($e1-$s1)+1;
					$exons{$j+1}+=($e1-$s1)+1;
					my $temp_s=$meta_models[$j][2]+($s1-$s2);
					push @{$exons_pos{$j+1}}, $temp_s-$meta_models[0][2];
					last;
				}
				elsif($e1 > $e2){
					push @{$exons_blocks{$j+1}}, ($e2-$s1)+1;
					my $temp_s=$meta_models[$j][2]+($s1-$s2);
					push @{$exons_pos{$j+1}}, $temp_s-$meta_models[0][2];
					$exons{$j+1}+=($e2-$s1)+1;
					for (my $n=$j+1; $n<= $#meta_models; $n++){
						if($e1 >= $meta_models[$n][1]){
							push @{$exons_blocks{$n+1}}, ($meta_models[$n][1]-$meta_models[$n][0])+1;
							push @{$exons_pos{$n+1}}, $meta_models[$n][2]-$meta_models[0][2];
							$exons{$n+1}+=($meta_models[$n][1]-$meta_models[$n][0])+1;
						}
						elsif($e1 >= $meta_models[$n][0] && $e1 <= $meta_models[$n][1]){
							push @{$exons_blocks{$n+1}}, ($e1-$meta_models[$n][0])+1;
							push @{$exons_pos{$n+1}}, $meta_models[$n][2]-$meta_models[0][2];
							$exons{$n+1}+=($e1-$meta_models[$n][0])+1;
							last;
						}
					}
					last;
				}
				last;
			}
		}
	}

	my @pat;
	foreach my $nb (sort{$a<=>$b}keys %exons){
		if($exons{$nb}/(($meta_models[$nb-1][1]-$meta_models[$nb-1][0])+1)>= $min_exon){
			$exon_coverage{$nb}++;
			if(!$non_coding{$nb}){
				push @pat, $nb;
				for my $z (0 .. $#{$exons_blocks{$nb}}){
					push @blocks, $exons_blocks{$nb}[$z];
					push @genome_pos, $exons_pos{$nb}[$z];
				}
			}
		}
	}


	print OUT join ("\t", $name2, join ("_", @pat)), "\n";
	my $name3= join ("_", @pat);
	push @{$splicing_pattern{$name3}}, $name2; 

	@pat=();
	my @g_pos2;
	for my $m (0 .. $#genome_pos){
		push @g_pos2, $genome_pos[$m]-$genome_pos[0];
	}
	if(defined $genome_pos[0]){
		my $g_s=($meta_models[0][2]-1)+$genome_pos[0];
		my $g_e=($meta_models[0][2]-1)+($genome_pos[$#genome_pos]+$blocks[$#blocks]);
		print OUT2 join ("\t", "chr12", $g_s, $g_e, $name2, 0, "+", $g_s, $g_e, 0, scalar @blocks, join(",", @blocks), join(",", @g_pos2)),"\n";
	}
	else{
		print "$name2\n";
	}
	@blocks=();
	@genome_pos=();
}
