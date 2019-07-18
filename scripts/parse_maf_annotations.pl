#!/usr/bin/perl
use strict;
use warnings;

my ($gtf, $maf, $prelim, $dist, $id, $outpath)=@ARGV[0 .. $#ARGV];

########################################################
#
# $gtf: annotation file in GTF format
# $maf: alignment file in MAF format
# $prelim: output file ("potential_missing_exons.out") from 2_parse_maf_transcripts_04.pl
# $dist: number of nucleotide to expand gene model definition
# $id: experiment name
# $outpath: path to the output directory
#
########################################################



my %models;
my %transcripts;

open(OUT4, ">$outpath/exons.bed");
FLATTEN_GTF($gtf);
print "DONE GTF\n";
close OUT4;

my %mapping;
open(IN, "$prelim")||die"IN $prelim\n";
while(<IN>){
	chomp;
	my @split=split /\t/, $_;
	if($models{$split[1]}){
		@{$mapping{$split[0]}}=@{$models{$split[1]}};
	}
}
close IN;

open(OUT, ">$outpath/reads_genome.bed");
PARSE_MAF($maf);
close OUT;

system("intersectBed -a $outpath/reads_genome.bed -b $outpath/exons.bed -v > $outpath/reads_missing_exons.bed");

########################################################
#Collect exonic positions and genomic positions of genes to parse the alignments and identify potential novel exons
########################################################
sub FLATTEN_GTF{
	my ($file)=(@_);
	my %temp_trans;
	my %temp_chr;
	open(FILE, "$file")||die"FILE $file \n";
	while (<FILE>){
		if($_=~/\texon\t/ || $_=~/\tgene\t/){
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
			if($split[2]=~/exon/){
				$transcripts{$transid}=$geneid;
				push @{$temp_chr{$split[0]}}, [$split[0], $split[3], $split[4], $geneid, $transid];
			}
			elsif($split[2]=~/gene/){
				my $start=$split[3]-$dist;
				if($split[3] <= $dist){
					$start=1;
				}
				@{$models{$geneid}}=($split[0], $start, $split[4]+$dist);
			}
		}
	}
	close FILE;
	foreach my $ele (keys %temp_chr){
		@{$temp_chr{$ele}}=sort{$a->[1]<=>$b->[1]||$a->[2]<=>$b->[2]}@{$temp_chr{$ele}};
		for my $i (0 .. $#{$temp_chr{$ele}}){
			print OUT4 join ("\t", @{$temp_chr{$ele}[$i]}),"\n";
		}
	}
}
########################################################
# Parse the alignments to the genomes, identify the gene present in the genomic region to which a read is aligned, 
########################################################
sub PARSE_MAF{
	my ($file)=(@_);
	my $score;
	my $c=0;
	my $first=0;
	my $test=0;
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
				$s1=$temp1[2];
				$e1=($s1+$temp1[3])-1;
				$size1=$temp1[5];
				$in[1]=~s/[\s]+/\t/g;
				my @temp2=split /\t/, $in[1];
				$name=$temp2[1];
				$s2=$temp2[2];
				$e2=($s2+$temp2[3])-1;
				$size2=$temp2[5];
# assess the location of the genomic alignment, test if conform to the mapping to the transcriptome: the genomic location of the alignment has to overlap with the genomic position of the gene
				if($mapping{$name}){
					if($mapping{$name}[0]=~/$trans_name\b/ && ($s1 >= $mapping{$name}[1] && $s1<=$mapping{$name}[2])){ # test if correct chromosome and positions within the gene coordinates
						print OUT join ("\t", $trans_name, $s1, $e1, $name, $s2, $e2, $id),"\n";
					}
				}
				else{
					print "$name\n";
				}
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
}
