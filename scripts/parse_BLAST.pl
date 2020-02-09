#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
##############################################################################
##	Jose Fco. Sanchez Herrero, 25/09/2018 jfsanchezherrero@ub.edu			##
##	Updated: February 2020
##############################################################################
##
## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3820096/
##
## An introduction to sequence similarity ("homology") searching:
##
## BLAST, FASTA, SSEARCH, and other commonly used similarity searching programs
## produce accurate statistical estimates that can be used to reliably infer homology. Searches
## with protein sequences (BLASTP, FASTP, SSEARCH,) or translated DNA sequences
## (BLASTX, FASTX) are preferred because they are 5- to 10-fold more sensitive than
## DNA:DNA sequence comparison. The 30% identity rule-of-thumb is too conservative;
## statistically significant [E() < 10−6 – 10−3] protein homologs can share less than 20%
## identity. E()-values and bit scores (bits >50) are far more sensitive and reliable than
## percent identity for inferring homology.
##
#######################################################################################

my $blast_file = $ARGV[0];
my $output_name = $ARGV[1];
my $similarity = $ARGV[2];
my $aln = $ARGV[3];

my $stop_signal = $ARGV[4];
#my $evalue=$ARGV[4];

if (!@ARGV) {print "Usage:\nperl $0 blast_results output [aln_len similarity] [stop_signal]\n";exit();}

## default
if (!$similarity) {$similarity=80;}
if (!$aln) {$aln=80;}

my $evalue = 1e-05;
## protein:protein aln = evalue 1e-03
## DNA:DNA aln = evalue 1e-10
my $bit_score = 50;

my %relations;
my $blast_file_parsed = $output_name.".BLAST_parsed.txt\n";
my $results_file_parsed = $output_name.".duplicate_relations.txt\n";
open (OUT, ">$blast_file_parsed");
open (BLAST, $blast_file);
while (<BLAST>) {
#
#CBG27754.1     CBG27899.1      100.00  285     0       0       1       285     1       285     0.0     582     285     285     N/A
#0				1				2		3		4		5		6		7		8		9		10		11		12		13		14
# 10: E-value
# 11: Bit-Score
# 12: query length
# 13: subject length

	#print $_."\n";
	chomp;
	
	my @line = split("\t", $_);
	next if $line[0] eq $line[1]; #discard autohits
	next if $line[2] < $similarity; #similarity
	next if $line[10] > $evalue; #evalue
	next if $line[11] < $bit_score; #bit_score	

	#print "##############\n";
	#print $_."\n";
	#print "······················\n";

	#my $perc_len = ($line[3]/$line[13])*100;	
	my $perc_len_sub = ($line[3]/$line[13])*100; ## subject alignment 
	my $perc_len_query = ($line[3]/$line[12])*100; ## subject alignment 
	if ($perc_len_sub > $aln) { ## len
		if ($perc_len_query > $aln) { ## len
	
		#if ($perc_len > $aln) { ## Alignment length
		print OUT $_."\n";
		my $flag=0;
		foreach my $keys (sort keys %relations) {
			if ($relations{$keys}{$line[0]}) {$flag++; last;}
	        	if ($relations{$keys}{$line[1]}) {$flag++; last;}
        	}
		if ($flag == 0) {  
			my @array = ($line[2], $line[3], $line[12], $line[13]); 
			#	     similarity aln-length qlen       slen
			#push (@ { $relations{$line[0]}{$line[1]} }, @array);
			$relations{$line[0]}{$line[1]}++;
}}}}
close (BLAST); close (OUT);

## if no further processing
if ($stop_signal){ exit(); }

##print Dumper \%relations;
my $file2grep = "tmp.txt";
open (F, ">$file2grep");
foreach my $keys (keys %relations) {
	my @arr = keys %{ $relations{$keys} };
	print F $keys."\t".join(",", @arr)."\n";
}
close (F);

my @array = keys %relations;
my %better_relations;
for (my $i=0; $i < scalar @array; $i++) {
	##print "Searching $array[$i]\n";
	system("grep $array[$i] $file2grep > tmp2.txt");
	my @sub_subkeys = keys %{ $relations{$array[$i]} };
	
	for (my $j=0; $j < scalar @sub_subkeys; $j++) {
		system("grep $sub_subkeys[$j] $file2grep >> tmp2.txt");
	}
	system("cat tmp2.txt | sort | uniq > tmp3.txt");
	my $initial; my $lines=0;
	open (IN, "tmp3.txt");
	while(<IN>) {
		chomp;
		$lines++;
		my @line = split("\t", $_);
		if (!$initial) {$initial = $line[0];}
		my @ids = split(",", $line[1]);
		for (my $h=0; $h < @ids; $h++) {
			if ($initial eq $ids[$h]) {next;}
			$better_relations{$initial}{$ids[$h]}++;
		}
		if ($lines > 1) {
			$better_relations{$initial}{$line[0]}++;
	}}
	close(IN);
}

## remove temp file
system("rm *tmp*");

#print "######### Initial ###########\n";
#print Dumper \%relations;
#print "######### Final ###########\n";
#print Dumper \%better_relations;

open (RES, ">$results_file_parsed");
foreach my $keys (keys %better_relations) {
	my @subkey;
	foreach my $subkeys (keys %{ $better_relations{$keys} }) {
		push(@subkey,$subkeys);
	}
	print RES $keys."\t".join(",", @subkey)."\n";
}
close(RES);





