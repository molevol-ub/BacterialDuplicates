#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
##############################################################################
##	Jose Fco. Sanchez Herrero, 25/09/2018 jfsanchezherrero@ub.edu			##
##############################################################################

my $blast_file = $ARGV[0];
my $output_name = $ARGV[1];
my $similarity = $ARGV[2];
my $aln = $ARGV[3];
#my $evalue=$ARGV[4];

if (!@ARGV) {print "Usage:\nperl $0 blast_results output [aln_len similarity]\n";exit();}

## default
if (!$similarity) {$similarity=80;}
if (!$aln) {$aln=80;}
#if (!$evalue) {$evalue=1e-05;}

my %relations;
my $blast_file_parsed = $output_name.".BLAST_parsed.txt\n";
my $results_file_parsed = $output_name.".duplicate_relations.txt\n";
open (OUT, ">$blast_file_parsed");
open (BLAST, $blast_file);
while (<BLAST>) {
#
#CBG27754.1     CBG27899.1      100.00  285     0       0       1       285     1       285     0.0     582     285     285     N/A
#0				1				2		3		4		5		6		7		8		9		10		11		12		13		14
# 12: query length
# 13: subject length

	chomp;
	my @line = split("\t", $_);
	next if $line[0] eq $line[1]; #discard autohits
	next if $line[2] < $similarity; #similarity
	#next if $evalue > $line[10]; #evalue
	
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
#exit();

open (RES, ">$results_file_parsed");
foreach my $keys (keys %better_relations) {
	my @subkey;
	foreach my $subkeys (keys %{ $better_relations{$keys} }) {
		push(@subkey,$subkeys);
	}
	print RES $keys."\t".join(",", @subkey)."\n";
}
close(RES);





