#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
##############################################################################
##	Jose Fco. Sanchez Herrero, 03/02/2020 jfsanchezherrero@ub.edu			##
##############################################################################

my $fasta = $ARGV[0];
if (!@ARGV) {print "Usage:\nperl $0 fasta_file\n";exit();}

my $fasta1 = $fasta."_tmp1";
my $fasta2 = $fasta."_tmp";
my $fasta3 = $fasta."_tmp_1";
system("sed 's/.>/./g' $fasta > $fasta1");
system("sed 's/<././g' $fasta1 > $fasta3");
system("sed 's/=</=/g' $fasta3 > $fasta2");

open(FILE, $fasta2) || die "Could not open the file $fasta\n";
$/ = ">"; ## Telling perl where a new line starts
while (<FILE>) {		
	next if /^#/ || /^\s*$/;
	chomp;
   	my ($titleline, $sequence) = split(/\n/,$_,2);
   	next unless ($sequence && $titleline);
   	chop $sequence;
	#if ($titleline =~ /.*pseudo=true.*/) { next; }
	$titleline =~ s/\-/\_/g;
	$titleline =~ s/\|/\_/g;
	print ">".$titleline."\n".$sequence."\n";
}
close(FILE); $/ = "\n";

## remove temp file
system("rm *_tmp*");
