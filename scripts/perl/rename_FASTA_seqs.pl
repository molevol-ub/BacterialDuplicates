#!/usr/bin/perl
use strict;
use warnings;

my $fasta = $ARGV[0];
my $file_name = $ARGV[1];
my $name = $ARGV[2];
my $option = $ARGV[3];

if (scalar @ARGV != 4) {
	print "\n########\n";
	print "Usage:\n########\n\n";
	print "Please provide the next arguments:\n";
	print "perl ".$0." fasta_file name_file name2add ADD|REPLACE|BEGIN|ADD_BEGIN\n\n";
	print "ADD: will add the id plus a counter\nREPLACE: will discard the previous name and add this unique id and counter\nBEGIN: will keep the first split of the id and add at the beginning the given name\n\n";
	exit();
} 
if ($option eq "ADD") {
} elsif ($option eq "REPLACE") {
} elsif ($option eq "BEGIN") {
} elsif ($option eq "ADD_BEGIN") {
} else { exit(); }
my $counter_seqs = 0;

open(FILE, $fasta) || die "Could not open the $fasta ...\n";
open (OUT, ">$file_name");
$/ = ">"; ## Telling perl where a new line starts
while (<FILE>) {		
	next if /^#/ || /^\s*$/;
	chomp;
	my ($titleline, $sequence) = split(/\n/,$_,2);
	next unless ($sequence && $titleline);
    	$sequence =~ s/\n//g;
	#$titleline =~ s/ /_/g;
	$counter_seqs++;
	
	if ($option eq "ADD") {
	        $titleline =~ s/ /_/g;
		print OUT ">".$titleline."-".$name."_".$counter_seqs."\n".$sequence."\n";
	} elsif ($option eq "REPLACE") {
	        $titleline =~ s/ /_/g;
		print OUT ">".$name."_".$counter_seqs."\n".$sequence."\n";
		print $titleline."\t".$name."_".$counter_seqs."\n";
        } elsif ($option eq "BEGIN") {
                my @array = split("-", $titleline);
               #if ($array[0] =~ /.*prot_(.*)/) {$protein_name=$1;} else {$protein_name=$array[0];}
#               print $protein_name."\n";
                print OUT ">$array[0]\n$sequence\n";
                #print $name."_".$protein_name."\t".$titleline."\n";
	} elsif ($option eq "ADD_BEGIN") {
        	my @array = split(" ", $titleline);
		my $protein_name;
		if ($array[0] =~ /.*prot_(.*)/) {$protein_name=$1;} else {$protein_name=$array[0];}
#		print $protein_name."\n";
	        print OUT ">".$name."_".$protein_name."\n".$sequence."\n";
                print $name."_".$protein_name."\t".$titleline."\n";

	}
}
close (FILE); close (OUT);
