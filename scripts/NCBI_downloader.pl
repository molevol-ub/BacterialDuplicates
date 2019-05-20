#!/usr/bin/perl
use strict;
use warnings;

##############################################################################
##	Jose Fco. Sanchez Herrero, 25/09/2018 jfsanchezherrero@ub.edu			##
##############################################################################

if (!@ARGV) {
	print "Usage: perl $0 csv_file option\n";
	print "\t- csv_file: containing ftp site and name\n";
	print "\t- option: gff,protein,feature,CDS,genome,ALL\n";
	exit();
}

my $ftp_file = $ARGV[0];
my $id = $ARGV[1];

my %strains;
open (FILE, "<$ftp_file");
while (<FILE>) {
	chomp;
	next if ($_ =~ /#.*/);
	next if ($_ =~ /^$/);
	my @split=split(",",$_);
	print $_."\n";
	$strains{$split[1]} = $split[0];
}
close(FILE);

foreach my $keys (keys %strains) {
	my @ftp_folder_path = split("/", $strains{$keys});
	my $name = $ftp_folder_path[-1];
	my $file = $strains{$keys}."/".$name;
	
	if (-d $keys) {
		print "DIR $keys exists...\n";
	} else {
		mkdir $keys,0755;
	}	
	chdir $keys;
	
	my @files;
	my $temp_GFF = $file."*gff.gz";
	my $temp_PRO = $file."*protein.faa.gz";
	my $temp_FEA = $file."*feature_table.txt.gz";
	my $temp_CDS = $file."*cds*";
	my $temp_geno = $file."*genomic.fna.gz";
	
	if ($id eq "gff") {
		push (@files, $temp_GFF);
	} elsif ($id eq "protein") {
		push (@files, $temp_PRO);
	} elsif ($id eq "feature") {
		push (@files, $temp_FEA);
	} elsif ($id eq "CDS") {
		push (@files, $temp_CDS);
        } elsif ($id eq "genome") {
                push (@files, $temp_geno);
	} elsif ($id eq "ALL") {
		push (@files, $temp_GFF);
		push (@files, $temp_PRO);
		push (@files, $temp_FEA);
		push (@files, $temp_CDS);
		push (@files, $temp_geno);
	}
	
	print "+ Downloading now:\n\n";
	print "\n\n+ Download file(s) for: $name\n\n";
	
	for (my $f=0; $f < scalar @files;$f++) {
		print "wget --passive-ftp $files[$f]\n";
		system("wget --passive-ftp $files[$f]");
	
		# Gunzip file
		my @array_tmp2 = split("/", $files[$f]);
		my $file_downloaded = $array_tmp2[-1];
		system("gunzip -f $file_downloaded");
	}
	chdir "..";
}
print "+ Finish Donwloading of files\n";



