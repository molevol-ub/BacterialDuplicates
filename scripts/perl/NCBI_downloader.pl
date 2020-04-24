#!/usr/bin/perl
use strict;
use warnings;

##############################################################################
##	Jose Fco. Sanchez Herrero, 25/09/2018 jfsanchezherrero@ub.edu			##
##############################################################################

if (!@ARGV) {
	print "Usage: perl $0 csv_file option\n";
	print "\t- csv_file: containing ftp site and name\n";
	print "\t- option: gff,gbf,protein,feature,CDS,genome,ALL\n";
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
		next;
	} else {
		mkdir $keys,0755;
	}	
	chdir $keys;
	
	my %files;
	my $temp_GFF = $file."*gff*";
	my $temp_PRO = $file."*protein.faa*";
	my $temp_FEA = $file."*feature_table.txt*";
	my $temp_CDS = $file."*cds*";
	my $temp_geno = $file."*genomic.fna*";
	my $temp_gbf = $file."*genomic.gbff*";
	
	if ($id eq "gff") {
		$files{"gff"} = $temp_GFF;
	} elsif ($id eq "protein") {
		$files{"protein.faa"} = $temp_PRO;
	} elsif ($id eq "feature") {
		$files{"feature_table.txt"} = $temp_FEA;
	} elsif ($id eq "CDS") {
		$files{"cds"} = $temp_CDS;
        } elsif ($id eq "genome") {
                $files{"genomic.fna"} = $temp_geno;
        } elsif ($id eq "gbf") {
                $files{"genomic.gbf"} = $temp_gbf;
	} elsif ($id eq "ALL") {
                $files{"gff"} = $temp_GFF;
                $files{"protein.faa"} = $temp_PRO;
                $files{"feature_table.txt"} = $temp_FEA;
                $files{"cds"} = $temp_CDS;
                $files{"genomic.fna"} = $temp_geno;
		$files{"genomic.gbff"} = $temp_gbf;
	}
	
	print "+ Downloading now:\n\n";
	print "\n\n+ Download file(s) for: $name\n\n";
	
	foreach my $keys2 (keys %files) {
		
		print "File $keys2: $files{$keys2}\n";
		print "wget --passive-ftp $files{$keys2}\n";
		system("wget --passive-ftp $files{$keys2}");
	
		# Gunzip file
		my @array_tmp2 = split("/", $files{$keys2});
		my $file_downloaded = $array_tmp2[-1];
		system("gunzip -f $file_downloaded");
	}
	chdir "..";
}
print "+ Finish downloading of files\n";


sub filesExist { 
	return scalar ( my @x = `ls -1a "$_[0]" 2> /dev/null` ) 
}
