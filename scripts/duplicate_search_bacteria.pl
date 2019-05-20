#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
##############################################################################
##	Jose Fco. Sanchez Herrero, 25/09/2018 jfsanchezherrero@ub.edu			##
##############################################################################

my ($fasta, $cpus, $sim, $len, $help, $blast_path, $name, $script_path);

GetOptions(
	"fasta=s" => \$fasta,
	"sim=i" => \$sim,
	"len=i" => \$len,
	"name=s" => \$name,
	"CPU=i" => \$cpus,
	"h|help" => \$help,
	"script_path=s" => \$script_path,
	"BLAST_path=s" => \$blast_path
);

## Controlling options are provided
if (!$fasta) { &print_help; exit(); }
if (!$blast_path) { &print_help; exit(); }
if (!$name) { &print_help; exit(); }
if (!$script_path) { &print_help; exit(); }

if (!$len) {$len=85;}
if (!$sim) {$sim=85;}
if (!$cpus) {$cpus=2;}

# Lets start
# 1. Generated folder for files
# 2. makeblastdb of proteins
# 3. blastp proteins vs itself
# 4. parse results
# 5. generate plot

print "############################\n";
print "Step 1: Make folder\n";
print "############################\n";
mkdir $name, 0755;
chdir $name;

print "############################\n";
print "Step 2: Discard pseudogenes\n";
print "############################\n";

my $fasta1 = $fasta."_tmp1";
my $fasta2 = $fasta."_tmp";
system("sed 's/.>/./g' $fasta > $fasta1");
system("sed 's/=</=/g' $fasta1 > $fasta2");
open(FILE, $fasta2) || die "Could not open the file $fasta\n";
my $clean_fasta = $fasta."_clean.fasta";
open (CLEAN, ">$clean_fasta");
$/ = ">"; ## Telling perl where a new line starts
while (<FILE>) {		
	next if /^#/ || /^\s*$/;
	chomp;
   	my ($titleline, $sequence) = split(/\n/,$_,2);
   	next unless ($sequence && $titleline);
   	chop $sequence;
	if ($titleline =~ /.*pseudo=true.*/) { next; }
	$titleline =~ s/\-/\_/g;
	$titleline =~ s/\|/\_/g;
	print CLEAN ">".$titleline."\n".$sequence."\n";
}
close(FILE); $/ = "\n";
close(CLEAN);

print "############################\n";
print "Step 3: Makeblastdb\n";
print "############################\n";
my $db_name = $name."_DB";
my $blastdb_command = $blast_path."/makeblastdb -in $clean_fasta -input_type fasta -dbtype nucl -out $db_name\n";
print "System call: $blastdb_command\n";
system($blastdb_command);

print "############################\n";
print "Step 4: BLAST proteins\n";
print "############################\n";
my $output_name = $name."_BLAST.out";
my $blastp_command = $blast_path."/blastn -query $clean_fasta -db $db_name -outfmt '6 std qlen slen' -num_threads $cpus -out $output_name";
print "System call: $blastp_command\n";
system($blastp_command);

print "############################\n";
print "Step 5: Parse results\n";
print "############################\n";
my $out_parsed = $name."_parsed.txt";
my $parse_command = "perl ".$script_path." $output_name	$clean_fasta $out_parsed $sim $len";
print "System call: $parse_command\n";
system($parse_command);

print "############################\n";
print "Step 6: Generate Plot\n";
print "############################\n";
print "Open RStudio and run script with the output coordinates generated\n";

sub print_help {
	print "############################\n";
	print "\tHELP Message:\n";
	print "############################\n";
	print "\nThis script generates a blast database of the provide protein fasta and searches for putative duplicates.\n";
	print "Usage:\nperl $0 -fasta proteins.fasta -script_path /path/to/script/parse_BLAST.pl -name example -BLAST_path /path/to/BLAST/bin [-n CPUs -sim 85 -len 85]\n\n";

	print "\nMandatory parameters:\n";
	print "############################\n";
	print "fasta: proteins in fasta format\n";
	print "name: name to add to identify files\n";
	print "BLAST_path: binary path contain blastp and makeblastdb\n\n\n";
	print "script_path: path for parse_BLAST.pl\n";
	print "Default parameters [in brakets]:\n";
	print "############################\n";
	print "CPUs: 2\n";
	print "sim: 85 (BLAST similarity threshold)\n";
	print "len: 85 (BLAST similarity threshold)\n\n";
	
	print "Tips:\n";
	print "############################\n";
	print "- Provide absolute path for files\n";
	print "- Names can not contain spaces\n\n\n";
}

