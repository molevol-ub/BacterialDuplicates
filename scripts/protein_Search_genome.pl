#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use Cwd 'abs_path';
use FindBin '$Bin';

my ($fasta, $cpus, $sim, $len, $evalue, $help, $blast_path, $name, @strain);

GetOptions(
	"proteins=s" => \$fasta,
	"sim=i" => \$sim,
	"len=i" => \$len,
	"strain=s" => \@strain,
	"evalue=i" => \$evalue,
	"name=s" => \$name,
	"CPU=i" => \$cpus,
	"h|help" => \$help,
	"BLAST_path=s" => \$blast_path
);

## Controlling options are provided
if (!$fasta) { &print_help; exit(); }
if (!$blast_path) { &print_help; exit(); }
if (!$name) { &print_help; exit(); }
if (!@strain) {&print_help; exit(); }

if (!$sim) {$sim=85;}
if (!$len) {$len=85;}
if (!$evalue) {$evalue=1e-10;}
if (!$cpus) {$cpus=2}

# Lets start
# 1. Generated folder for files
# 2. makeblastdb of proteins
# 3. blastp proteins vs itself
# 4. parse results
# 5. generate plot

## get path for scripts
my $script_paths = $Bin;
my $parse_blast_script = $script_paths."/parse_BLAST.pl";

print "############################\n";
print "Step 0: Check proteins\n";
print "############################\n";
my @proteins_file = &readFASTA_hash($fasta);
my %proteins_ids = %{ $proteins_file[0] }; 
my @order_proteins = @{ $proteins_file[1] }; 

print "############################\n";
print "Step 1: Make folder\n";
print "############################\n";
mkdir $name, 0755;
chdir $name;

my %files;
for (my $i=0; $i < scalar @strain; $i++) {	
	my @strain_Provided = split(",",$strain[$i]);
	my $strain_path = $strain_Provided[0];
	my $strain_name = $strain_Provided[1];
	
	print "####################################\n";
	print "Step 2: Makeblastdb $strain_name\n";
	print "####################################\n";
	
	my $db_name = $strain_name."_DB";
	my $blastdb_command = $blast_path."/makeblastdb -in $strain_path -input_type fasta -dbtype prot -out $db_name\n";
	
	print "Checking if already exists...\n";
	if (-f -r -s $db_name.".phr") { 
		print "Database $db_name exists in the folder provided. Do not generate it again...\n";
	} else {
		print "System call: $blastdb_command\n";
		system($blastdb_command);
	}
	
	print "#########################################\n";
	print "Step 3: BLAST proteins: $strain_name\n";
	print "#########################################\n";
	my $output_name = $db_name."_BLAST.out";
	my $blastp_command = $blast_path."/blastp -query $fasta -db $db_name -outfmt '6 std qlen slen' -num_threads $cpus -out $output_name";
	print "System call: $blastp_command\n";
	system($blastp_command);
	
	## parse_results
	print "############################\n";
	print "Step 4: Parse duplicate results\n";
	print "############################\n";
	my $out_parsed = $output_name;
	my $parse_command = "perl ".$parse_blast_script." $output_name $out_parsed $sim $len";
	print "System call: $parse_command\n";
	system($parse_command);

	##	
	$files{$strain_name} = $output_name.".duplicate_relations.txt"
}

## read duplicate relations & generate table
my $relations = "relations_proteins.tsv";
open (REL,">$relations");
print REL "#taxa\tQuery_Protein\tSubject_Proteins\n";

my %proteins;
foreach my $strain_name (keys %files) {
	my $file2read = $files{$strain_name};	
	open(RESULTS, "<$file2read");
	while (<RESULTS>) {
		chomp;
		my @items = split("\t", $_);
		my @subitems = split(",", $items[1]);
		$proteins{$strain_name}{$items[0]} = scalar @subitems;
		print REL $strain_name."\t".$items[0]."\t".join(",",@subitems)."\n";
	}
	close (RESULTS);
}

close (REL);

sleep 3;
print Dumper \%proteins;
sleep 3;

my $table_count = "table.csv";
open (TAB,">$table_count");
print TAB "strain,".join(",",@order_proteins)."\n";
foreach my $taxa (keys %proteins) {
	my $string = "$taxa";
	for (my $i=0; $i < scalar @order_proteins; $i++) {
		$string .= ",".$proteins{$taxa}{$order_proteins[$i]};
		print $string."\n";
		sleep 0.5;
	}
	print TAB $string."\n";
}
close (TAB);


sub print_help {
	print "############################\n";
	print "\tHELP Message:\n";
	print "############################\n";
	print "\nThis script generates a blast database for each strain provided and search the given proteins.\n";
	print "Usage:\nperl $0 -proteins proteins.fasta -name example -BLAST_path /path/to/BLAST/bin/\n-strain /mydirectory/DATA/ecoli/strain_BL21_proteins.faa,BL21 -strain /mydirectory/DATA/ecoli/strain_042_proteins.faa,Ecoli042 [-CPU nCPUs -sim 85 -len 85]\n\n";

	print "\nMandatory parameters:\n";
	print "############################\n";
	print "proteins: proteins in fasta format\n";
	print "name: name to add to identify files\n";
	print "BLAST_path: binary path contain blastp and makeblastdb\n";
	print "strain: several strains can be provided. Provide a comma separated value with path to protein file and id for the strain\n";
	print "Example: -strain /mydirectory/DATA/ecoli/strain_BL21_proteins.faa,BL21 -strain /mydirectory/DATA/ecoli/strain_042_proteins.faa,Ecoli042\n\n\n";
	
	print "Default parameters [in brakets]:\n";
	print "############################\n";
	print "CPU: 2\n";
	print "sim: 85 (BLAST similarity threshold)\n";
	print "len: 85 (BLAST similarity threshold)\n\n";
	
	print "Tips:\n";
	print "############################\n";
	print "- Provide absolute path for files\n";
	print "- Names can not start with numbers\n";
	print "- Names can not contain spaces neither commas\n\n\n";
}


sub readFASTA_hash {

	my $file = $_[0];
	my %hash; my @order;
	my $counter=0;
	open(FILE, $file) || die "Could not open the file $file [DOMINO.pm: readFASTA_hash]\n";
	$/ = ">"; ## Telling perl where a new line starts
	while (<FILE>) {		
		next if /^#/ || /^\s*$/;
		chomp;
    	my ($titleline, $sequence) = split(/\n/,$_,2);
    	next unless ($sequence && $titleline);
    	chop $sequence;
    	my @title=split(" ",$titleline);
    	$counter++;
    	$hash{$title[0]} = $titleline;
    	push (@order, $title[0]);
	}
	close(FILE); $/ = "\n";
	my $hashRef = \%hash;
	my $order_ref = \@order;
	return ($hashRef, $order_ref);
}
