#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

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

my %proteins;
my $relations = "relations_proteins.csv";
open (REL,">$relations");
print REL "#taxa\tQuery_Protein\tSubject_Proteins\n";
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
	my $blastp_command = $blast_path."/blastp -query $fasta -db $db_name -outfmt '6 std qlen slen staxids' -num_threads $cpus -out $output_name";
	print "System call: $blastp_command\n";
	system($blastp_command);
	
	## parse_results
	my %relations;
	my $blast_file_parsed = $output_name."_parsed.txt";
	open (OUT, ">$blast_file_parsed");
	open (BLAST, $output_name);
	while (<BLAST>) {
		chomp;
		my @line = split("\t", $_);
		next if $line[0] eq $line[1]; #discard autohits
		next if $line[2] < $sim; ## Similarity > 85%
		next if $line[10] > $evalue;

		my $perc_len_sub = ($line[3]/$line[13])*100; ## subject alignment 
		my $perc_len_query = ($line[3]/$line[12])*100; ## subject alignment 
		if ($perc_len_sub > $len) { ## len
			if ($perc_len_query > $len) { ## len
				print OUT $_."\n";
				my $flag=0;
				$relations{$line[0]}{$line[1]}++;
	}}}
	close (BLAST); close (OUT);
	#print Dumper \%relations;
	
	## Initialize hash
	foreach my $ids (keys %proteins_ids) { $proteins{$strain_name}{$ids}=0; }
	foreach my $keys (keys %relations) {
		my @ids = keys %{ $relations{$keys} };
		$proteins{$strain_name}{$keys} = scalar @ids;
		print REL $strain_name."\t".$keys."\t".join(",",@ids)."\n";
}}
close (REL);

#print Dumper \%proteins;

my $table_count = "table.csv";
open (TAB,">$table_count");
print TAB "strain,".join(",",@order_proteins)."\n";
foreach my $taxa (keys %proteins) {
	my $string = "$taxa";
	for (my $i=0; $i < scalar @order_proteins; $i++) {
		$string .= ",".$proteins{$taxa}{$order_proteins[$i]};
	}
	print TAB $string."\n";
}
close (TAB);


sub print_help {
	print "############################\n";
	print "\tHELP Message:\n";
	print "############################\n";
	print "\nThis script generates a blast database for each strain provided and search the given proteins.\n";
	print "Usage:\nperl $0 -proteins proteins.fasta -name example -BLAST_path /path/to/BLAST/bin/\n-strain /mydirectory/DATA/ecoli/strain_BL21_proteins.faa,BL21 -strain /mydirectory/DATA/ecoli/strain_042_proteins.faa,Ecoli042 [-n CPUs -sim 85 -len 85]\n\n";

	print "\nMandatory parameters:\n";
	print "############################\n";
	print "proteins: proteins in fasta format\n";
	print "name: name to add to identify files\n";
	print "BLAST_path: binary path contain blastp and makeblastdb\n";
	print "strain: several strains can be provided. Provide a comma separated value with path to protein file and id for the strain\n";
	print "Example: -strain /mydirectory/DATA/ecoli/strain_BL21_proteins.faa,BL21 -strain /mydirectory/DATA/ecoli/strain_042_proteins.faa,Ecoli042\n\n\n";
	
	print "Default parameters [in brakets]:\n";
	print "############################\n";
	print "CPUs: 2\n";
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