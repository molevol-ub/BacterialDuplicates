#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Cwd 'abs_path';
use FindBin '$Bin';

##############################################################################
##	Jose Fco. Sanchez Herrero, 25/09/2018 jfsanchezherrero@ub.edu	    ##
##	Updated: January, 2020						    ##
##############################################################################

my ($fasta, $cpus, $sim, $len, $help, $blast_path, $name, $CARD_blast_db, $vfdb_blast_db);

GetOptions(
	"fasta=s" => \$fasta,
	"sim=i" => \$sim,
	"len=i" => \$len,
	"name=s" => \$name,
	"CPU=i" => \$cpus,
	"h|help" => \$help,
	"BLAST_path=s" => \$blast_path,
	"CARD_db=s" => \$CARD_blast_db,
	"vfdb_db=s" => \$vfdb_blast_db
);

## get path for scripts
my $script_paths = $Bin;

## Controlling options are provided
if (!$fasta) { &print_help; exit(); }
if (!$blast_path) { &print_help; exit(); }
if (!$name) { &print_help; exit(); }
if (!$CARD_blast_db) { &print_help; exit(); }
if (!$vfdb_blast_db) { &print_help; exit(); }

if (!$len) {$len=85;}
if (!$sim) {$sim=85;}
if (!$cpus) {$cpus=2;}

# Lets start
# 1: Make folder
# 2: Clean sequences
# 3: Makeblastdb
# 4: BLAST proteins
# 5: Parse BLAST
# 6: Generate duplicate results
# 7: Get fasta seqs results
# 8: BLAST duplicated proteins vs. CARD
# 9: Parse BLAST duplicated CARD results
# 10: BLAST ALL proteins vs. CARD
# 11: Parse BLAST ALL CARD results
# 12: BLAST proteins vs. VFDB
# 13: Parse BLAST VFDB results
# 14: BLAST ALL proteins vs. VFDB
# 15: Parse BLAST ALL VFDB results
# 16: Generate Plot

print "############################\n";
print "Step 1: Make folder\n";
print "############################\n";
mkdir $name, 0755;
chdir $name;

print "############################\n";
print "Step 2: Clean sequences\n";
print "############################\n";
my $clean_fasta_script = $script_paths."/clean_fasta.pl";
my $clean_fasta = $name."_clean.fasta";
my $clean_command = "perl ".$clean_fasta_script." $fasta > $clean_fasta";
print "System call: $clean_command\n";
system($clean_command);

print "############################\n";
print "Step 3: Makeblastdb\n";
print "############################\n";
my $db_name = $name."_DB";
my $blastdb_command = $blast_path."/makeblastdb -in $clean_fasta -input_type fasta -dbtype prot -out $db_name\n";
print "System call: $blastdb_command\n";
system($blastdb_command);

print "############################\n";
print "Step 4: BLAST proteins\n";
print "############################\n";
my $output_name = $name."_BLAST.out";
my $blastp_command = $blast_path."/blastp -query $clean_fasta -db $db_name -outfmt '6 std qlen slen' -num_threads $cpus -out $output_name";
print "System call: $blastp_command\n";
system($blastp_command);

print "############################\n";
print "Step 5: Parse BLAST\n";
print "############################\n";
my $parse_blast_script = $script_paths."/parse_BLAST.pl";
my $out_parsed = $name."_parsed.txt";
my $parse_command = "perl ".$parse_blast_script." $output_name	$out_parsed $sim $len";
print "System call: $parse_command\n";
system($parse_command);

print "############################\n";
print "Step 6: Generate duplicate results\n";
print "############################\n";
my $duplicate_results = $script_paths."/generate_results_duplicates.pl";
my $duplicates_relations = $name."_parsed.txt.duplicate_relations.txt";
my $results_command = "perl ".$duplicate_results." $duplicates_relations $clean_fasta $output_name";
print "System call: $results_command\n";
system($results_command);

print "############################\n";
print "Step 7: Get fasta seqs results\n";
print "############################\n";
my $get_ids_script = $script_paths."/get-seq_ids.pl";
my $out_duplicates = $output_name.".allseqs_duplicated.fasta";
my $ids_duplicated = abs_path()."/".$output_name.".allseqs_duplicated.ids.txt";
my $parse_command_duplicates= "perl ".$get_ids_script." $ids_duplicated $clean_fasta $out_duplicates";
print "System call: $parse_command_duplicates\n";
system($parse_command_duplicates);

print "############################\n";
print "Step 8: BLAST duplicated proteins vs. CARD\n";
print "############################\n";
my $output_name_CARD = $name."-BLAST_CARD_duplicates.out";
my $blastp_command_card = $blast_path."/blastp -query $out_duplicates -db $CARD_blast_db -outfmt '6 std qlen slen' -num_threads $cpus -out $output_name_CARD";
print "System call: $blastp_command_card\n";
system($blastp_command_card);

### relaxed similarity and length thresholds
my $sim_Relaxed = 30;
my $len_Relaxed = 85;

print "############################\n";
print "Step 9: Parse BLAST duplicated CARD results\n";
print "############################\n";
my $out_parsed_CARD = $output_name_CARD;
my $parse_command_CARD = "perl ".$parse_blast_script." $output_name_CARD $out_parsed_CARD $sim_Relaxed $len_Relaxed";
print "System call: $parse_command_CARD\n";
system($parse_command_CARD);
## It returns the number of duplicated groups of duplicates within CARD db

print "############################\n";
print "Step 10: BLAST ALL proteins vs. CARD\n";
print "############################\n";
my $output_name_CARD_all = $name."-BLAST_CARD_all.out";
my $blastp_command_card_all = $blast_path."/blastp -query $clean_fasta -db $CARD_blast_db -outfmt '6 std qlen slen' -num_threads $cpus -out $output_name_CARD_all";
print "System call: $blastp_command_card_all\n";
system($blastp_command_card_all);

print "############################\n";
print "Step 11: Parse BLAST ALL CARD results\n";
print "############################\n";
my $out_parsed_CARD_all = $output_name_CARD_all;
my $parse_command_CARD_all = "perl ".$parse_blast_script." $output_name_CARD_all $out_parsed_CARD_all $sim_Relaxed $len_Relaxed stop";
print "System call: $parse_command_CARD_all\n";
system($parse_command_CARD_all);

print "############################\n";
print "Step 12: BLAST proteins vs. VFDB\n";
print "############################\n";
my $output_name_VFDB = $name."-BLAST_VFDB_duplicates.out";
my $blastp_command_vfdb = $blast_path."/blastp -query $out_duplicates -db $vfdb_blast_db -outfmt '6 std qlen slen' -num_threads $cpus -out $output_name_VFDB";
print "System call: $blastp_command_vfdb\n";
system($blastp_command_vfdb);

print "############################\n";
print "Step 13: Parse BLAST VFDB results\n";
print "############################\n";
my $out_parsed_VFDB = $output_name_VFDB;
my $parse_command_VFDB = "perl ".$parse_blast_script." $output_name_VFDB $out_parsed_VFDB $sim_Relaxed $len_Relaxed";
print "System call: $parse_command_VFDB\n";
system($parse_command_VFDB);
## It returns the number of duplicated groups of duplicates within VFDB db

print "############################\n";
print "Step 14: BLAST ALL proteins vs. VFDB\n";
print "############################\n";
my $output_name_VFDB_all = $name."-BLAST_VFDB_all.out";
my $blastp_command_VFDB_all = $blast_path."/blastp -query $clean_fasta -db $vfdb_blast_db -outfmt '6 std qlen slen' -num_threads $cpus -out $output_name_VFDB_all";
print "System call: $blastp_command_VFDB_all\n";
system($blastp_command_VFDB_all);

print "############################\n";
print "Step 15: Parse BLAST ALL VFDB results\n";
print "############################\n";
my $out_parsed_VFDB_all = $output_name_VFDB_all;
my $parse_command_VFDB_all = "perl ".$parse_blast_script." $output_name_VFDB_all $out_parsed_VFDB_all $sim_Relaxed $len_Relaxed stop";
print "System call: $parse_command_VFDB_all\n";
system($parse_command_VFDB_all);

print "############################\n";
print "Step 16: Generate Plot\n";
print "############################\n";

sub print_help {
	print "############################\n";
	print "\tHELP Message:\n";
	print "############################\n";
	print "\nThis script generates a blast database of the provide protein fasta and searches for putative duplicates.\n";
	print "Usage:\nperl $0 -fasta proteins.fasta -name example 
		-BLAST_path /path/to/BLAST/bin 
		-CARD_db /path/to/CARD_databases/blast_id_name
		-vfdb_db /path/to/VFDB_databases/blast_id_name
		[-CPU nCPU -sim 85 -len 85]\n\n";

	print "\nMandatory parameters:\n";
	print "############################\n";
	print "fasta: proteins in fasta format\n";
	print "name: name to add to identify files\n";
	print "BLAST_path: binary path contain blastp and makeblastdb\n\n\n";
	print "CARD_db: path for CARD protein databases indexed by makeblast db\n";
	print "vfdb_db: path for VFDB protein databases indexed by makeblast db\n";
	
	print "Default parameters [in brakets]:\n";
	print "############################\n";
	print "CPU: 2\n";
	print "sim: 85 (BLAST similarity threshold)\n";
	print "len: 85 (BLAST similarity threshold)\n\n";
	
	print "Tips:\n";
	print "############################\n";
	print "- Provide absolute path for files\n";
	print "- Names can not contain spaces\n\n\n";
}

