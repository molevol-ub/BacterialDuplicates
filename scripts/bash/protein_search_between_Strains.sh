
## initiate command file
echo "perl ./BacterialDuplicates/scripts/perl/protein_Search_genome.pl -BLAST_path /path/to/ncbi-blast-2.4.0/bin/ -CPU X -name name_desired " > command.sh

## create fasta file with selected proteins
echo "-proteins selected_proteins.fasta " >> command.sh

## retrieve all strains clean fasta and ID
for i in `dir results`;
do
	file=`readlink -f results/$i/*clean.fasta`
        echo "-strain "$file,$i;
done > strains.txt

## add to command
cat strains.txt | tr '\n' ' ' >> command.sh

## execute command.sh
sh command.sh
