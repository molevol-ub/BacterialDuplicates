# BacterialDuplicates git repository

## Description

This scripts and data collected in this repo correspond to the analysis performed to identify bacterial duplicated genes that led to the following paper. 

https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-019-5683-4

## Citation

Gene duplications in the E. coli genome: common themes among pathotypes. Manuel Bernabeu, José Francisco Sánchez-Herrero, Pol Huedo, Alejandro Prieto, Mário Hüttener, Julio Rozas and Antonio Juárez. BMC Genomics 2019 20:313, https://doi.org/10.1186/s12864-019-5683-4


## Usage

In order to provide an example of the procedure we followed we have generated an example set in the example folder.

1) Retrieve data

Data can be downloaded from NCBI or could be provided by the user. Necessary data is basically protein annotated genes, gff file and genome fasta file for each strain. 

If strains are deposited on GenBank, data can be downloaded using a script we provide here: [NCBI_downloader.pl](https://github.com/molevol-ub/BacterialDuplicates/blob/master/scripts/NCBI_downloader.pl).

```
Usage: 
perl BacterialDuplicates/scripts/NCBI_downloader.pl csv_file option
    csv_file: csv file containing ftp site and name
    option: gff,protein,feature,CDS,genome,ALL
```
Please provide a csv file and option to download data. Use ALL for the download of all information available for each strain specified. 

Please select from https://www.ncbi.nlm.nih.gov/genome your strains of interest and generate a comma separated (csv) table containing the ftp site for each strain and the strain name you would like to add. Please do not use any spaces or special characters.

Example data file ([example_strains2download.csv](https://github.com/molevol-ub/BacterialDuplicates/blob/master/example/example_strains2download.csv)): 

```
## Commensal
ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/482/265/GCF_000482265.1_EC_K12_MG1655_Broad_SNP,GCF_000482265.1_K-12_MG1655
## EAEC
ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/027/125/GCA_000027125.1_ASM2712v1,GCA_000027125.1_Ecoli042
```


To reproduce the example provided here:

`perl BacterialDuplicates/scripts/NCBI_downloader.pl BacterialDuplicates/example/strains2download.csv ALL`


Two folders should be generated named as the strain provided containing several files: GFF, protein and nucleotide sequences of the protein-coding genes, genomic sequence and other information.


2) Gene duplication among strains analysis

The script [duplicate_search_bacteria.pl](https://github.com/molevol-ub/BacterialDuplicates/blob/master/scripts/duplicate_search_bacteria.pl) generates a blast database of the provide protein fasta and searches for putative duplicates.

Please provide absolute path for files and remember that names can not contain spaces or special characters.

```
Usage:
perl BacterialDuplicates/scripts/duplicate_search_bacteria.pl 
    -fasta proteins.fasta 
    -script_path /path/to/script/parse_BLAST.pl 
    -name example 
    -BLAST_path /path/to/BLAST/bin 
    [-n CPUs -sim 85 -len 85]
```

##### Mandatory parameters:
fasta: proteins in fasta format E.g. *.protein.faa

name: name to add to identify files

BLAST_path: binary path contain blastp and makeblastdb. E.g. /usr/bin/, /software/ncbi-blast/bin, etc.

script_path: path for [parse_BLAST.pl](https://github.com/molevol-ub/BacterialDuplicates/blob/master/scripts/parse_BLAST.pl)

##### Default parameters [in brakets]:
CPUs: 2

sim: 85

len: 85



3) Gene duplication between strains analysis
