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

If strains are deposited on GenBank, data can be downloaded using a script we provide here: NCBI_downloader.pl. Please provide a file and option. Use ALL for the download of all information available for each strain specified. 

```
Usage: 
perl BacterialDuplicates/scripts/NCBI_downloader.pl ftp_folder option
    file: csv file containing ftp site and name
    option: gff,protein,feature,CDS,genome,ALL
```

Please select from https://www.ncbi.nlm.nih.gov/genome your strains of interest and generate a comma separated (csv) table containing the ftp site for each strain and the strain name you would like to add. Please do not use any spaces or special characters.

To reproduce the example provided here:

`perl BacterialDuplicates/scripts/NCBI_downloader.pl BacterialDuplicates/example/strains2download.csv ALL`

2) Analysis
