This is an example workflow and additional details for the bioinformatic analysis generated for the analysis of duplicated genes within the genomes of gram positive cocci: *Staphylococcus aureus* and *Enterococcus faecalis*/*faecium*

## Citation

xx

## Usage

In order to provide an example of the procedure we followed we have generated an example set in the example folder.

#### 1) Retrieve data

Data would be downloaded from NCBI. Necessary data is basically translated cds protein for each strain, genomic fasta files and genomic gene feature files (gff). 

```
ATTENTION: DO NOT use file *protein.faa
   
Translated CDS contain directly translated coding sequence regions and sometimes proteins that are identically the same are collapsed into 1 entry into database, so a duplicated gene that has two different positions in the genome, two translated cds would only have one protein.
```

If strains are deposited on GenBank, data can be downloaded using a script we provide here: [NCBI_downloader.pl](https://github.com/molevol-ub/BacterialDuplicates/blob/master/scripts/perl/NCBI_downloader.pl).

```
Usage: 
perl BacterialDuplicates/scripts/NCBI_downloader.pl csv_file option
    csv_file: csv file containing ftp site and name
    option: gff,protein,feature,CDS,genome,ALL
```
Please provide a csv file and option to download data. Use ALL for the download of all information available for each strain specified. 

Please select from https://www.ncbi.nlm.nih.gov/genome your strains of interest and generate a comma separated (csv) table containing the ftp site for each strain and the strain name you would like to add. Please do not use any spaces or special characters.

Example data file ([example_strains2download.csv](https://github.com/molevol-ub/BacterialDuplicates/blob/master/Gram_positive/example/example_strains2download.csv)): 
```
## Saureus
ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/013/425/GCF_000013425.1_ASM1342v1,GCA_000013425.1_Saureus_NCTC8325

## Efaecalis
ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/391/485/GCF_000391485.2_ASM39148v2,GCA_000391485.2_B594

## Efaecium
ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/720/945/GCF_001720945.1_ASM172094v1,GCA_001720945.1_ISMMS_VRE_1

## Sepidermidis
ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/006/094/375/GCF_006094375.1_ASM609437v1,GCA_006094375.1_ATCC_14990

## Scarnosus
ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/405/GCF_000009405.1_ASM940v1,GCA_000009405.1_TM300

## Sxylosus
ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/709/415/GCF_000709415.1_ASM70941v1,GCA_000709415.1_SMQ-121
```

See an example of this data retrieval process in the previous [Ecoli analysis example](https://github.com/molevol-ub/BacterialDuplicates/blob/master/Ecoli/Ecoli_genome.md#example)

We will retrieve all information for each strain but we basically rely on the header information supplied by translated_cds.faa file

```
>lcl|FN554766.1_prot_CBG32835.1_1 [gene=thrA] [locus_tag=EC042_0001] [db_xref=GOA:D3H385,InterPro:IPR001048,InterPro:IPR001341,InterPro:IPR001342,InterPro:IPR002912,InterPro:IPR005106,InterPro:IPR011147,InterPro:IPR016040,InterPro:IPR018042,InterPro:IPR019811,InterPro:IPR027795,UniProtKB/TrEMBL:D3H385] [protein=bifunctional aspartokinase I/homoserine dehydrogenase I] [protein_id=CBG32835.1] [location=336..2798] [gbkey=CDS]
MRVLKFGGTSVANAERFLRVADILESNARQGQVATVLSAPAKITNHLVAMIEKTISGQDALPNISDAERIFAELLTGLAA
....
```

NOTE: It is possible to generate the same results from a gff file, but it is not implemented neither intended right now.
Please contact us for further explanation or clarification or to show interest as we might by thinking of implementing this feature. 
