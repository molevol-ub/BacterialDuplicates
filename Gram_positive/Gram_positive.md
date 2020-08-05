# Duplicates in Gram positive cocci

This is an example workflow and additional details for the bioinformatic analysis generated for the analysis 
of duplicated genes within the genomes of gram positive cocci: *Staphylococcus aureus* and *Enterococcus faecalis*/*faecium*

## Table of Contents

- [Usage](#usage)
  * [Retrieve data](#retrieve-data)
     * [Example](#example)
  * [Gene duplication analysis](#gene-duplication-among-strains-analysis)
     * [Explanation of the workflow](#explanation-of-the-workflow)
     * [Virulence and Resistance Databases retrieval](#virulence-and-resistance-databases-retrieval)
     * [Parameters](#parameters)
     * [Example of use](#example-of-use)
     * [Generate plot](#generate-plot)
  * [Gene duplication between strains analysis](#gene-duplication-between-strains-analysis)
     * [Parameters](#parameters)
     * [Example](#example)
- [Supplementary information](#supplementary-information)
- [Citation](#citation)

## Usage

In order to provide an example of the procedure we followed we have generated an example set in the example folder 
under the `Gram_positive` [folder](https://github.com/molevol-ub/BacterialDuplicates/tree/master/Gram_positive).

```
ATTENTION: To replicate the analysis and the example we provide here, we encourage to download the GitHub 
release corresponding [v2.0](https://github.com/molevol-ub/BacterialDuplicates/releases/tag/v2.0). Some scripts, options 
and parameters might have changed since that moment and we do not assure 100% resemblance with latest Github release.
```

Here, we show different steps for the retrieval of the data, analysis and visualization.

### Retrieve data

The data retrieval remains the same in the previous version. Data would be downloaded from NCBI. Necessary data is 
basically translated cds protein for each strain, genomic fasta files and genomic gene feature files (gff). We would also
need genbank files (gbk) for different analysis.

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
We will retrieve all information for each strain but we basically rely on the header information supplied by translated_cds.faa file

```
>lcl|FN554766.1_prot_CBG32835.1_1 [gene=thrA] [locus_tag=EC042_0001] [db_xref=GOA:D3H385,InterPro:IPR001048,InterPro:IPR001341,InterPro:IPR001342,InterPro:IPR002912,InterPro:IPR005106,InterPro:IPR011147,InterPro:IPR016040,InterPro:IPR018042,InterPro:IPR019811,InterPro:IPR027795,UniProtKB/TrEMBL:D3H385] [protein=bifunctional aspartokinase I/homoserine dehydrogenase I] [protein_id=CBG32835.1] [location=336..2798] [gbkey=CDS]
MRVLKFGGTSVANAERFLRVADILESNARQGQVATVLSAPAKITNHLVAMIEKTISGQDALPNISDAERIFAELLTGLAA
....

```
ATTENTION: DO NOT use file *protein.faa
   
Translated CDS contain directly translated coding sequence regions and sometimes proteins that are identically the same are collapsed into 1 entry into database, so a duplicated gene that has two different positions in the genome, two translated cds would only have one protein.
```

```

NOTE: It is possible to generate the same results from a gff file, but it is not implemented neither intended right now.
Please contact us for further explanation or clarification or to show interest as we might by thinking of implementing this feature. 

#### Example
See additional details and an example of this data retrieval process in the previous [Ecoli analysis example](https://github.com/molevol-ub/BacterialDuplicates/blob/master/Ecoli/Ecoli_genome.md#example)

## Gene duplication among strains analysis

As for the previous analysis in E.coli, we analyzed for each of the samples we retrieved from NCBI the duplicated 
genes within each genome.

During the development of this version, we include several new features and control for putative confounding results. E.g.

- We modified the parsing blast result parameters and set, not only lenght and similarity, 
but also e-values (<1e-05)and bit-score (>50) parameters as default cutoffs.

- We take into account that duplicates might be between any pair of available sequences. 

	* Chrm -> Chrm; 
	* Chrm -> Plasmid; 
	* Plasmid -> Plasmid; 
	* Contig1 -> Contig23214;
	* ...

- We also annotated accordingly pseudogenes to make some statistics later on.

- We additionally cleaned fasta headers that might producing errors due to character encoding.

- We implemented a virulence and resistance search for all proteins and also only among duplicated proteins.

### Explanation of the workflow

Most of these new implementations and corrections were done in script `generate_results_duplicates.pl` and in `parse_BLAST.pl`.
We also set the script `duplicate_search_bacteria.pl` to accomodate and generate all results resulting it as the
main pipeline of the project, generating up to 15 different steps

- 1: Make folder
- 2: Clean sequences
- 3: Makeblastdb
- 4: BLAST proteins
- 5: Parse BLAST
- 6: Generate duplicate results
- 7: Get fasta seqs results
- 8: BLAST duplicated proteins vs. CARD
- 9: Parse BLAST duplicated CARD results
- 10: BLAST ALL proteins vs. CARD
- 11: Parse BLAST ALL CARD results
- 12: BLAST proteins vs. VFDB
- 13: Parse BLAST VFDB results
- 14: BLAST ALL proteins vs. VFDB
- 15: Parse BLAST ALL VFDB results

Basically, for each strain a folder is generated containing all results. Sequences in **translated_cds.faa** file are indexed using
blast and search discarding autohits. A similarity and length cutoff of 85% is employed and duplicated genes are reported and grouped
accordingly. Then, duplicated proteins identified are retrieved from the total set. Now, only for the duplicated proteins identified and 
later for all the proteins, virulence and resistance databases are searched. The similarity cutoff is automatically relaxed to 35% following
homology search strategies ([PMID: 23749753](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3820096/)).

### Example of usage

```
perl scripts/perl/duplicate_search_bacteria.pl -fasta proteins.fasta -name example 
		-BLAST_path /path/to/BLAST/bin 
		-CARD_db /path/to/CARD_databases/blast_id_name
		-vfdb_db /path/to/VFDB_databases/blast_id_name
		[-CPU nCPU -sim 85 -len 85]

```

### Virulence and Resistance Databases retrieval

We downloaded the Virulence Factor Database ([VFDB](http://www.mgc.ac.cn/cgi-bin/VFs/search.cgi)) 
and the Comprehensive Antibiotic Resistance Database ([CARD](https://card.mcmaster.ca/)) databases. 

See additional details
on version in file `db_info.csv` available in the Supplementary information [section](#supplementary-information)


```
perl /users/jfsanchez/jfsh_scripts/rename_FASTA_seqs.pl protein_fasta_protein_homolog_model.fasta protein_fasta_protein_homolog_model_renamed.fasta card REPLACE > protein_fasta_protein_homolog_model_renamed.conversion.txt 
makeblastdb -in protein_fasta_protein_homolog_model_renamed.fasta -dbtype prot -input_type fasta -out card_db
```

```
perl /users/jfsanchez/jfsh_scripts/rename_FASTA_seqs.pl VFDB_setB_pro.fas VFDB_setB_pro.renamed.fas VFDB REPLACE > VFDB_setB_pro.renamed.conversion.txt 
makeblastdb -in VFDB_setB_pro.renamed.fas -dbtype prot -input_type fasta -out VFDB_db
```


```
for i in `dir ../1.data/NCBI_refseq_data/`; 
do 
	file=`readlink -f ../1.data/NCBI_refseq_data/$i/*translated_cds.faa`; 
	echo "Processing file: " $file; echo "###############"; 
	echo "qsub -cwd -V -q h0107.q -pe ompi8 2 -b y perl BacterialDuplicates/scripts/duplicate_search_bacteria.pl \
		-fasta $file -script_path BacterialDuplicates/scripts/parse_BLAST.pl -name $i -BLAST_path /soft/ncbi-blast-2.4.0/bin  \
		-CPU 2"; echo "#################"; 
	echo "";  
done > commands_sent.txt
```



## Supplementary information

Strains analyzed in the corresponding analysis and paper are deposited in files `Gram_positive/data/*samples.csv`.
See additional details [here](https://github.com/molevol-ub/BacterialDuplicates/blob/master/Gram_positive/data/)

## Citation

Gene Duplications in the Genomes of Staphylococci and Enterococci. Sanchez-Herrero JF., Bernabeu M., Prieto A., Hüttener M. and Juárez A. **Front. Mol. Biosci.** *2020 7:160*. https://doi.org/10.3389/fmolb.2020.00160
