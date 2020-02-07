#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
##############################################################################
##	Jose Fco. Sanchez Herrero, 06/02/2020 jfsanchezherrero@ub.edu			##
##############################################################################

my $duplicates_file = $ARGV[0];
my $interpro_bin = $ARGV[1];
my $name = $ARGV[2];

if (!@ARGV) {print "Usage:\nperl $0 duplicate_relations translated_cds output\n";exit();}

my $intepro_command = $interpro_bin "-i ".$duplicates_file." -appl TIGRFAM,SFLD,SUPERFAMILY,Gene3D,Hamap,Coils,ProSiteProfiles,SMART,PRINTS,ProSitePatterns,Pfam,ProDom,MobiDBLite,PIRSF -b $name -pa -dp -goterms";
print "System call: $intepro_command\n";
system($intepro_command);




__END__
06/02/2020 19:02:25:252 Welcome to InterProScan-5.31-70.0
usage: java -XX:+UseParallelGC -XX:ParallelGCThreads=2 -XX:+AggressiveOpts
            -XX:+UseFastAccessorMethods -Xms128M -Xmx2048M -jar
            interproscan-5.jar


Please give us your feedback by sending an email to

interhelp@ebi.ac.uk

 -appl,--applications <ANALYSES>            Optional, comma separated list
                                            of analyses.  If this option
                                            is not set, ALL analyses will
                                            be run.
 -b,--output-file-base <OUTPUT-FILE-BASE>   Optional, base output filename
                                            (relative or absolute path).
                                            Note that this option, the
                                            --output-dir (-d) option and
                                            the --outfile (-o) option are
                                            mutually exclusive.  The
                                            appropriate file extension for
                                            the output format(s) will be
                                            appended automatically. By
                                            default the input file
                                            path/name will be used.
 -cpu,--cpu <CPU>                           Optional, number of cores for
                                            inteproscan.
 -d,--output-dir <OUTPUT-DIR>               Optional, output directory.
                                            Note that this option, the
                                            --outfile (-o) option and the
                                            --output-file-base (-b) option
                                            are mutually exclusive. The
                                            output filename(s) are the
                                            same as the input filename,
                                            with the appropriate file
                                            extension(s) for the output
                                            format(s) appended
                                            automatically .
 -dp,--disable-precalc                      Optional.  Disables use of the
                                            precalculated match lookup
                                            service.  All match
                                            calculations will be run
                                            locally.
 -dra,--disable-residue-annot               Optional, excludes sites from
                                            the XML, JSON output
 -f,--formats <OUTPUT-FORMATS>              Optional, case-insensitive,
                                            comma separated list of output
                                            formats. Supported formats are
                                            TSV, XML, JSON, GFF3, HTML and
                                            SVG. Default for protein
                                            sequences are TSV, XML and
                                            GFF3, or for nucleotide
                                            sequences GFF3 and XML.
 -goterms,--goterms                         Optional, switch on lookup of
                                            corresponding Gene Ontology
                                            annotation (IMPLIES -iprlookup
                                            option)
 -help,--help                               Optional, display help
                                            information
 -i,--input <INPUT-FILE-PATH>               Optional, path to fasta file
                                            that should be loaded on
                                            Master startup. Alternatively,
                                            in CONVERT mode, the
                                            InterProScan 5 XML file to
                                            convert.
 -iprlookup,--iprlookup                     Also include lookup of
                                            corresponding InterPro
                                            annotation in the TSV and GFF3
                                            output formats.
 -ms,--minsize <MINIMUM-SIZE>               Optional, minimum nucleotide
                                            size of ORF to report. Will
                                            only be considered if n is
                                            specified as a sequence type.
                                            Please be aware of the fact
                                            that if you specify a too
                                            short value it might be that
                                            the analysis takes a very long
                                            time!
 -o,--outfile <EXPLICIT_OUTPUT_FILENAME>    Optional explicit output file
                                            name (relative or absolute
                                            path).  Note that this option,
                                            the --output-dir (-d) option
                                            and the --output-file-base
                                            (-b) option are mutually
                                            exclusive. If this option is
                                            given, you MUST specify a
                                            single output format using the
                                            -f option.  The output file
                                            name will not be modified.
                                            Note that specifying an output
                                            file name using this option
                                            OVERWRITES ANY EXISTING FILE.
 -pa,--pathways                             Optional, switch on lookup of
                                            corresponding Pathway
                                            annotation (IMPLIES -iprlookup
                                            option)
 -t,--seqtype <SEQUENCE-TYPE>               Optional, the type of the
                                            input sequences (dna/rna (n)
                                            or protein (p)).  The default
                                            sequence type is protein.
 -T,--tempdir <TEMP-DIR>                    Optional, specify temporary
                                            file directory (relative or
                                            absolute path). The default
                                            location is temp/.
 -version,--version                         Optional, display version
                                            number
 -vtsv,--output-tsv-version                 Optional, includes a TSV
                                            version file along with any
                                            TSV output (when TSV output
                                            requested)
Copyright © EMBL European Bioinformatics Institute, Hinxton, Cambridge,
UK. (http://www.ebi.ac.uk) The InterProScan software itself is provided
under the Apache License, Version 2.0
(http://www.apache.org/licenses/LICENSE-2.0.html). Third party components
(e.g. member database binaries and models) are subject to separate
licensing - please see the individual member database websites for
details.

Available analyses:
                      TIGRFAM (15.0) : TIGRFAMs are protein families based on Hidden Markov Models or HMMs
                         SFLD (4) : SFLDs are protein families based on Hidden Markov Models or HMMs
                  SUPERFAMILY (1.75) : SUPERFAMILY is a database of structural and functional annotation for all proteins and genomes.
                       Gene3D (4.2.0) : Structural assignment for whole genes and genomes using the CATH domain structure database
                        Hamap (2018_03) : High-quality Automated and Manual Annotation of Microbial Proteomes
                        Coils (2.2.1) : Prediction of Coiled Coil Regions in Proteins
              ProSiteProfiles (2018_02) : PROSITE consists of documentation entries describing protein domains, families and functional sites as well as associated patterns and profiles to identify them
                        SMART (7.1) : SMART allows the identification and analysis of domain architectures based on Hidden Markov Models or HMMs
                          CDD (3.16) : Prediction of CDD domains in Proteins
                       PRINTS (42.0) : A fingerprint is a group of conserved motifs used to characterise a protein family
              ProSitePatterns (2018_02) : PROSITE consists of documentation entries describing protein domains, families and functional sites as well as associated patterns and profiles to identify them
                         Pfam (31.0) : A large collection of protein families, each represented by multiple sequence alignments and hidden Markov models (HMMs)
                       ProDom (2006.1) : ProDom is a comprehensive set of protein domain families automatically generated from the UniProt Knowledge Database.
                   MobiDBLite (2.0) : Prediction of disordered domains Regions in Proteins
                        PIRSF (3.02) : The PIRSF concept is being used as a guiding principle to provide comprehensive and non-overlapping clustering of UniProtKB sequences into a hierarchical order to reflect their evolutionary relationships.

Deactivated analyses:
        SignalP_GRAM_NEGATIVE (4.1) : Analysis SignalP_GRAM_NEGATIVE is deactivated, because the resources expected at the following paths do not exist: bin/signalp/4.1/signalp
                  SignalP_EUK (4.1) : Analysis SignalP_EUK is deactivated, because the resources expected at the following paths do not exist: bin/signalp/4.1/signalp
                        TMHMM (2.0c) : Analysis TMHMM is deactivated, because the resources expected at the following paths do not exist: bin/tmhmm/2.0c/decodeanhmm, data/tmhmm/2.0c/TMHMM2.0c.model
        SignalP_GRAM_POSITIVE (4.1) : Analysis SignalP_GRAM_POSITIVE is deactivated, because the resources expected at the following paths do not exist: bin/signalp/4.1/signalp
                      PANTHER (12.0) : Analysis Panther is deactivated, because the resources expected at the following paths do not exist: data/panther/12.0/panther.hmm, data/panther/12.0/names.tab
                      Phobius (1.01) : Analysis Phobius is deactivated, because the resources expected at the following paths do not exist: bin/phobius/1.01/phobius.pl

