## load modules
library(BioCircos)
library(colorspace)

## get sequence lengths
###awk '{if ($3=="region"){print $1","$5}}' ../../../1.data/Enterococcus/GCA_000007785.1_V583/GCF_000007785.1_ASM778v1_genomic.gff

########################
## Functions
########################
check_dot_refName <- function(name) {
  ifelse(grepl('.', name),return(unlist(strsplit(name,"\\."))[1]),return(name))
}
##
check_seqs <- function(list_items) {
  l <- NULL
  for (i in list_items) l<-c(l, check_dot_refName(i))
  return(unlist(l))
}
##
to_list <- function(df) {
  l <- list()
  for (i in rownames(df)) l[[ check_dot_refName(i) ]] <- df[i, ]    
  return(l)
}
##
duplicate_group_tracks <- function(dataF) {
  ##########################
  ## for each duplicated group
  ##########################
  
  ## sort by chr & start position first
  dataF <- dataF[order(dataF$chr, dataF$start ),]
  
  ## sequence location
  seqs_vector <- check_seqs(as.character(dataF$chr))
  links_chromosomes_1 = head(seqs_vector, -1)
  links_chromosomes_2 = seqs_vector[-1]
  
  ## start position
  starts_vector <- dataF$start
  links_pos_1 = head(starts_vector, -1)
  links_pos_2 = starts_vector[-1]
  
  links_labels = paste0("dup", dataF$group[1])

  ##########################
  ## create a track with the links
  ##########################
  DupTrack = BioCircosLinkTrack(trackname = paste0('myDupTrack_',links_labels), 
                                
                                opacities=0.3,
                                
                                ## link position 1
                                gene1Chromosomes = links_chromosomes_1, 
                                gene1Starts = links_pos_1, 
                                gene1Ends = links_pos_1 + 1, #we set the end position of each link: start+1
                                
                                ## link position 2
                                gene2Chromosomes = links_chromosomes_2, 
                                gene2Starts = links_pos_2, 
                                gene2Ends = links_pos_2 + 1, #we set the end position of each link: start+1
                                
                                ## aesthetics
                                color = '#eeb242', #dataF$color[1], ## add the color from palette generated
                                displayAxis = T,
                                
                                ## Labels: add name but do not show in track, only when mouse hovers
                                displayLabel = F, 
                                labels = links_labels, 
                                
                                ## geneXNames
                                #gene1Names = c("test"), 
                                #gene2Names = c("test2"),
                                
                                maxRadius = 0.855, width = 0.5)

  ##
  return(DupTrack)
}
##
parse_data <- function(bed_info) {
  ## add header
  colnames(bed_info) <- c("group", "chr", "start", "end", "gene_name", "description", "strand")
  
  ## create more informative:
  bed_info$Gene <- 'Gene: '
  cols <- c("Gene", "gene_name")
  bed_info$new_gene_name <- apply(bed_info[, cols ], 1, paste, collapse = "" )
  
  cols <- c("new_gene_name", "description")
  bed_info$new_gene_name <- apply(bed_info[, cols ], 1, paste, collapse = ": " )
  
  cols <- c("new_gene_name", "strand")
  bed_info$new_gene_name <- apply(bed_info[, cols ], 1, paste, collapse = "; Strand: (" )
  
  cols <- c("new_gene_name", "group")
  bed_info$new_gene_name <- apply(bed_info[, cols ], 1, paste, collapse = "); Duplicate Group: " )
  
  ## create color palette
  rbPal <- colorRampPalette(rainbow_hcl(length(unique(bed_info$group))))
  bed_info$color <- 'test'
  
  for (i in 1:length(unique(bed_info$group)))
    bed_info[bed_info$group==i,]$color <- rbPal(length(bed_info$group))[i]
  
  return(bed_info)
}
## 
add_phages <- function(phages_phispy) {
  ##########################
  #### add PhiSpy regions
  ##########################
  head(phages_phispy)
  phages_phispy$Phage <- 'Phage Region: '
  cols_phispy <- c("Phage", "X")
  phages_phispy$new_col <- apply( phages_phispy[ , cols_phispy] , 1 , paste , collapse = "" )
  
  # details
  arcs_chromosomes_phages = check_seqs(as.character(phages_phispy$Contig))  # Chromosomes on which the arcs should be displayed
  arcs_begin_phages = phages_phispy$Start
  arcs_end_phages = phages_phispy$End
  
  # create tracklist
  phage_tracklist = BioCircosArcTrack('myPhageTrack', 
                                      labels = as.character(phages_phispy$new_col),
                                      arcs_chromosomes_phages, arcs_begin_phages, arcs_end_phages,
                                      minRadius = 0.955, maxRadius = 0.99, colors = "blue")
  
  background_tracklist = BioCircosBackgroundTrack("myBackgroundTrack_phispy", 
                                                  minRadius = 0.95, maxRadius = 0.995,
                                                  borderColors = "black", borderSize = 0.3, 
                                                  fillColors = "#FF000000")  
  tracks = background_tracklist + phage_tracklist
  return(tracks)
}
##
add_dups <- function(bed_info) {
  ##########################
  ## add duplicate information
  ##########################
  ## connection between duplicates is done using BioCircosLinkTrack
  ## and duplicate genes are show using genes_tracklist
  
  ##########################
  ## add as many tracks as groups
  ##########################
  ## split by group
  bed_info_split <- split(bed_info, bed_info$group)
  
  tracks_duplicates <- NULL
  for (i in bed_info_split)
    tracks_duplicates = tracks_duplicates + duplicate_group_tracks(i)
  
  return(tracks_duplicates)
}
##
add_genes_strand <- function(bed_info) {
  ##########################
  ## add genes per strand
  ##########################
  bed_info_strand <- split(bed_info, bed_info$strand)
  
  arcs_chromosomes_genes_strand_pos = check_seqs(as.character(bed_info_strand$`+`$chr))  # Chromosomes on which the arcs should be displayed
  arcs_begin_genes_strand_pos = bed_info_strand$`+`$start
  arcs_end_genes_strand_pos = bed_info_strand$`+`$end
  
  ## Positive strand genes
  gene_pos_tracklist = BioCircosArcTrack('myGenePosStrandTrack', 
                                         labels = as.character(bed_info_strand$`+`$new_gene_name),
                                         arcs_chromosomes_genes_strand_pos, 
                                         arcs_begin_genes_strand_pos, 
                                         arcs_end_genes_strand_pos,
                                         minRadius = 0.94, maxRadius = 0.905, colors = "#ee42d4")
  pos_background_tracklist = BioCircosBackgroundTrack("myBackgroundTrack_posGenes", 
                                                                         minRadius = 0.90, maxRadius = 0.945,
                                                                         borderColors = "black", borderSize = 0.3, 
                                                                         fillColors = "#FF000000")  
  
  arcs_chromosomes_genes_strand_neg = check_seqs(as.character(bed_info_strand$`-`$chr))  # Chromosomes on which the arcs should be displayed
  arcs_begin_genes_strand_neg = bed_info_strand$`-`$start
  arcs_end_genes_strand_neg = bed_info_strand$`-`$end
  
  ## Negative strand genes
  gene_neg_tracklist = BioCircosArcTrack('myGeneNegStrandTrack', 
                                         labels=as.character(bed_info_strand$`-`$new_gene_name),
                                         arcs_chromosomes_genes_strand_neg, 
                                         arcs_begin_genes_strand_neg, 
                                         arcs_end_genes_strand_neg,
                                         minRadius = 0.86, maxRadius = 0.89, colors = "#42eeb2")
  
  neg_background_tracklist = BioCircosBackgroundTrack("myBackgroundTrack_islandphat", 
                                                                         minRadius = 0.855, maxRadius = 0.895,
                                                                         borderColors = "black", borderSize = 0.3, 
                                                                         fillColors = "#FF000000")  
  
  tracks = gene_pos_tracklist + pos_background_tracklist + gene_neg_tracklist + neg_background_tracklist
  return(tracks)
}
##
create_BioCircos <- function(seq_lengths, bed_info_file, phages_phispy) {
  
  ## parse information
  myGenome= to_list(seq_lengths)
  myGenome
  bed_info <- parse_data(bed_info_file)
  
  ###############################
  ## add tracks
  ###############################
  
  ## phages
  tracks <- add_phages(phages_phispy)
  
  ## duplicated gene connections
  tracks <- tracks + add_dups(bed_info)
  
  ## duplicated genes bands
  tracks <- tracks + add_genes_strand(bed_info)
  ##########################
  
  ##########################
  ## plot it
  ##########################
  BioCircos(genome = myGenome, 
            tracklist = tracks, 
            
            ## Aesthetics
            genomeFillColor = c('#ee5c42', '#42d4ee', '#42d4ee', '#42d4ee', '#42d4ee'), #"Spectral", 
            chrPad = 0.1, ## distance between seqs 
            displayGenomeBorder = T, 
            
            ## Ticks
            genomeTicksScale = 1e+6, #genomeTicksLen = 10, #genomeTicksTextSize = 12, 
            genomeTicksDisplay = TRUE,
            genomeTicksLen = 5, 
            genomeTicksColor = "#000",
            genomeTicksTextSize = "0.6em", 
            genomeTicksTextColor = "#000",
            
            ## Labels
            genomeLabelTextSize = 12, 
            genomeLabelOrientation = 60, 
            genomeLabelDy = 45, 
            genomeLabelDx = 0.1, 
            
            zoom = TRUE, #TEXTModuleDragEvent = TRUE,
  )
}
########################

###############################
## retrieve info & plot
###############################
supData_folder <- "/home/jfsanchez/DATA/BacterialDuplicates/SupData"
## V583
seq_lengths_V583 <- read.csv(paste0(supData_folder, "/SupData/Efaecalis/ref_strains/V583-seq_lengths.csv"), header=FALSE, row.names = 1)
phages_phispy_V583 <- read.csv(paste0(supData_folder, "/SupData/Efaecalis/ref_strains/GCF_000007785.1_Efaecalis_V583-prophage-coordinates.tsv"), sep="\t", header=T)
bed_info_file_V583 <- read.table(paste0(supData_folder, "/SupData/Efaecalis/ref_strains/GCA_000007785.1_V583_BLAST.out.coordinates.csv"),sep=",", header=FALSE)
##
create_BioCircos(seq_lengths = seq_lengths_V583, 
                 phages_phispy = phages_phispy_V583, 
                 bed_info_file = bed_info_file_V583)

## 6E6
seq_lengths_6E6 <- read.csv(paste0(supData_folder, "/SupData/Efaecium/ref_strains/Efaecium_6E6.csv"), header=FALSE, row.names = 1)
phages_phispy_6E6 <- read.csv(paste0(supData_folder, "/SupData/Efaecium/ref_strains/GCF_001518735_Efaecium_6E6-prophage-coordinates.tsv"), sep="\t", header=T)
bed_info_file_6E6 <- read.table(paste0(supData_folder, "/SupData/Efaecium/ref_strains/GCA_001518735.1_Strain_6E6_BLAST.out.coordinates.csv"),sep=",", header=FALSE)
##
create_BioCircos(seq_lengths =  seq_lengths_6E6, 
                 phages_phispy = phages_phispy_6E6, 
                 bed_info_file = bed_info_file_6E6)

# Newman
seq_lengths_Newman <- read.csv(paste0(supData_folder, "/SupData/Saureus/ref_strains/Newman_lengths.csv"), header=FALSE, row.names = 1)
phages_phispy_Newman <- read.csv(paste0(supData_folder, "/SupData/Saureus/ref_strains/GCF_000010465.1_Saureus_Newman-prophage-coordinates.tsv"), sep="\t", header=T)
bed_info_file_Newman <- read.table(paste0(supData_folder, "/SupData/Saureus/ref_strains/GCA_000010465.1_Saureus_Newman_BLAST.out.coordinates.csv"),sep=",", header=FALSE)
##
create_BioCircos(seq_lengths = seq_lengths_Newman, 
                 phages_phispy = phages_phispy_Newman, 
                 bed_info_file = bed_info_file_Newman)

# ST398
seq_lengths_ST398 <- read.csv(paste0(supData_folder, "/SupData/Saureus/ref_strains2/ST398.csv"), header=FALSE, row.names = 1)
phages_phispy_ST398 <- read.csv(paste0(supData_folder, "/SupData/Saureus/ref_strains2/GCF_000009585.1_Saureus_ST398-prophage-coordinates.tsv"), sep="\t", header=T)
bed_info_file_ST398 <- read.table(paste0(supData_folder, "/SupData/Saureus/ref_strains2/GCA_000009585.1_Saureus_typestrain_ST398_BLAST.out.coordinates.csv"),sep=",", header=FALSE)
##
create_BioCircos(seq_lengths = seq_lengths_ST398, 
                 phages_phispy = phages_phispy_ST398, 
                 bed_info_file = bed_info_file_ST398)
