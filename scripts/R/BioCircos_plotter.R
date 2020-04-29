## http://www.nanocell.cl/making-a-simply-circos-representation-of-a-bacterial-genome/
## http://circos.ca/tutorials/course/tunis/handouts/session-3.pdf
## https://dbsloan.github.io/TS2019/exercises/circos.html
## https://cran.rstudio.com/web/packages/BioCircos/vignettes/BioCircos.html

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

check_seqs <- function(list_items) {
  l <- NULL
  for (i in list_items) l<-c(l, check_dot_refName(i))
  return(unlist(l))
}

to_list <- function(df) {
  l <- list()
  for (i in rownames(df)) l[[ check_dot_refName(i) ]] <- df[i, ]    
  return(l)
}

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
  
  print (list_colors[dataF$group[1]])
  
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
                                color = list_colors[dataF$group[1]], ## add the color from palette generated
                                displayAxis = T,
                                
                                ## Labels: add name but do not show in track, only when mouse hovers
                                displayLabel = F, 
                                labels = links_labels, 
                                
                                ## geneXNames
                                #gene1Names = c("test"), 
                                #gene2Names = c("test2"),
                                
                                maxRadius = 0.83)
  ##########################
  
  ##########################
  ## add points for each gene
  ##########################
  ## For each point, add a point with same color as link
  ## circle: strand +
  ## rect: strand -
  ## triangle: 
  ## x: pseudo
  #points_values = 0:4
  #genes_tracklist = BioCircosSNPTrack(trackname = paste0('myGeneTrack_',links_labels), 
  #                                    chromosomes = seqs_vector, 
  #                                    positions = starts_vector, 
  #                                   
  #                                   values = points_values, 
  #                                   colors = list_colors[dataF$group[1]], ## add the color from palette generated
  #                                   shape = "circle", 
  #                                   minRadius = 0.87, maxRadius = 0.87, size = 4)
  ##################
  
  ## add tracks and return them
  #tracks = DupTrack + genes_tracklist
  return(DupTrack)
}
########################

###############################
## retrieve info
###############################

## V583
#seq_lengths <- read.csv("/Users/josed/Desktop/Developer/BacterialDuplicates/BioCircos/V583-seq_lengths.txt", header=FALSE, row.names = 1)
#phages_phispy <- read.csv("/Users/josed/Desktop/Developer/BacterialDuplicates/dataFromHercules/PhiSpy/GCF_000007785.1_Efaecalis_V583-prophage-coordinates.tsv", sep="\t", header=T)
#bed_info <- read.table("/Users/josed/Desktop/Developer/BacterialDuplicates/search_Resutls/V583/GCA_000007785.1_V583_BLAST.out.coordinates.csv",sep=",",header=FALSE)
#regions_islandpath <- read.csv("/Users/josed/Desktop/Developer/BacterialDuplicates/dataFromHercules/IslandPath/Efaecalis_V583-islandpath.txt", sep="\t", header=T)

## 6E6
seq_lengths <- read.csv("/Users/josed/Desktop/Developer/BacterialDuplicates/BioCircos/Efaecium_6E6.csv", header=FALSE, row.names = 1)
phages_phispy <- read.csv("/Users/josed/Desktop/Developer/BacterialDuplicates/dataFromHercules/PhiSpy/GCF_001518735_Efaecium_6E6-prophage-coordinates.tsv", sep="\t", header=T)
bed_info <- read.table("/Users/josed/Desktop/Developer/BacterialDuplicates/search_Resutls/6E6/GCA_001518735.1_Strain_6E6_BLAST.out.coordinates.csv",sep=",",header=FALSE)
regions_islandpath <- read.csv("/Users/josed/Desktop/Developer/BacterialDuplicates/dataFromHercules/IslandPath/Efaecium_6E6-islandpath.txt", sep="\t", header=T)

# Newman
#seq_lengths <- read.csv("/Users/josed/Desktop/Developer/BacterialDuplicates/BioCircos/Saures_Newman.csv.txt", header=FALSE, row.names = 1)
#phages_phispy <- read.csv("/Users/josed/Desktop/Developer/BacterialDuplicates/dataFromHercules/PhiSpy/GCF_000010465.1_Saureus_Newman-prophage-coordinates.tsv", sep="\t", header=T)
#bed_info <- read.table("/Users/josed/Desktop/Developer/BacterialDuplicates/search_Resutls/Newman/GCA_000010465.1_Saureus_Newman_BLAST.out.coordinates.csv",sep=",",header=FALSE)
#regions_islandpath <- read.csv("/Users/josed/Desktop/Developer/BacterialDuplicates/dataFromHercules/IslandPath/Saureus_Newman-islandpath.txt", sep="\t", header=T)


myGenome= to_list(seq_lengths)
myGenome


## add header
colnames(bed_info) <- c("group", "chr", "start", "end", "gene_name", "description", "strand")

## create more informative:
bed_info$Gene <- 'Gene: '
cols <- c("Gene", "gene_name")
bed_info$new_gene_name <- apply( bed_info[ , cols ] , 1 , paste , collapse = "" )

cols <- c("new_gene_name", "description")
bed_info$new_gene_name <- apply( bed_info[ , cols ] , 1 , paste , collapse = ": " )

cols <- c("new_gene_name", "strand")
bed_info$new_gene_name <- apply( bed_info[ , cols ] , 1 , paste , collapse = "; Strand: (" )

cols <- c("new_gene_name", "group")
bed_info$new_gene_name <- apply( bed_info[ , cols ] , 1 , paste , collapse = "); Duplicate Group: " )


## split by group
bed_info_split <- split(bed_info, bed_info$group)

## create color palette
rbPal <- colorRampPalette(rainbow_hcl(length(bed_info$group)))
list_colors <- as.vector(rbPal(length(bed_info$group)))

###############################
## add tracks
###############################

###############################
## title track
###############################
#title_name <- "Enterococcus faecalis V583"
#title_track = BioCircosTextTrack('Title', title_name, size = "2em", weight = "bold", opacity = 0.5)

###############################
## add Arc Tracks

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

##########################

##########################
## add IslandPath regions
##########################
##head(regions_islandpath)
##
##regions_islandpath$ID <- 'Pathogenic Island: '
##cols_Islandpath <- c("ID", "X..GI")
##regions_islandpath$new_col <- apply( regions_islandpath[ , cols_Islandpath] , 1 , paste , collapse = "" )
##
# details
##arcs_chromosomes_islandpath = check_seqs(as.character(regions_islandpath$seq)) # Chromosomes on which the arcs should be displayed
##arcs_begin_islandpath = regions_islandpath$start
##arcs_end_islandpath = regions_islandpath$end
##
##pathogenic_tracklist = BioCircosArcTrack('Pathogenic_Islands', 
##                                       labels=as.character(regions_islandpath$new_col),
##                                         arcs_chromosomes_islandpath, arcs_begin_islandpath, arcs_end_islandpath, 
##                                         minRadius = 0.94, maxRadius = 0.905, colors = "red")
##background_tracklist = background_tracklist + BioCircosBackgroundTrack("myBackgroundTrack_islandphat", 
##                                                                       minRadius = 0.90, maxRadius = 0.945,
##                                                                       borderColors = "black", borderSize = 0.3, 
##                                                                       fillColors = "#FF000000")  
#########################

##########################
## add duplicate information
##########################
## connection between duplicates is done using BioCircosLinkTrack
## and duplicate genes are show using genes_tracklist

##########################
## add as many tracks as groups
##########################
#bed_info_split

tracks_duplicates <- NULL
for (i in bed_info_split)
  tracks_duplicates = tracks_duplicates + duplicate_group_tracks(i)
##########################

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
                                       minRadius = 0.94, maxRadius = 0.905, colors = "orange")
background_tracklist = background_tracklist + BioCircosBackgroundTrack("myBackgroundTrack_posGenes", 
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
                                       minRadius = 0.865, maxRadius = 0.895, colors = "green")

background_tracklist = background_tracklist + BioCircosBackgroundTrack("myBackgroundTrack_islandphat", 
                                                                       minRadius = 0.86, maxRadius = 0.9,
                                                                       borderColors = "black", borderSize = 0.3, 
                                                                       fillColors = "#FF000000")  
##########################


##########################
## add tracks
##########################
##tracklist2display = title_track + phage_tracklist + pathogenic_tracklist + tracks_duplicates
#tracklist2display = phage_tracklist + pathogenic_tracklist + gene_neg_tracklist + gene_pos_tracklist 
tracklist2display = phage_tracklist + 
  gene_neg_tracklist + 
  gene_pos_tracklist + background_tracklist + 
  tracks_duplicates

##########################
## plot it
##########################
BioCircos(genome = myGenome, 
          tracklist = tracklist2display, 
          
          ## Aesthetics
          genomeFillColor = "Spectral", 
          chrPad = 0.1, ## distance between seqs 
          displayGenomeBorder = T, 
          
          ## Ticks
          genomeTicksScale = 1e+6, #genomeTicksLen = 10, #genomeTicksTextSize = 12, 
          genomeTicksDisplay = TRUE,
          genomeTicksLen = 5, genomeTicksColor = "#000",
          genomeTicksTextSize = "0.6em", genomeTicksTextColor = "#000",
          
          ## Labels
          genomeLabelTextSize = 14, 
          genomeLabelOrientation = 90, genomeLabelDy = 45, genomeLabelDx = 0, 
          
          zoom = TRUE, #TEXTModuleDragEvent = TRUE,
          
)
