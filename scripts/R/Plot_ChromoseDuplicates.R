## load modules
library(colorspace)

data_split <- function(info_df) {
  
  ## split data by duplicate group id
  bed_info_split <- split(info_df, info_df$V1)
  info_df$V8 <- 'test'
  info_df$V9 <- 'test'
  
  ## check if a duplicate is distributed in >1 sequences (== plasmids)
  for (i in bed_info_split)
    if (length(unique(i$V2)) > 1) {
      info_df[info_df$V1==as.factor(i$V1[1]),]$V8 <- 'yes'
      info_df[info_df$V1==as.factor(i$V1[1]),]$V9 <- 17
    } else {
      info_df[info_df$V1==as.factor(i$V1[1]),]$V8 <- 'no'
      info_df[info_df$V1==as.factor(i$V1[1]),]$V9 <- 16
    }

  ## assign pch symbol to each gene coordinate according to
  ## strand and if in plasmid or not.
  # chromosome strand +:  16 (filled circle)
  # chromosome strand -: 1 (empty circle)
  # plasmid strand +: 17 (filled triangle)
  # plasmid strand -: 2 (empty triangle)
  fun1 <- function(x, y) if (x == 17) { if (y=='-') { 2 } else { 17} } else { if (y=='-') { 1 } else { 16} }
  info_df$V10 <- 'test'
  info_df$V10 <- mapply(fun1, info_df$V9, info_df$V7)
  
  ## split data by sequence id
  bed_info_split2 <- split(info_df, info_df$V2)
  
  return(bed_info_split2)
}

plot_duplicates <- function(df) {
  # create variables
  list_Ids <- unique(as.numeric(df$V1))
  list_inte <- max(list_Ids) ##length(list_Ids)
  group <- df$V1
  position_chr <- df$V3
  strand <- as.character(df$V7)
  form <- as.numeric(df$V10)

  ## create color palette for duplicates groups
  rbPal <- colorRampPalette(rainbow_hcl(list_inte))
  list <- rbPal(list_inte)
  for (ids in list) {
    df$col <- list[df$V1]
  }
  
  ## plot duplicates
  plot(
    position_chr, 
    group,
    pch=form,
    col=ifelse(grepl("transpos",df$V6),"black",df$col)
  )
  for (ids in list_Ids) {
    lines(position_chr[group==ids], group[group==ids], type="l", lty=4, lwd=0.5)
  }
  text(position_chr, group, labels=group, pos=2, cex=0.7)
}

## supplementary data folder path (modify it as necessary)
supData_folder <- "/home/jfsanchez/DATA/BacterialDuplicates/SupData"

#####################################################
## Saureus Newman
#####################################################
# the coordinates file desired
coordinates_duplicates_file_Saures <- paste0(supData_folder, "/SupData/Saureus/ref_strains/GCA_000010465.1_Saureus_Newman_BLAST.out.coordinates.csv")

## split data by sequence
bed_info_Saures <- read.table(coordinates_duplicates_file_Saures, sep=",",header=FALSE)
bed_info_split_Saures <- data_split(bed_info_Saures)

## select the sequence to plot: main chromosome only
plot_duplicates(bed_info_split_Saures$NC_009641.1)
#####################################################

#####################################################
## Saureus ST398
#####################################################
# the coordinates file desired
coordinates_duplicates_file_Saures_ST398 <- paste0(supData_folder, "/SupData/Saureus/ref_strains2/GCA_000009585.1_Saureus_typestrain_ST398_BLAST.out.coordinates.csv")

## split data by sequence
bed_info_Saures_ST398 <- read.table(coordinates_duplicates_file_Saures_ST398, sep=",",header=FALSE)
bed_info_split_Saures_ST398 <- data_split(bed_info_Saures_ST398)

## select the sequence to plot: main chromosome only
plot_duplicates(bed_info_split_Saures_ST398$NC_017333.1)
#####################################################



#####################################################
## Efaecalis V583
#####################################################
# the coordinates file desired
coordinates_duplicates_file_Efaecalis <-  paste0(supData_folder, "/SupData/Efaecalis/ref_strains/GCA_000007785.1_V583_BLAST.out.coordinates.csv")

## split data by sequence
bed_info_Efaecalis <- read.table(coordinates_duplicates_file_Efaecalis, sep=",", header=FALSE)
bed_info_split_Efaecalis <- data_split(bed_info_Efaecalis)

## select the sequence to plot: main chromosome only
plot_duplicates(bed_info_split_Efaecalis$NC_004668.1)
#####################################################

#####################################################
## Efaecium 6E6
#####################################################
# the coordinates file desired
  coordinates_duplicates_file_Efaecium <-  paste0(supData_folder, "/SupData/Efaecium/ref_strains/GCA_001518735.1_Strain_6E6_BLAST.out.coordinates.csv")

## split data by sequence
bed_info_Efaecium <- read.table(coordinates_duplicates_file_Efaecium, sep=",",header=FALSE)
bed_info_split_Efaecium <- data_split(bed_info_Efaecium)

## select the sequence to plot: main chromosome only
plot_duplicates(bed_info_split_Efaecium$NZ_CP013994.1)
#####################################################
