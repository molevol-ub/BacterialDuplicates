## load modules
library(colorspace)
library(Sushi)

## https://www.bioconductor.org/packages/release/bioc/vignettes/Sushi/inst/doc/Sushi.pdf

# the coordinates file desired
coordinates_duplicates_file <-"C:/Users/josed/Desktop/Developer/BacterialDuplicates/search_Resutls/Newman/GCA_000010465.1_Saureus_Newman_BLAST.out.coordinates.csv"
phages_phispy <- read.csv("/Users/josed/Desktop/Developer/BacterialDuplicates/dataFromHercules/PhiSpy/GCF_000010465.1_Saureus_Newman-prophage-coordinates.tsv", sep="\t", header=T)

## split data by sequence
bed_info <- read.table(coordinates_duplicates_file,sep=",",header=FALSE)
bed_info_split <- split(bed_info,bed_info$V2)

# create groups 
list_Ids <- unique(as.numeric(bed_info_split$NC_009641.1$V1))
list_inte <- length(list_Ids)
group <- bed_info_split$NC_009641.1$V1
position_chr <- bed_info_split$NC_009641.1$V3

## create colors for duplicates
rbPal <- colorRampPalette(rainbow_hcl(list_inte))
list <- rbPal(list_inte)
for (ids in list) {
  bed_info_split$NC_009641.1$col <- list[bed_info_split$NC_009641.1$V1]
}

## create data for phages track
chrom = "NC_009641.1"
chromstart = 0
chromend = max(bed_info_split$NC_009641.1$V3)

phages_bed <- NULL
phages_bed <- phages_phispy[c("Contig")]
phages_bed['chrom'] <- phages_bed['Contig']
phages_bed['start'] <- as.character(phages_phispy$Start)
phages_bed['end'] <- as.character(phages_phispy$End)
phages_bed['strand'] <- c(1,-1,1,-1,1)
phages_bed$Contig <- NULL

## create layout
layout(matrix(c(1,2,2,2), nrow = 4, ncol = 1, byrow = TRUE))

# Add boxplots to a scatterplot
par(fig=c(0,1,0,0.9), new=TRUE)
## plot duplicates
plot(
  position_chr, 
  group,
  pch=ifelse(bed_info_split$NC_009641.1$V7=="-",1,16),
  col=ifelse(grepl("transpos",bed_info_split$NC_009641.1$V6),"black",bed_info_split$NC_009641.1$col)
)
for (ids in list_Ids) {
  lines(position_chr[group==ids], group[group==ids], type="l",lty=4,lwd=0.5)
}
text(position_chr, group, labels=group, pos=2, cex=0.7)

par(fig=c(0,1,0.68,1), new=TRUE)
plotBed(beddata=phages_bed,
        chrom = chrom, 
        chromstart=chromstart, 
        chromend=chromend, 
        type = 'region'
)
#labelgenome('',chromstart,chromend,n=3, side = 3, scale = 'Mb')
mtext("Test duplicated plot", side=3, outer=TRUE, line=-3)


phages_bed

#par(mfrow=c(2,1))
## plot track
dev.off()
plotBed(beddata=phages_bed,
        chrom = chrom, 
        chromstart=chromstart, 
        chromend=chromend, 
        type = 'region'
        )

labelgenome('',chromstart,chromend,n=3, side = 3, scale = 'Mb')

## plot duplicates
plot(
  position_chr, 
  group,
  pch=ifelse(bed_info_split$NC_009641.1$V7=="-",1,16),
  col=ifelse(grepl("transpos",bed_info_split$NC_009641.1$V6),"black",bed_info_split$NC_009641.1$col)
)

for (ids in list_Ids) {
  lines(position_chr[group==ids], group[group==ids], type="l",lty=4,lwd=0.5)
}
text(position_chr, group, labels=group, pos=2, cex=0.7)




