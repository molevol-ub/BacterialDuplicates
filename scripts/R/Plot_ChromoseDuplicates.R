library(colorspace)
#

bed_info <- read.table("coordinates.csv",sep=",",header=FALSE)
bed_info_split <- split(bed_info,bed_info$V2)

bed_info

# 
list_Ids <- unique(as.numeric(bed_info_split$chrm$V1))
list_inte <- length(list_Ids)
group <- bed_info_split$chrm$V1
position_chr <- bed_info_split$chrm$V3

bed_info_split$chrm
rbPal <- colorRampPalette(rainbow_hcl(list_inte))
list <- rbPal(list_inte)

for (ids in list) {
  bed_info_split$chrm$col <- list[bed_info_split$chrm$V1]
}

## plot for chrm
plot(
  position_chr, 
  group,
  pch=ifelse(bed_info_split$chrm$V7=="-",1,16),
  col=ifelse(grepl("transpos",bed_info_split$chrm$V6),"black",bed_info_split$chrm$col)
)

for (ids in list_Ids) {
  lines(position_chr[group==ids], group[group==ids], type="l",lty=4,lwd=0.5)
}
text(position_chr, group, labels=group, pos=2, cex=0.7)

