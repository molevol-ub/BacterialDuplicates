library(ggplot2)
library(tidyr)
library(dplyr)
library(stringr)
library(tidyverse)
library(ggpubr)

######################
##### Functions ######
######################
## parse data
split_into_multiple <- function(column, pattern = ",", into_prefix){
  cols <- str_split_fixed(column, pattern, n = Inf)
  # Sub out the ""'s returned by filling the matrix to the right, with NAs which are useful
  cols[which(cols == "")] <- NA
  cols <- as.tibble(cols)
  # name the 'cols' tibble as 'into_prefix_1', 'into_prefix_2', ..., 'into_prefix_m' 
  # where m = # columns of 'cols'
  m <- dim(cols)[2]
  
  names(cols) <- paste(into_prefix, 1:m, sep = "_")
  return(cols)
}
order_data <- function(summary_results, order_nwk_file, PhiSpy_results, resistance_results) {
  
  splitted <- summary_results %>% bind_cols(split_into_multiple(.$strain, "_", "fields"))
  summary_results_new <- splitted %>% separate(fields_2, c("Refseq", "version"), "\\.") %>% mutate(new_ID=paste0(fields_1, "_", Refseq)) 
  summary_results2 <- summary_results_new[1:9]
  summary_results2["index"] <- summary_results_new$new_ID
  row.names(summary_results2) <- summary_results_new$new_ID
  summary_results2["strain"] <- NULL
  summary_results2["fields_1"] <- NULL
  #print(head(summary_results_new))
  #print(head(summary_results2))
  
  ## PhiSpy results
  if (PhiSpy_results == "na") {
    print()
  } else {
    PhiSpy_results_df <-read.csv(PhiSpy_results, header=FALSE, stringsAsFactors = FALSE)
    rownames(PhiSpy_results_df) <- PhiSpy_results_df$V1
    summary_results2["phages"] <- PhiSpy_results_df[rownames(summary_results2),]$V2
  }
  #print(head(summary_results2))
  
  ## resistance results
  if (resistance_results == "na") {
    print()
  } else {
    resistance_results_df <-read.csv(resistance_results, header=FALSE, stringsAsFactors = FALSE)
    rownames(resistance_results_df) <- resistance_results_df$V1
    summary_results2["resistance"] <- resistance_results_df[rownames(summary_results2),]$V2
  }
  
  row.names(summary_results2) <- summary_results2$Refseq
  #print(head(summary_results2))
  
  new_order <- read.csv(order_nwk_file, header=FALSE, stringsAsFactors = FALSE)
  splitted_order <- new_order %>% bind_cols(split_into_multiple(.$V1, "_", "fields"))
  new_order1 <- splitted_order %>% 
    mutate(Refseq=fields_2) %>% 
    mutate(new_ID=paste0(fields_1, "_", Refseq)) 
  
  row.names(new_order1) <- new_order1$Refseq
  
  #print (new_order1)
  summary_results3 <- summary_results2[match(rownames(new_order1), rownames(summary_results2)),]
  #print (summary_results3)
  rownames(summary_results3) <- 1:nrow(summary_results3)
  
  
  summary_results3['phylo'] <- rownames(summary_results3)
  summary_results3_dups <- summary_results3[order(summary_results3$group),] 
  rownames(summary_results3_dups) <- 1:nrow(summary_results3_dups)

  return(summary_results3_dups)
}

## plot regression
ggplotRegression <- function (fit, name) {
  ## https://sejohnston.com/2012/08/09/a-quick-and-easy-function-to-plot-lm-results-in-r/
  require(ggplot2)
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5))) + 
    theme(plot.title = element_text(size=10)) +
    theme_classic() +
    labs(y = "Number of duplicates") + labs(x = name)
  
}

## plot lines of counts
ggplot_line <- function(df, sep, xlabel) {
  df['id'] <- rownames(df)
  df = df %>% select(., id, group, dups, transpo) %>% 
    pivot_longer(., -id,names_to = "Line", values_to = "Value")
  
 ggline(df, "id", y="Value", color="Line", 
         palette = "jco", add="jitter", point.size = 1) + 
    scale_x_discrete(breaks=seq(0,1000, as.integer(sep))) + 
    labs(y = "Counts") + labs(x = xlabel)
}

## linea model representation
lm_plotter_transpo <- function(data, name) {
  fit1 <- lm(dups ~ transpo, data)
  return(ggplotRegression(fit1, name))
}
lm_plotter_hypo <- function(data, name) {
  fit1 <- lm(dups ~ hypo , data)
  return(ggplotRegression(fit1, name))
}
lm_plotter_phages <- function(data, name) {
  fit1 <- lm(dups ~ phages , data)
  return(ggplotRegression(fit1, name))
}
lm_plotter_DUF <- function(data, name) {
  fit1 <- lm(dups ~ DUF, data)
  return(ggplotRegression(fit1, name))
}

# check variables
check_variables <- function(df, name){
  
  ## create binary variable for phages
  df$phages_bin <- ifelse(df$phages > 4, ">4", "<=4")
  
  ## phages
  my_comparisons3 <- combn(x = as.character(unique(df$phages_bin)), m = 2)
  my_comparisons3_1 <- split(my_comparisons3, rep(1:ncol(my_comparisons3), each = nrow(my_comparisons3)))
  ggboxplot_phages <- ggboxplot(df, x="phages_bin", y="dups",  fill = "phages_bin",
                                palette = "jco", add="jitter", point.size = 2, 
                                title = "Phages inserted") + 
    stat_compare_means(method="t.test", comparisons = my_comparisons3_1) +
    scale_x_discrete(limits=c("<=4",">4")) + labs(y = "Number of duplicates") + labs(x = "Phages inserted")
  
  
  ## resistance
  my_comparisons4 <- combn(x = as.character(unique(df$resistance)), m = 2)
  my_comparisons4_1 <- split(my_comparisons4, rep(1:ncol(my_comparisons4), each = nrow(my_comparisons4)))
  ggboxplot_res <- ggboxplot(df, x="resistance", y="dups", fill = "resistance", 
                             palette = "jco", add="jitter", point.size = 2, 
                             title = "Resistance")+
    stat_compare_means(method = "t.test", comparisons = my_comparisons4_1) + 
    labs(y = "Number of duplicates") + labs(x = "Resistance")
  
  ## create pdf
  transpo <- lm_plotter_transpo(df, "Transposases")
  hypo <- lm_plotter_hypo(df, "Hypothetical")
  DUF <- lm_plotter_DUF(df, "Unknown function")
  phages <- ggboxplot_phages
  res <- ggboxplot_res
  
  arrange_plots <- ggarrange(transpo, hypo, DUF, 
            labels = c("A", "B", "C"),
            ncol = 1, nrow = 3)
  arrange_plots2 <- ggarrange(phages, res, 
                             labels = c("D", "E"),
                             ncol = 1, nrow = 2)

    ## create pdf
  pdf(name, width = 8.27, height = 11.69)
  print(arrange_plots)
  print(arrange_plots2)
  dev.off()
  
  #pdf(name)
  #par(mfrow=c(3,2))
  #print(lm_plotter_transpo(df, "Transposases"))
  #print(lm_plotter_hypo(df, "Hypothetical"))
  #print(lm_plotter_DUF(df, "Unknown function"))
  #print(ggboxplot_phages)
  #print(ggboxplot_res)
  #dev.off()
}
######################

######################
## Efaecium
######################
#Efaecium_stats_files <- "./SupData_paper/Efaecium"
Efaecium_stats_files <- "C:/Users/josed/Desktop/Developer/BacterialDuplicates/dataFromHercules/summary_results_final/"

Efaecium_stats <- read.csv(paste0(Efaecium_stats_files, "/20200423_Efaecium_summary_results.csv"), 
                           header = FALSE, col.names = c("strain", "group", "dups", "all_seqs", "transpo", "hypo", "DUF"))
order_Efaecium_nwk <- paste0(Efaecium_stats_files,"/20200423_Efaecium_order-tree.txt")
Efaecium_phispy <- paste0(Efaecium_stats_files,"/20200423_Efaecium_phispy_results.csv")
Efaecium_resistance <- paste0(Efaecium_stats_files,"/20200423_Efaecium_resistance.csv")

Efaecium_summary_results_ordered <- order_data(Efaecium_stats, order_Efaecium_nwk, Efaecium_phispy, Efaecium_resistance)
head(Efaecium_summary_results_ordered)
dim(Efaecium_summary_results_ordered)

## create new resistance column
Efaecium_summary_results_ordered["resistance_original"] <- Efaecium_summary_results_ordered$resistance
Efaecium_summary_results_ordered$resistance <- ifelse(Efaecium_summary_results_ordered$resistance_original =="na", "None", "VRE")

## check variables
#out_pdf_name_Efaecium <- "Efaecium_variables.pdf"
out_pdf_name_Efaecium <- "C:/Users/josed/Desktop/Developer/BacterialDuplicates/process_data/Efaecium_variables.pdf"
#check_variables(Efaecium_summary_results_ordered, out_pdf_name_Efaecium)

ggplot_line(Efaecium_summary_results_ordered, 10, "Strains")

######################

######################
## Efaecalis
######################
#Efaecalis_stats_files <- "./SupData_paper/Efaecalis"
Efaecalis_stats_files <- "C:/Users/josed/Desktop/Developer/BacterialDuplicates/dataFromHercules/summary_results_final/"
Efaecalis_stats <- read.csv(paste0(Efaecalis_stats_files, "/20200423_Efaecalis_summary_results.csv"),
  header = FALSE, col.names = c("strain", "group", "dups", "all_seqs", "transpo", "hypo", "DUF"))
order_Efaecalis_nwk <- paste0(Efaecalis_stats_files,"/20200423_Efaecalis_order-tree.txt")
Efaecalis_phispy <- paste0(Efaecalis_stats_files,"/20200423_Efaecalis_phispy_results.csv")
Efaecalis_resistance <- paste0(Efaecalis_stats_files,"/20200423_Efaecalis_resistance.csv")

## order results
Efaecalis_summary_results_ordered <- order_data(Efaecalis_stats, order_Efaecalis_nwk, Efaecalis_phispy, Efaecalis_resistance)

## create new resistance column
Efaecalis_summary_results_ordered["resistance_original"] <- Efaecalis_summary_results_ordered$resistance
Efaecalis_summary_results_ordered$resistance <- ifelse(Efaecalis_summary_results_ordered$resistance_original =="na", "None", "VRE")

## check variables
#out_pdf_name_Efaecalis <- "Efaecalis_variables.pdf"
out_pdf_name_Efaecalis <- "C:/Users/josed/Desktop/Developer/BacterialDuplicates/process_data/Efaecalis_variables.pdf"
#check_variables(Efaecalis_summary_results_ordered, )
ggplot_line(Efaecalis_summary_results_ordered, 10, "Strains")

######################

######################
## Saureus
######################
#Saureus_stats_files <- "./SupData_paper/Saureus"

Saureus_stats_files <- "C:/Users/josed/Desktop/Developer/BacterialDuplicates/dataFromHercules/summary_results_final/"
Saureus_stats <- read.csv(paste0(Saureus_stats_files, "20200423_Saureus_summary_results.csv"),
                          header = FALSE, col.names = c("strain", "group", "dups", "all_seqs", "transpo", "hypo", "DUF"))
order_Saureus_nwk <- paste0(Saureus_stats_files, "/20200423_Saureus_order-tree.txt")
Saureus_phispy <- paste0(Saureus_stats_files, "/20200423_Saureus_phispy_results.csv")
Saureus_resistance <- paste0(Saureus_stats_files, "/20200423_Saureus_resistance.csv")

Saureus_summary_results_ordered <- order_data(Saureus_stats, order_Saureus_nwk, Saureus_phispy, Saureus_resistance)
Saureus_summary_results_ordered

## check variables
#out_pdf_name_Saureus <- "Saureus_variables.pdf"
out_pdf_name_Saureus <- "C:/Users/josed/Desktop/Developer/BacterialDuplicates/process_data/Saureus_variables_phages2.pdf"
#check_variables(Saureus_summary_results_ordered, out_pdf_name_Saureus)
ggplot_line(Saureus_summary_results_ordered, 25, "Strains")
######################

