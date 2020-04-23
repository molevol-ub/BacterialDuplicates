library(ggplot2)
ggplotRegression <- function (fit, name) {
  ## https://sejohnston.com/2012/08/09/a-quick-and-easy-function-to-plot-lm-results-in-r/
  require(ggplot2)
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste(name, ": Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)))
}

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
generate_bin <- function(df) {
  df$dups_CARD_bin <- ifelse(df$dups_CARD > 1, ">1", "0")
  df$dups_VFDB_bin <- ifelse(df$dups_VFDB > 1, ">1", "0")
  df$phages_bin <- ifelse(df$phages > 2, ">2", "<2")
  return(df)
}

check_variables <- function(df, name){
  
  df <- generate_bin(df)
  
  ## phages
  my_comparisons3 <- combn(x = as.character(unique(df$phages_bin)), m = 2)
  my_comparisons3_1 <- split(my_comparisons3, rep(1:ncol(my_comparisons3), each = nrow(my_comparisons3)))
  ggboxplot_phages <- ggboxplot(df, x="phages_bin", y="dups",  fill = "phages_bin",
                                palette = "jco", add="jitter", point.size = 2, 
                                title = "Phages inserted") + stat_compare_means(method="t.test", comparisons = my_comparisons3_1)
  
  
  ## resistance
  my_comparisons4 <- combn(x = as.character(unique(df$Resistance)), m = 2)
  my_comparisons4_1 <- split(my_comparisons4, rep(1:ncol(my_comparisons4), each = nrow(my_comparisons4)))
  ggboxplot_res <- ggboxplot(df, x="Resistance", y="dups", fill = "Resistance", 
                             palette = "jco", add="jitter", point.size = 2, 
                             title = "Resistance")+
    stat_compare_means(method = "t.test", comparisons = my_comparisons4_1)
  
  
  pdf(name)
  print(lm_plotter_transpo(df, "Transposases"))
  print(lm_plotter_hypo(df, "Hypothetical"))
  print(ggboxplot_phages)
  print(ggboxplot_res)
  dev.off()
  
}

### load data
Efaecium_stats <- read.csv("C:/Users/josed/Desktop/Developer/BacterialDuplicates/process_data/Efaecium_commas.csv")
Efaecalis_stats <- read.csv("C:/Users/josed/Desktop/Developer/BacterialDuplicates/process_data/Efaecalis_commas2.csv")
Saureus_stats <- read.csv("C:/Users/josed/Desktop/Developer/BacterialDuplicates/process_data/Saures_comma.csv")

check_variables(Saureus_stats, "C:/Users/josed/Desktop/Developer/BacterialDuplicates/Saureus_variables.pdf")
check_variables(Efaecalis_stats, "C:/Users/josed/Desktop/Developer/BacterialDuplicates/Efaecalis_variables.pdf")
check_variables(Efaecium_stats, "C:/Users/josed/Desktop/Developer/BacterialDuplicates/Efaecium_variables.pdf")

