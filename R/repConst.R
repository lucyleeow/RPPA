#' Replicate consistency
#' 
#' Determine replicate consistency. by calculating the mean, standard deviation
#' and coefficient of variation (CV).
#' 
#' 
#' 
#' @param tidydf Tidy dataframe of RFI data.
#' @param reps Vector of replicate names as string, to calculate consistency of.
#' @param title Title of the faceted barplot, as string.
#' @param filename Filename of the replicate consistency data (include .tsv at 
#'     the end).
#' 
#' 
#' @importFrom assertthat assert_that
#' @importFrom magrittr %>%
#' @import ggplot2
#' 

#' @describeIn plot_repConst
#' @export
plot_repConst <- function(tidydf, reps, title) {
  
  # check inputs
  assert_that(sum(c("Condition", "Antibody.Name", "RFI") %in% colnames(tidydf))
              == 3, msg = "Check there are 'Condition', 'Antibody.Name' and 
                           'RFI' columns in 'tidydf'")
  
  assert_that(is.character(reps), msg = "Check that 'reps' is a string vector")
  
  assert_that(is.character(title), length(title) == 1,
              msg = "Check that 'title' is a single string")
  
  
  
  tidydf %>%
    dplyr::filter(Condition %in% reps) %>%
    dplyr::group_by(Condition, Antibody.Name) %>%
    dplyr::summarise(mean = mean(RFI), sd = sd(RFI)) %>%
    ggplot(aes(y = mean, x = Condition)) + 
    geom_bar(stat = "identity", position = "dodge") + 
    geom_errorbar(aes(ymin=mean-sd, ymax = mean+sd, width = 0.5)) + 
    facet_wrap(Antibody.Name~.) + 
    coord_flip() +
    labs(title = title, y = "Mean RFI for each group (error bar=SD)")
  
}

#' @describeIn plot_repConst
#' @export 
df_repConst <- function(tidydf, reps, filename) {
  
  # check inputs
  assert_that(sum(c("Condition", "Antibody.Name", "RFI") %in% colnames(tidydf))
              == 3, msg = "Check there are 'Condition', 'Antibody.Name' and 
                           'RFI' columns in 'tidydf'")
  
  assert_that(is.character(reps), msg = "Check that 'reps' is a string vector")
  
  
  const_df <- tidydf %>%
    dplyr::filter(Condition %in% reps) %>%
    dplyr::group_by(Condition, Antibody.Name) %>%
    dplyr::summarise(Mean = mean(RFI), SD = sd(RFI), CV = SD/Mean) %>%
    
    
  write.table(const_df, file = filename, sep = "\t", quote = FALSE,
              row.names = FALSE)
  
}






