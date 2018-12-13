#' Replicate consistency
#' 
#' Determine replicate consistency by calculating the mean, standard deviation
#' and coefficient of variation (CV).
#' 
#' 
#' 
#' 
#' @param tidydf
#' @param reps Vector of replicate names as string, to calculate consistency of.
#' @param title Title of the faceted barplot, as string.
#' 
#' 
#' 
#' @importFrom assertthat assert_that
#' @importFrom magrittr %>%
#' @import ggplot2
#' 
#' 



plot_repConst <- function(tidydf, reps, title) {
  
  # check inputs
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


df_repConst <- function(tidydf, reps) {
  
  # check inputs
  assert_that(is.character(reps), msg = "Check that 'reps' is a string vector")
  
  
  tidydf %>%
    dplyr::filter(Condition %in% reps) %>%
    dplyr::group_by(Condition, Antibody.Name) %>%
    dplyr::summarise(Mean = mean(RFI), SD = sd(RFI), CV = SD/Mean) %>%
    
  
  
}






